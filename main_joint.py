from matplotlib.colors import BoundaryNorm
import numpy as np
from pyhmc.hmc import HMCSample
from mpi4py import MPI
from model.model_rf import ReceiverFunc
from model.model_surf import SurfWD
from model.model_rf_swd_vs_thk import Joint_RF_SWD
import matplotlib.pyplot as plt 
import time


seed = 19951128
def main():
    # mpi information
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()

    # surf parameters
    tRc = np.linspace(5,40,36)
    tRg = tRc.copy()
    model_swd = SurfWD(tRc=tRc,tRg=tRg)

    # receiver function parameters
    rayp = 0.045; dt = 0.4; 
    nt = 125; gauss = 1.5; time_shift = 5.0
    water = 0.001; rf_type = 'P'
    model_rf = ReceiverFunc(rayp,nt,dt,gauss,time_shift,
                            water,rf_type)
    
    
    
    
    # set model parameters
    thk =  np.array([6,6,13,5,10,30,0]) 
    vs = np.array([3.2,2.8,3.46,3.3,3.9,4.5,4.7])
    

    model_swd.set_thk(thk)
    model_rf.set_thk(thk)  

    # joint inversion model
    sigma1 = 1.0; sigma2 = 1.0
    model = Joint_RF_SWD(sigma1,sigma2,model_rf,model_swd)

    # obs data
    dobs = None
    x = None
    if rank == 0:
        # true model
        x = np.hstack((vs,thk))
        drsyn,dssyn,_ = model.forward(x)
        dobs = np.zeros((model.ndata))
        
        

        
        dobs[:model.rfmodel.nt] = drsyn
        dobs[model.rfmodel.nt:] = dssyn  
        np.savetxt("./real_syn",dobs)


        
    if ncores > 1:
        dobs = comm.bcast(dobs)
        x = comm.bcast(x)
        comm.barrier()
    model.set_obsdata(dobs[:nt],dobs[nt:])

    # parameters
    save_folder = "chain_joint"
    n = len(x)
    boundaries = np.ones((n,2))




    
    # set the search range of inversion
    for i in range(len(thk)):
        boundaries[i,0] = vs[i] - vs[i]*0.8
        boundaries[i,1] = vs[i] + vs[i]*0.8       
        if boundaries[i,0] < 1.5:
            boundaries[i,0] = 1.5
        if boundaries[i,1] > 5:
            boundaries[i,1] = 5
        
        boundaries[i + len(thk),0] = thk[i] - thk[i]*0.2
        boundaries[i + len(thk),1] = thk[i] + thk[i]*0.2     

    boundaries[-1,:] = 0.0,2.0
              

    # set the search range of inversion
    Lrange = [5,10] 
    delta = 0.05
    nsamples = 200
    ndraws = 0

    
    if ncores == 1:
        misfit = HMCSample(model,nsamples,ndraws,boundaries,
                        delta,Lrange,seed,myrank=rank,
                        save_folder=save_folder)
    else:
        misfit = None
        tmp = HMCSample(model,nsamples,ndraws,boundaries,
                            delta,Lrange,seed+rank,myrank=rank,
                            save_folder=save_folder)
        if rank == 0:
            misfit = np.zeros((ncores,nsamples))
        comm.Gather(tmp,misfit,root=0)
    
    if rank == 0:
        np.savetxt("misfit.dat",misfit)

if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print("time elapse: {}".format(toc-tic))
