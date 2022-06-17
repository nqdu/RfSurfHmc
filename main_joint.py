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
    sigma1 = 1.0; sigma2 = 1.0
    # set thickness
    
    

    thk_real = np.array([20,20,20,0]) 
    vs_real = np.array([3.0,3.5,4,4.5])
    
    #merge the second layers
    thk_inv =  np.array([20,20,20,0]) 
    vs_inv = np.array([3.0,3.5,4,4.5])
    
    model_swd.set_thk(thk_real)
    model_rf.set_thk(vs_real)
    model_real = Joint_RF_SWD(sigma1,sigma2,model_rf,model_swd)
    x = np.hstack((vs_real,thk_real))
    drsyn,dssyn,_ = model_real.forward(x)  

    # joint inversion model
    model_swd.set_thk(thk_inv)
    model_rf.set_thk(vs_inv)
    model = Joint_RF_SWD(sigma1,sigma2,model_rf,model_swd)

    # obs data
    dobs = None
    x = None
    if rank == 0:
        # true model
        x = np.hstack((vs_inv,thk_inv))
        dobs = np.zeros((model.ndata))
        
        
        dobs[:model.rfmodel.nt] = drsyn
        dobs[model.rfmodel.nt:] = dssyn  
        np.savetxt("real_data",dobs)

        
    if ncores > 1:
        dobs = comm.bcast(dobs)
        x = comm.bcast(x)
        comm.barrier()
    model.set_obsdata(dobs[:nt],dobs[nt:])

    # parameters
    save_folder = "chain_joint"
    n = len(x)
    boundaries = np.ones((n,2))

    for i in range(len(thk_inv)):
        
        boundaries[i,0] = 3
        boundaries[i,1] = 5
           
        boundaries[i + len(thk_inv),0] = 15
        boundaries[i + len(thk_inv),1] = 25
                        
    boundaries[-1,:] = 0.0,2.0
              


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
