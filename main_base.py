from matplotlib.colors import BoundaryNorm
import numpy as np
from pyhmc.hmc import HamitonianMC
from mpi4py import MPI
from model.model_rf import ReceiverFunc
from model.model_surf import SurfWD
from model.model_rf_swd_vs_thk import Joint_RF_SWD
import matplotlib.pyplot as plt 
import time
import yaml
import os


def main():
    # mpi information
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()

    # load yaml files
    with open("param.yaml","r") as f:
        param = yaml.safe_load(f)

    # surf parameters
    sparam = param['swd']
    model_swd = SurfWD.init(**sparam)

    # receiver function parameters
    rparam = param['rf']
    model_rf = ReceiverFunc.init(**rparam)
    
    # set model parameters
    thk = np.asarray(param['true_model']['thk'])
    vs = np.asarray(param['true_model']['vs'])
    model_swd.set_thk(thk)
    model_rf.set_thk(thk)  

    # joint inversion model
    sigma1 = 1.0; sigma2 = 1.0
    model = Joint_RF_SWD(sigma1,sigma2,model_rf,model_swd)

    # get output dir
    outdir = param['hmc']['OUTPUT_DIR']
    os.makedirs(outdir,exist_ok=True)

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
        np.save(f"{outdir}/real_syn.npy",dobs)

    # bcast to other procs
    dobs = comm.bcast(dobs)
    x = comm.bcast(x)
    nt = model.rfmodel.nt
    model.set_obsdata(dobs[:nt],dobs[nt:])

    # set the search range of inversion
    n = len(x)
    boundaries = np.ones((n,2))
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

    # initialize HMC sampler
    hparam = param['hmc']
    chain = HamitonianMC.init(model,boundaries,rank,**hparam)
    nsamples = chain.nsamples
    tmp = chain.sample()

    # gather misfits
    if rank == 0:
        misfit = np.zeros((ncores,nsamples))
    else:
        misfit = None
    comm.Gather(tmp,misfit,root=0)

    if rank == 0:
        np.save(f"{outdir}/misfit.npy",misfit)

if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print("time elapse: {}".format(toc-tic))
