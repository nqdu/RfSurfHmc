from matplotlib.colors import BoundaryNorm
import numpy as np
from model.model_rf import ReceiverFunc
from model.model_surf import SurfWD
from model.model_rf_swd_vs_thk import Joint_RF_SWD
import matplotlib.pyplot as plt 
import time

def main():


    # surf parameters
    tRc = np.linspace(5,40,36)
    tRg = tRc.copy()
    model_swd = SurfWD(tRc=tRc,tRg=tRg)

    # receiver function parameters
    rayp = 0.045; dt = 0.1; 
    nt = 500; gauss = 1.0; time_shift = 5.0
    water = 0.001; rf_type = 'P'; rf_method = "time"
    model_rf = ReceiverFunc(rayp,nt,dt,gauss,time_shift,
                            water,rf_type, rf_method)
    
    
    
    
    # set model parameters
    thk =  np.array([6,6,13,5,10,30,0]) 
    # vs = np.array([3.2,2.8,3.46,3.3,3.9,4.5,4.7])
    vs = np.array([3.2,3.4,3.46,3.7,3.9,4.5,4.7])
    

    model_swd.set_thk(thk)
    model_rf.set_thk(thk)  

    # joint inversion model
    sigma1 = 1.0; sigma2 = 1.0
    model = Joint_RF_SWD(sigma1,sigma2,model_rf,model_swd)

    # obs data
    dobs = None
    x = None

    # true model
    x = np.hstack((vs,thk))
    drsyn,dssyn,_ = model.forward(x)


    t_rf = np.linspace(0,(nt-1)*dt, nt) - time_shift
    
    plt.figure(1,figsize=(14,30))
    plt.subplot(3,1,1)
    plt.plot(t_rf,drsyn)
    plt.title("rf_syn")
    plt.subplot(3,1,2)
    plt.plot(tRc,dssyn[:len(tRc)])
    plt.title("Rc_syn")
    plt.subplot(3,1,3)
    plt.plot(tRg,dssyn[len(tRc):])
    plt.title("Rg_syn")

    plt.savefig("./syn_test.png")


if __name__ == "__main__":
    tic = time.time()
    main()
    toc = time.time()
    print("time elapse: {}".format(toc-tic))
