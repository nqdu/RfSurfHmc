from model.model_rf import ReceiverFunc
from model.model_surf import SurfWD
import numpy as np 

class Joint_RF_SWD:
    def __init__(self,sigma1,sigma2,rfmodel:ReceiverFunc,swdmodel:SurfWD) -> None:
        self.sigma1 = sigma1
        self.sigma2 = sigma2 
        self.rfmodel = rfmodel 
        self.swdmodel = swdmodel

        self.ndata = rfmodel.nt + swdmodel.nt
    
    def set_obsdata(self,rfobs:np.ndarray, swdobs:np.ndarray):
        self.dobs = np.zeros((self.ndata))
        self.rfobs = rfobs 
        self.swdobs = swdobs



        self.rfmodel.set_obsdata(rfobs)
        self.swdmodel.set_obsdata(swdobs)

        self.dobs[:self.rfmodel.nt] = self.rfobs * 1. 
        self.dobs[self.rfmodel.nt:] = self.swdobs * 1.

    def forward(self,x:np.ndarray):
        """
        compute surface wave dispersion for a given model x   
        Parameters:
        ==================================================
        x   : np.ndarray, shape(n,4)
            input model, [thk,vp,vs,rho]
        
        Returns:
        ==================================================
        d   : np.ndarray,shape(self.nt)
            output dispersion dat
        """
        n1 = self.rfmodel.nt
        n2 = self.swdmodel.nt
        drf = np.zeros((n1))
        dswd = np.zeros((n2))

        

    
        drf = self.rfmodel.forward(x)
        #print(drf)
        dswd,flag = self.swdmodel.forward(x)

        return drf,dswd,flag 
    
    def misfit(self,x:np.ndarray):
        drf,dswd,flag = self.forward(x)
        if flag:
            n1 = drf.size 
            n2 = dswd.size 
            wt = (self.sigma1 / self.sigma2)**2 * n1 / n2 
            
            
            misfit = 0.5 * np.sum((drf - self.rfobs)**2) + \
                    0.5 * np.sum((dswd - self.swdobs)**2) * wt               
                    
            return misfit,True
        else:
            return 0.0,flag
    
    def misfit_and_grad(self,x:np.ndarray):
        # compute misfit and grad for each dataset
        
        misfitr,gradr,dr = self.rfmodel.misfit_and_grad(x)
        misfits,grads,ds,flag = self.swdmodel.misfit_and_grad(x)

        # check flag
        if flag == False:
            return 0.0,np.zeros((gradr.shape)),self.dobs,flag 
        
        # compute weighted misfit
        n1 = dr.size 
        n2 = ds.size 
        wt = (self.sigma1 / self.sigma2)**2 * n1 / n2 
        misft = misfitr + wt * misfits
        grad = gradr + wt * grads 
        dobs = np.zeros((n1 + n2))
        dobs[:n1] = dr 
        dobs[n1:] = ds 

        return misft,grad,dobs,True
