import os 
os.environ['OPENBLAS_NUM_THREADS'] = "1"
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np
import sys


class HamitonianMC:
    def __init__(self,UserDefinedModel):
        self.invert_Mass = None 
        self.model = UserDefinedModel
        self.boundaries = np.zeros((2,2)) 
        self.dt = None 
        self.Lrange = [0,0]

        self.seed = None
        self.myrank = None
        self.save_folder = None

        self.cache = {}
        self.ii = 0

    def _set_initial_model(self):
        n = self.boundaries.shape[0]
        xcur = np.zeros((n))
        for i in range(n):
            bdl = self.boundaries[i,0]
            bdr = self.boundaries[i,1]
            xcur[i] = bdl + (bdr - bdl) * np.random.rand()

        # make the vel array ordered
        
        
        tmp = xcur[:int(len(xcur)/2)]
        idx = np.argsort(tmp)
        tmp_order = tmp[idx]
        xcur[:int(len(xcur)/2)] = tmp_order
        
        tmp = xcur[int(len(xcur)/2):]
        tmp_order = tmp[idx]
        xcur[int(len(xcur)/2):] = tmp_order        
        
        #print("initial model :",xcur)
        self._save_init_model(xcur)
        return xcur
        
    def _check_init_is_in_boundary(self,xcur):
        for i in range(len(xcur)-1):
            if xcur[i] < self.boundaries[i,0] or xcur[i] > self.boundaries[i,1]:
                return False
        return True
            

    def _kinetic(self,p):
        """
        Kinetic energy
        """ 
        K = np.dot(self.invert_Mass @ p,p) * 0.5 

        return K 
    
    def _forward(self,x):
        return self.model.forward(x)
    
    def _misfit_and_grad(self,x):
        """
        compute misfit function and corresponding gradient
        """
        misfit,grad,dsyn,flag = self.model.misfit_and_grad(x)

        return misfit,grad,dsyn,flag     

    def _leapfrog(self,xcur,dt,L,i):
        """
        leap frog scheme
        """
        n = len(xcur)
        #np.random.seed(self.seed + self.myrank + i)
        pcur = np.random.randn(n) * 0.5

        # initialize xnew and pnew
        pnew = pcur * 1.0 
        xnew = xcur * 1.0

        # compute current Hamiltonian
        K = self._kinetic(pnew)
        U,grad,dsyn,flag = self._misfit_and_grad(xnew)
        if flag == False or np.sum(np.isnan(dsyn)) > 0:
            return xcur,np.inf,self.model.dobs,False
        Hcur = K + U

        # save current potential and synthetics
        dsyn_new = dsyn.copy()
        Unew = U

        # update
        pnew -= dt * grad * 0.5 
        for i in range(L):
            pp = pnew.copy()
            xnew += dt * pnew # update xnew
            xx = xnew.copy()

            # check boundaries
            xtmp = xnew.copy()
            ptmp = pnew.copy()
            high = self.boundaries[:,1]
            low = self.boundaries[:,0]
            idx1 = xtmp > high
            idx2 = xtmp < low 
            while np.sum(np.logical_or(idx1,idx2))>0:
                xtmp[idx1] = 2 * high[idx1] - xtmp[idx1]
                ptmp[idx1] = -ptmp[idx1]
                xtmp[idx2] = 2 * low[idx2] - xtmp[idx2]
                ptmp[idx2] = -ptmp[idx2]
                idx1 = xtmp > high 
                idx2 = xtmp < low 
                
            # end while 
            pnew = ptmp.copy()
            xnew = xtmp.copy()

            # update pnew
            #if self.ii / 5000 * 100 > 4.3:
            #    print(self.ii,xnew,xx,pp)
            if np.sum(np.isnan(xnew)) > 0:
                return xcur,np.inf,self.model.dobs,False
            Unew,grad,dsyn_new,flag = self._misfit_and_grad(xnew)
            if np.sum(np.isnan(grad)) > 0:
                # the grad might be NaN. check it.
                return xcur,np.inf,self.model.dobs,False            
            if flag == False or np.sum(np.isnan(dsyn_new)) > 0:
                return xcur,np.inf,self.model.dobs,False
            if i < L-1:
                pnew -= dt * grad
            #if self.ii / 5000 * 100 > 30.01:
            #    print(grad,dsyn_new)
        # end for loop
        pnew = -pnew 

        # update Hamiltonian
        Knew = self._kinetic(pnew)
        Hnew = Knew + Unew


        # accept or not
        AcceptFlag = False
        u = np.random.rand()
        if u < np.exp(-(Hnew - Hcur)):
            xcur = xnew
            U = Unew 
            dsyn = dsyn_new
            AcceptFlag = True
        
        return xcur,U,dsyn,AcceptFlag
    
    def _save_results(self,x,dsyn,idx):
        # save model
        f = open(self.save_folder + "/" + "model" + str(idx) + ".dat","w")
        for i in range(x.shape[0]):
            f.write("%f\n"%(x[i]))
        f.close()

        # save synthetics
        f = open(self.save_folder + "/" + "synthetic" + str(idx) + ".dat","w")
        n = dsyn.shape[0]
        for i in range(n):
            f.write("%f %f\n"%(self.model.dobs[i],dsyn[i]))
        f.close()

    def _save_init_model(self,x):
        # save model
        if not os.path.exists(self.save_folder):
            os.mkdir(self.save_folder)
        f = open(self.save_folder + "/" + "initmodel" + ".dat","w")
        for i in range(x.shape[0]):
            f.write("%f\n"%(x[i]))
        f.close()

               

    def sample(self,nsamples,ndraws,**kwargs):
        # set random seed
        np.random.seed(self.seed)

        # init model
        x = self._set_initial_model()
        
        # check the init model 
        while not self._check_init_is_in_boundary(x):
            x = self._set_initial_model()

        # sample posterior distributions
        misfit = np.zeros((nsamples))
        x_cache = np.zeros((nsamples,len(x)))
        syndata = np.zeros((nsamples,self.model.dobs.shape[0]))
        ncount = 0
        i = 0
        while i < ndraws + nsamples:
            #np.random.seed(self.seed + self.myrank + i)
            L = np.random.randint(self.Lrange[0],self.Lrange[1]+1)
            x,U,dsyn,AcceptFlag = self._leapfrog(x,self.dt,L,i)
            if AcceptFlag: # accept this new sample
                if i>=ndraws:
                    #self._save_results(x,dsyn,i-ndraws)
                    misfit[i-ndraws] = U
                    x_cache[i-ndraws,:] = x.copy()
                    syndata[i-ndraws,:] = dsyn
                i += 1   
                self.ii += 1         
            ncount += 1
            if i > -1:
                msg = "chain {}: {:.2%}, misfit={} -- accept ratio {:.2%}".  \
                    format(self.myrank,i/(ndraws + nsamples),U,i/ncount)
                    
                print(msg)
                #print(x)
                sys.stdout.flush()
        # end i-loop

        # compute average of the best n models
        nbests = 10
        idx = np.argsort(misfit)
        xmean = np.mean(x_cache[idx[:nbests],:],axis=0)
        _,_,dsyn,_ = self._misfit_and_grad(xmean)

        # save synthetics and each sample
        if not os.path.exists(self.save_folder):
            os.mkdir(self.save_folder)
        self._save_results(xmean,dsyn,"mean")
        for i in range(nsamples):
            self._save_results(x_cache[i,:],syndata[i,:],i)

        return misfit

def HMCSample(model,nsamples,ndraws,boundaries,
            delta,Lrange,seed,myrank=0,
            save_folder = "mychain"):
    """
    HMC sampling function
    """
    chain = HamitonianMC(model)
    chain.myrank = myrank
    chain.save_folder = save_folder + str(myrank)
    chain.seed = seed + myrank
    nt = boundaries.shape[0]
    chain.boundaries = boundaries
    chain.Lrange = Lrange
    chain.dt = delta
    chain.invert_Mass = np.eye(nt)

    misfit = chain.sample(nsamples,ndraws)

    return misfit
