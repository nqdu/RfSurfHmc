import os 
os.environ['OPENBLAS_NUM_THREADS'] = "1"
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np
import sys
import h5py

"""
@class HMC + Dual Averaging
"""
class HMCDualAveraging:
    def __init__(self,UserDefinedModel,
                boundaries,dt:float,L0:int,nbest_model:int,
                target_ratio:float,seed:int,
                nsamples:int,ndraws:int,myrank=0,
                name = "mychain",outdir = "./"):
        """
        HMC parameters

        Parameters
        ----------
        UserDefinedModel: class
            user defined model, model.misfit_and_grad is required
        boundaries: np.ndarray, shape(n,2)
            search boundaries, (low,high)
        dt: float
            dt in Hamiltonian mechnics
        L0 : int
            initial no. of steps
        target_ratio: float
            target sampling ratio
        seed: int
            random seed
        nsamples: int
            no. of samples required
        ndraws: int
            no. of warmup samples
        nbest_model: int
            how many best models to average, to get best estimation
        myrank: int
            current rank
        name: str
            output name, output as {name}.{myrank}.h5
        """

        nt = boundaries.shape[0]
        self.invert_Mass = np.eye(nt)
        self.model = UserDefinedModel
        self.boundaries = boundaries 
        self.dt = dt
        self.L = L0
        self.nbest_model = nbest_model
        self.nsamples = nsamples
        self.ndraws = ndraws

        # check ndraws > 0.1 * nsamples
        if ndraws < 0.1 * nsamples:
            print("in dual averaging, ndraws should > nsamples * 0.1")
            print(f"ndraws = {ndraws}, nsamples = {nsamples}")
            exit(1)

        self.seed = seed + myrank
        self.myrank = myrank
        self.name = name

        self.cache = {}
        self.ii = 0

        # parameters for dual averaging
        self.delta = target_ratio
        self.dtbar = dt 
        self._h0 = 0.
        self._gamma = 0.05 
        self._t0 = 10.
        self._kappa = 0.75
        self._lambda = L0 * self.dt

        # open hdf5 to save results
        os.makedirs(outdir,exist_ok=True)
        self.fio = h5py.File(f"{outdir}/{name}.{self.myrank}.h5","w")

        # init global random seed
        np.random.seed(self.seed)
    
    @classmethod
    def init(self,UserDefinedModel,boundaries,rank,**kargs):
        return HMCDualAveraging(
            UserDefinedModel,boundaries,
            kargs['dt'],kargs['L0'],
            kargs['nbest'],
            kargs['target_ratio'],
            kargs['seed'],
            kargs['nsamples'],
            kargs['ndraws'],
            rank,kargs['name'],
            kargs['OUTPUT_DIR']
        )

    def sort_initial_model(self,xcur:np.ndarray):
        # make the vel array ordered        
        tmp = xcur[:int(len(xcur)/2)]
        idx = np.argsort(tmp)
        tmp_order = tmp[idx]
        xcur[:int(len(xcur)/2)] = tmp_order
        
        tmp = xcur[int(len(xcur)/2):]
        tmp_order = tmp[idx]
        xcur[int(len(xcur)/2):] = tmp_order 

        return xcur

    def set_initial_model(self):
        n = self.boundaries.shape[0]
        xcur = np.zeros((n))
        for i in range(n):
            bdl = self.boundaries[i,0]
            bdr = self.boundaries[i,1]
            xcur[i] = bdl + (bdr - bdl) * np.random.rand()

        # sort velocity profile
        xcur = self.sort_initial_model(xcur)
        
        #print("initial model :",xcur)
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
    
    def _mirror(self,xcur,pcur):
        # check boundaries
        xtmp = xcur.copy()
        ptmp = pcur.copy()
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
        
        return xtmp,ptmp
    
    def _find_initial_dt(self,dt0,x):
        xcur = x.copy()
        n = len(xcur)
        dt = dt0
        pcur = np.random.randn(n) * 0.5

        # compute gradient and potential energy
        U,grad,_,flag = self._misfit_and_grad(xcur)
        K = self._kinetic(pcur)
        Hcur = U + K

        # leap frog to change dt
        a = 0.
        pcur -= 0.5 * dt * grad 
        for it in range(20):
            xcur += dt * pcur

            # mirror
            xcur,pcur = self._mirror(xcur,pcur) 

            # compute grad for current x
            U,grad,_,flag = self._misfit_and_grad(xcur)

            if not flag:
                print('error in chain %d!\n' %(self.myrank))
                exit(1)

            # update pcur
            pcur -= 0.5 * dt * grad

            # compute new H 
            K = self._kinetic(pcur)
            Hnew = U + K 

            # factor to increase/decrease
            ediff = -(Hnew - Hcur)
            if it == 0:
                a = 2 * (ediff > np.log(0.5)) - 1
            
            if ediff < np.log(0.5):
                break
            else:
                # update P for another half step
                pcur -= 0.5 * dt * grad
                Hcur = Hnew * 1.
                dt = dt * 2**a
        if it == 20:
            print(f"dt estimation for chain {self.myrank} is not reliable!")
        print(f"chain {self.myrank}: change dt from {dt0} to {dt}")

        return dt 

    def _leapfrog(self,xcur,dt,L):
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
            return xcur,np.inf,self.model.dobs,0.
        Hcur = K + U

        # save current potential and synthetics
        dsyn_new = dsyn.copy()
        Unew = U

        # leap frog scheme
        pnew -= 0.5 * dt * grad 
        for it in range(L):
            # update x
            xnew += dt * pnew

            # check boundaries
            xnew,pnew = self._mirror(xnew,pnew)
            if np.sum(np.isnan(xnew)) > 0:
                return xcur,np.inf,self.model.dobs,0

            # comptue gradient
            Unew,grad,dsyn_new,flag = self._misfit_and_grad(xnew)
            if np.sum(np.isnan(grad)) > 0:
                # the grad might be NaN. check it.
                return xcur,np.inf,self.model.dobs,0           
            if flag == False or np.sum(np.isnan(dsyn_new)) > 0:
                return xcur,np.inf,self.model.dobs,0

            # update pnew
            if it < L - 1:
                pnew -= dt * grad
            else :
                pnew -= dt * grad * 0.5
        # end for loop

        # update Hamiltonian
        pnew = -pnew 
        Knew = self._kinetic(pnew)
        Hnew = Knew + Unew

        # ratio
        alpha = min(1.,np.exp(-(Hnew - Hcur)))
        return xnew,Unew,dsyn_new,alpha
    
    def sample(self):
        # init model
        x = self.set_initial_model()
        
        # check the init model 
        while not self._check_init_is_in_boundary(x):
            x = self.set_initial_model()
        self._save_init(x)

        # sample posterior distributions
        nsamples = self.nsamples
        ndraws = self.ndraws
        misfit = np.zeros((nsamples))
        x_cache = np.zeros((nsamples,len(x)))
        syndata = np.zeros((nsamples,self.model.dobs.shape[0]))

        # get initial dt 
        dt = self._find_initial_dt(self.dt,x)
        dtbar = dt * 1.
        h0 = self._h0
        mu = np.log(10 * self.dt)

        # sample begin
        ncount = 0
        i = 0
        while i < ndraws + nsamples:
            #np.random.seed(self.seed + self.myrank + i)
            L = max(1,int(self._lambda / dt))
            x1,U,dsyn,alpha = self._leapfrog(x,dt,L)

            # check if accept this model
            u = np.random.rand()
            AcceptFlag = False
            if u < alpha:
                x = x1.copy()
                AcceptFlag = True

            if AcceptFlag: # accept this new sample
                if i>=ndraws:
                    #self.save_results(x,dsyn,i-ndraws)
                    misfit[i-ndraws] = U
                    x_cache[i-ndraws,:] = x.copy()
                    syndata[i-ndraws,:] = dsyn

                # update i
                i += 1   
                self.ii += 1

            # dual averaging
            if ncount < self.ndraws:
                # update h
                m = ncount + 1
                fac = 1. / (m + self._t0)
                h0 = (1 - fac) * h0 + fac * (self.delta - alpha)
                
                # update dt
                logdt = mu - np.sqrt(m) / self._gamma * h0
                dt = np.exp(logdt)

                # update eps_cur
                fac = m**(-self._kappa)
                logdtbar = fac * logdt + (1 - fac) * np.log(dtbar)
                dtbar = np.exp(logdtbar)
            else:
                # set dt_prev
                dt = dtbar * 1.   

            # update counter     
            ncount += 1
            if i > -1 and (i % 50 == 0 or i == nsamples - 1):
                msg = "chain {}: {:.2%}, dt = {:.3},  misfit={:.3} -- accept ratio {:.2%}".  \
                    format(self.myrank,i/(ndraws + nsamples),dt,U,i/ncount)
                    
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
        self.save_results(xmean,dsyn,"mean")
        for i in range(nsamples):
            self.save_results(x_cache[i,:],syndata[i,:],i)

        return misfit
    
    def save_results(self,x,dsyn,idx):

        # create group
        gname = f"{idx}"
        self.fio.create_group(f"{gname}")

        # save model
        self.fio.create_dataset(
            f"{gname}/model",dtype='f8',shape=x.shape)
        self.fio[f"{gname}/model"][:] = x

        # save synthetics
        self.fio.create_dataset(
            f"{gname}/syn",dtype='f8',shape=dsyn.shape
        )
        
        self.fio[f"{gname}/syn"][:] = dsyn[:]

    def _save_init(self,x):
        # save model
        self.fio.create_dataset("initmodel",dtype='f8',shape=x.shape)
        self.fio['initmodel'][:] = x

        # save obs
        self.fio.create_dataset("obs",dtype='f8',shape=self.model.dobs.shape)
        self.fio["obs"][:] = self.model.dobs[:]

    def close(self):
        # close h5 file
        self.fio.close()