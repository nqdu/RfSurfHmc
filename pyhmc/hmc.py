import os 
os.environ['OPENBLAS_NUM_THREADS'] = "1"
os.environ['MKL_NUM_THREADS'] = '1'
import numpy as np
import sys
import h5py


class HamitonianMC:
    def __init__(self,UserDefinedModel,boundaries,
            dt,Lrange,nbest_model,seed,
            nsamples,ndraws,
            myrank=0,name = "mychain",
            outdir = "./"):
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
        Lrange : List[int], shape(2)
            min/max steps forward
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
        
        self.myrank = myrank
        self.seed = seed + myrank
        nt = boundaries.shape[0]
        self.boundaries = boundaries
        self.Lrange = Lrange
        self.dt = dt
        self.invert_Mass = np.eye(nt)
        self.model = UserDefinedModel
        self.nbest_model = nbest_model
        self.nsamples = nsamples
        self.ndraws = ndraws

        self.cache = {}
        self.ii = 0

        # open hdf5 to save results
        self.fio = h5py.File(f"{outdir}/{name}.{self.myrank}.h5","w")

        # set random seed
        np.random.seed(self.seed)

    @classmethod
    def init(self,UserDefinedModel,boundaries,rank,**kargs):
        return HamitonianMC(UserDefinedModel,boundaries,
                            kargs['dt'],kargs['Lrange'],
                            kargs['nbest'],kargs['seed'],
                            kargs['nsamples'],
                            kargs['ndraws'],
                            rank,kargs['name'],
                            kargs['OUTPUT_DIR']
                )

    def set_initial_model(self):
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
    
    def misfit_and_grad(self,x):
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
        U,grad,dsyn,flag = self.misfit_and_grad(xnew)
        if flag == False or np.sum(np.isnan(dsyn)) > 0:
            return xcur,np.inf,self.model.dobs,False
        Hcur = K + U

        # save current potential and synthetics
        dsyn_new = dsyn.copy()
        Unew = U

        # update
        pnew -= dt * grad * 0.5 
        for i in range(L):
            xnew += dt * pnew # update xnew

            # check boundaries
            xnew,pnew = self._mirror(xnew,pnew)

            # update pnew
            if np.sum(np.isnan(xnew)) > 0:
                return xcur,np.inf,self.model.dobs,False
            Unew,grad,dsyn_new,flag = self.misfit_and_grad(xnew)
            if np.sum(np.isnan(grad)) > 0:
                # the grad might be NaN. check it.
                return xcur,np.inf,self.model.dobs,False            
            if flag == False or np.sum(np.isnan(dsyn_new)) > 0:
                return xcur,np.inf,self.model.dobs,False
            if i < L-1:
                pnew -= dt * grad
            else:
                pnew -= dt * grad * 0.5 

        # end for loop
         
        # update Hamiltonian
        pnew = -pnew
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
        self.fio.create_dataset("initmodel",data= x)

        # save obs
        self.fio.create_dataset("obs",data = self.model.dobs)

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
        ncount = 0
        i = 0
        while i < ndraws + nsamples:
            #np.random.seed(self.seed + self.myrank + i)
            L = np.random.randint(self.Lrange[0],self.Lrange[1]+1)
            x,U,dsyn,AcceptFlag = self._leapfrog(x,self.dt,L)
            if AcceptFlag: # accept this new sample
                if i>=ndraws:
                    #self.save_results(x,dsyn,i-ndraws)
                    misfit[i-ndraws] = U
                    x_cache[i-ndraws,:] = x.copy()
                    syndata[i-ndraws,:] = dsyn
                i += 1   
                self.ii += 1         
            ncount += 1
            if i > -1 and (i % 50 == 0 or i == nsamples - 1):
                msg = "chain {}: {:.2%}, misfit={:.3} -- accept ratio {:.2%}".  \
                    format(self.myrank,i/(ndraws + nsamples),U,i/ncount)   
                print(msg)
                sys.stdout.flush()
        # end i-loop

        # compute average of the best n models
        nbests = self.nbest_model
        idx = np.argsort(misfit)
        xmean = np.mean(x_cache[idx[:nbests],:],axis=0)
        _,_,dsyn,_ = self.misfit_and_grad(xmean)

        # save synthetics and each sample
        self.save_results(xmean,dsyn,"mean")
        for i in range(nsamples):
            self.save_results(x_cache[i,:],syndata[i,:],i)

        return misfit