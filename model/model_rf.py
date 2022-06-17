import numpy as np 
import librf 
class ReceiverFunc:
    def __init__(self, ray_p, nt, dt, gauss,
                    time_shift, water_level, rf_type = "p"):
        """
        Initialize Rceiver Function
        """
        self.ray_p = ray_p
        self.nt = nt
        self.dt = dt
        self.gauss = gauss
        self.time_shift = time_shift
        self.water_level = water_level
        self.rf_type = rf_type
        self.t = np.arange(nt) * dt - time_shift

    def set_obsdata(self, dobs):
        self.dobs = dobs

        """
        if np.sum(self.t - t) > 0.1 * self.dt:
            print("warning! the time sequence is not match!")
        """

    def set_thk(self,thk):
        self.thk = thk * 1.0
    

    def empirical_relation(self,vs,deriv=False):
        """
        compute vp/rho and derivative from empirical relations
        Parameters:
        ===================================================
        vs  : np.ndarray
            vs model
        ==================================================
        vp  : np.ndarray
            vp
        rho : np.ndarray
            rho      
        """  
        vp = 0.9409 + 2.0947*vs - 0.8206*vs**2 + \
            0.2683*vs**3 - 0.0251*vs**4
        rho = 1.6612*vp- 0.4721*vp**2 + \
             0.0671*vp**3 - 0.0043*vp**4 + 0.000106*vp**5
        if deriv:
            drda = 1.6612 - 0.4721*2*vp + 0.0671*3*vp**2 - 0.0043*4*vp**3 + 0.000106*5*vp**4
            dadb = 2.0947 - 0.8206*2*vs + 0.2683*3 * vs**2 - 0.0251*4*vs**3
            return vp,rho,drda,dadb
        else:
            return vp, rho

    def forward(self,x:np.ndarray):
        """
        compute rf for a given model x   
        Parameters:
        ==================================================
        x   : np.ndarray, shape(n*2,1) 
            n means layers
            input model, [vs thk]
        Returns:
        ==================================================
        d   ： np.ndarray,shape(self.nt)
            output dispersion data
        """
        # allocate spac

        # prepare model
        layers = int(len(x) / 2)
        vs = x[:layers]
        thk = x[layers:]
        

        
        vp,rho = self.empirical_relation(vs,deriv=False)
    
  

        # compute rf
        d = librf.forward(thk, rho, vp, vs, self.ray_p, self.nt, self.dt,
                        self.gauss, self.time_shift, self.water_level, self.rf_type)

        # bug here
        # sometimes the d will be 0. 
        # It might cause by the initialization of fortran code
        if np.sum(d) == 0:
           d = librf.forward(thk, rho, vp, vs, self.ray_p, self.nt, self.dt,
                        self.gauss, self.time_shift, self.water_level, self.rf_type) 

        return d

    def misfit(self,x):
        """
        compute misfit function 
        Parameters:
        ==================================================
        x   : np.ndarray, shape(n*2,1) 
            n means layers
            input model, [vs thk]
        Returns:
        ==================================================
        out ： float
             misfit
        """
        d = self.forward(x)
        return 0.5 * np.sum((d - self.dobs)**2)

    def misfit_and_grad(self,x):
        """
        compute gradient for current model
        """
        n = int(x.shape[0] / 2)
        nt = self.nt
        kernel = np.zeros((nt,n))
        kernel_thk = np.zeros((nt,n))
        d = np.zeros((nt))
        grad = np.zeros((n))
        


        # prepare model
        thk = x[int(len(x)/2):]
        vs = x[:int(len(x)/2)]
        vp,rho,drda,dadb = self.empirical_relation(vs,True)
        drda = drda.reshape(len(vs),1)
        dadb = dadb.reshape(len(vs),1)

        # compute gradient
        d,kvs = librf.adjoint_kernel(thk,rho,vp,vs,
                    self.ray_p, self.nt, self.dt, self.gauss, self.time_shift,
                    self.water_level, self.rf_type, 'vs')
        _,krho = librf.adjoint_kernel(thk,rho,vp,vs,
                    self.ray_p, self.nt, self.dt, self.gauss, self.time_shift,
                    self.water_level, self.rf_type, 'rho')
        _,kvp = librf.adjoint_kernel(thk,rho,vp,vs,
                    self.ray_p, self.nt, self.dt, self.gauss, self.time_shift,
                    self.water_level, self.rf_type, 'vp')

        _,kthk = librf.adjoint_kernel(thk,rho,vp,vs,
                    self.ray_p, self.nt, self.dt, self.gauss, self.time_shift,
                    self.water_level, self.rf_type, 'thick')        

        kernel = kvs + dadb * kvp + drda * dadb * krho  
        kernel_thk = kthk

        #grad = kernel @ (d - self.dobs) 
        grad = np.hstack((  kernel @ (d - self.dobs) ,  kernel_thk @ (d - self.dobs) ))


        misfit = 0.5 * np.sum((d - self.dobs)**2)

        return misfit,grad,d