import numpy as np 
from .lib import libsurf

class SurfWD:
    def __init__(self,mode=0,sphere=False,tRc=None,tRg=None,tLc = None,tLg=None):
        # init swd params
        self.mode = mode # mode, =0 for fundamental
        self.sphere = sphere  # False: flat model, True: spherical mod

        # init period vector
        self.tRc,self.tRg,self.tLc,self.tLg = None,None,None,None  
        self.nt = 0
        self.ntRc,self.ntRg,self.ntLc,self.ntLg = 0,0,0,0
        if tRc is not None and len(tRc) > 0:  
            self.tRc = np.asarray(tRc)
            self.ntRc = len(tRc)
            self.nt += self.ntRc 
        if tRg is not None  and len(tRg) > 0:  
            self.tRg = np.asarray(tRg)
            self.ntRg = len(tRg)
            self.nt += self.ntRg 
        if tLc is not None  and len(tLc) > 0:
            self.tLc = np.asarray(tLc)
            self.ntLc = len(tLc)
            self.nt += self.ntLc
        if tLg is not None  and len(tLg) > 0: 
            self.tLg = np.asarray(tLg)
            self.ntLg = len(tLg)
            self.nt += self.ntLg 

    @classmethod
    def init(self,**kargs):
        return SurfWD(
            tRc=kargs['tRc'],
            tRg=kargs['tRg'],
            tLc=kargs['tLc'],
            tLg=kargs['tLg']
        )
    
    def set_obsdata(self,dobs):
        # obs data 
        self.dobs = dobs
    
    def set_thk(self,thk):
        self.thk = thk* 1.0
    
    def empirical_relation(self,vs:np.ndarray,deriv=False):
        """
        compute vp/rho and derivative from empirical relations
        Parameters:
        ===================================================
        vs  : np.ndarray
            vs model
        deriv : bool
            if True,also return derivatives
        Returns:
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
        dadb = np.zeros(vp.shape)
        drda = np.zeros((vp.shape))

        if deriv:
            drda = 1.6612 - 0.4721*2*vp + 0.0671*3*vp**2 - 0.0043*4*vp**3 + 0.000106*5*vp**4
            dadb = 2.0947 - 0.8206*2*vs + 0.2683*3 * vs**2 - 0.0251*4*vs**3
        
        # return
        if deriv == False:
            return vp,rho 
        else:
            return vp,rho,dadb,drda

    def forward(self,x:np.ndarray):
        """
        compute surface wave dispersion for a given model x   
        Parameters:
        ==================================================
        x   : np.ndarray, shape(n*2,1) 
            n means layers
            input model, [vs thk]
        Returns:
        ==================================================
        d   : np.ndarray,shape(self.nt)
            output dispersion data
        """
        # allocate space
        d = np.zeros((self.nt)) 

        # prepare model
        layers = int(len(x) / 2)
        vs = x[:layers]
        thk = x[layers:]
        vp,rho = self.empirical_relation(vs)
        
        # compute dispersion
        if self.ntRc >0 :
            k1 = 0
            k2 = k1 + self.ntRc 
            d[k1:k2],flag = libsurf.forward(thk,vp,vs,rho,
                            self.tRc,"Rc",self.mode,self.sphere)
            if flag == False:
                return d,flag
        if self.ntRg > 0:
            k1 = self.ntRc 
            k2 = k1 + self.ntRg 
            d[k1:k2],flag = libsurf.forward(thk,vp,vs,rho,
                            self.tRc,"Rg",self.mode,self.sphere)
            if flag == False:
                return d,flag

        if self.ntLc > 0:
            k1 = self.ntRc + self.ntRg  
            k2 = k1 + self.ntLc 
            d[k1:k2],flag = libsurf.forward(thk,vp,vs,rho,
                            self.tRc,"Lc",self.mode,self.sphere)
            if flag == False:
                return d,flag
        if self.ntLg > 0:
            k1 = self.ntRc + self.ntRg + self.ntLc   
            k2 = k1 + self.ntLg 
            d[k1:k2],flag = libsurf.forward(thk,vp,vs,rho,
                            self.tRc,"Lg",self.mode,self.sphere)
            if flag == False:
                return d,flag
        return d,True    
    
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
        out : float
             misfit
        """
        d,flag = self.forward(x)
        if flag:
            return 0.5 * np.sum((d - self.dobs)**2),True
        else:
            return 0.0,flag
    

    def misfit_and_grad(self,x):
        """
        compute misfit and gradient for current model
        """
        # allocate space
        n = int(x.shape[0] / 2)
        nt = self.nt
        kernel = np.zeros((nt,n))
        kernel_thk = np.zeros((nt,n))
        d = np.zeros((nt))
        grad = np.zeros((n))

        # prepare model
        thk = x[int(len(x)/2):]
        vs = x[:int(len(x)/2)]
        vp,rho,dadb,drda = self.empirical_relation(vs,True)



        # compute sensitivity and kernel
        if self.ntRc > 0:
            k1 = 0
            k2 = self.ntRc
            cg,dcda,dcdb,dcdr,dcdh,flag = libsurf.adjoint_kernel(
                thk,vp,vs,rho,self.tRc,"Rc",self.mode,self.sphere
            )
            if flag == False:
                return 0.0,grad,d,flag
            d[k1:k2] = cg 
            kernel[k1:k2,:] = dcdb + dcda * dadb + dcdr * drda * dadb
            kernel_thk[k1:k2,:] = dcdh 
        if self.ntRg >0:
            k1 = self.ntRc
            k2 = k1 + self.ntRg
            cg,dcda,dcdb,dcdr,dcdh,flag = libsurf.adjoint_kernel(
                thk,vp,vs,rho,self.tRg,"Rg",self.mode,self.sphere
            )
            if flag == False:
                return 0.0,grad,d,flag
            d[k1:k2] = cg 
            kernel[k1:k2,:] = dcdb +dcda * dadb + dcdr * drda * dadb
            kernel_thk[k1:k2,:] = dcdh 
        if self.ntLc > 0:
            k1 = self.ntRc + self.ntRg
            k2 = k1 + self.ntLc
            cg,dcda,dcdb,dcdr,dcdh,flag = libsurf.adjoint_kernel(
                thk,vp,vs,rho,self.tRc,"Lc",self.mode,self.sphere
            )
            if flag == False:
                return 0.0,grad,d,flag
            d[k1:k2] = cg 
            kernel[k1:k2,:] = dcdb +dcda * dadb + dcdr * drda * dadb  
            kernel_thk[k1:k2,:] = dcdh 
        if self.ntLg > 0:
            k1 = self.ntRc + self.ntRg + self.ntLc
            k2 = k1 + self.ntLg
            cg,dcda,dcdb,dcdr,dcdh,flag= libsurf.adjoint_kernel(
                thk,vp,vs,rho,self.tRc,"Lg",self.mode,self.sphere
            )
            if flag == False:
                return 0.0,grad,d,flag
            d[k1:k2] = cg 
            kernel[k1:k2,:] = dcdb +dcda * dadb + dcdr * drda * dadb   
            kernel_thk[k1:k2,:] = dcdh 

        
        # compute gradient
        #grad = (d - self.dobs) @ kernel
        #grad = (d - self.dobs) @ kernel_thk
        grad = np.hstack(( (d - self.dobs) @ kernel , (d - self.dobs) @ kernel_thk ))
        # compute misfit
        out = 0.5 * np.sum((d - self.dobs)**2)

        return out,grad,d,True
"""
TEST CODE
===================================
thk = np.array([2,4,6,10,12,0])
tRc = np.linspace(4,40,19)
model = SurfWD(thk,tRc=tRc,mode=0)
x = np.array([2.0,2.3,2.8,3.4,3.5,4.0])
dobs = model.forward(x)
dobs = dobs * (1.0+ 0.01 * np.random.randn(len(dobs)))
model.set_obsdata(dobs)
plt.plot(tRc,dobs)
plt.show()
exit()
"""
