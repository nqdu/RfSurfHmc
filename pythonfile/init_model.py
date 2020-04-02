import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d

def draw_model(thk,v,ax,color):
    n = len(thk)
    dep = np.zeros(n)
    dep[0] = 0.0
    for i in range(1,n):
        dep[i] = dep[i-1] + thk[i-1]
    
    func = interp1d(dep,v,'slinear')

    x0 = np.linspace(dep[0],dep[n-1],1000)
    y0 = func(x0)

    ax.plot(y0,x0,color=color)


def generate_synthetic_file():
    thk = [5.0,5.0,5.0,10.0,15.0,0.0]
    v = np.array([2.5,2.7,3.0,2.8,3.2,3.5])
    nz = len(thk)

    # interpolate function
    dep = np.zeros(nz)
    dep[0] = 0.0
    for i in range(1,nz):
        dep[i] = dep[i-1] + thk[i-1]
    
    func = interp1d(dep,v,'linear')

    #dep0 = dep.copy()
    nz0 = 30
    zmax = np.sum(thk)
    dep0 = np.linspace(0,zmax,nz0)
    thk0 = np.zeros(nz0)
    thk0[:nz0-1] = dep0[1:nz0] - dep0[0:nz0-1]
    thk0[nz0-1] = 0.0
    v1 = func(dep0) - 0.2
    v2 = func(dep0) + 0.2
    print(v1)

    v0 = np.linspace(2.5,3.55,nz0)

    # save true model
    f = open('input/true.dat','w')
    f.write('%d\n'%nz)
    for i in range(nz):
        f.write('%f %f\n'%(thk[i],v[i]))
    f.close()

    # 
    os.system('./bin/synthetic')
    
    # generate file
    np.random.seed(1000)
    f = open('input/hmc.in','w')
    f.write('%d %d\n'%(nz0,1))
    for i in range(nz0):
        f.write('%f %f %f\n'%(thk0[i],v1[i],v2[i]))
    for g in ['Rc','Rg','Lc','Lg']: 
        try:
            d = np.loadtxt(g + '.true')
            f.write('%d '% (d.shape[0]))
        except:
            f.write('%d '% (0))
    f.write('\n')
    for g in ['Rc','Rg','Lc','Lg']: 
        try:
            d = np.loadtxt(g + '.true')
            n = d.shape[0]
            d[:,-1] *= 1.0 +0.02 * np.random.randn(n)
            for j in range(n):
                f.write('%f %f\n'%(d[j,0],d[j,1]))
        except:
            pass
        
        #os.system('rm '+g +'.txt')

    f.write('0.01  0.02 50 100\n')
    f.write('100\n')
    f.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # show init model
    draw_model(thk,v,ax,'k')
    draw_model(thk0,v1,ax,'b')
    draw_model(thk0,v2,ax,'b')
    #plt.legend()
    plt.gca().invert_yaxis()
    plt.show()

generate_synthetic_file()
