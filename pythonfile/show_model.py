import numpy as np 
import matplotlib.pyplot as plt 
from glob import glob
from scipy.interpolate import interp1d

plt.switch_backend('agg')
def draw_model(thk,v,ax,color):
    n = len(thk)
    dep = np.zeros(n)
    dep[0] = 0.0
    for i in range(1,n):
        dep[i] = dep[i-1] + thk[i-1]
    
    func = interp1d(dep,v,'linear')

    x0 = np.linspace(dep[0],dep[n-1],1000)
    y0 = func(x0)

    ax.plot(y0,x0,color=color)

# read model and dispersion
with open('input/hmc.in','r') as f:
    lines = f.readlines()
n = lines[0].split()[0]
n = int(n)
model = np.zeros((n,3))
for i,line in enumerate(lines[1:n+1]):
    thk,v0,v1 = list(map(lambda x: float(x),line.split()))
    model[i,0] = thk
    model[i,1] = v0
    model[i,2] = v1
    print(thk,v0,v1)
nc,ng = lines[n+1].split()[:2]
nc = int(nc)
ng = int(ng)
tc = np.zeros(nc)
tg = np.zeros(ng)
vc = tc * 1.0
vg = tg * 1.0
n += 1
for i,line in enumerate(lines[n+1:n+1+nc]):
    t,v = line.split()
    t = float(t)
    v = float(v)
    tc[i] = t
    vc[i] = v

for i,line in enumerate(lines[n+1+nc:n+1+nc+ng]):
    t,v = line.split()
    t = float(t)
    v = float(v)
    tg[i] = t
    vg[i] = v

# draw model
fig = plt.figure(1)
ax = fig.add_subplot(111)
draw_model(model[:,0],model[:,1],ax,'k')
draw_model(model[:,0],model[:,2],ax,'k')

d_mean = 0.0
n = 0
filenames = glob('Models/iter*')
for f in filenames:
    d = np.loadtxt(f)
    d_mean += d[:,1]
    n += 1
    #draw_model(d[:,0],d[:,1],ax,'k')
    #ax.step(d[:,1],d[:,0],color='k')
d_mean = d_mean / n
draw_model(d[:,0],d_mean,ax,'b')

ax.invert_yaxis()
ax.set_xlabel('velocity,km/s')
ax.set_ylabel('depth,km')
plt.savefig("model.pdf")

#draw Rayleigh phase
fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(tc,vc)
filenames = glob('Models/Rc*')
for f in filenames:
    d = np.loadtxt(f)
    ax.plot(d[:,0],d[:,1],color='k')

ax.set_ylabel('velocity,km/s')
ax.set_xlabel('period,s')
ax.set_title('Rayleigh phase')
plt.savefig("phase.pdf")
#ax.set_ylim([2.0,3.0]) 

#draw Rayleigh group
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(tg,vg)
filenames = glob('Models/Rg*')
for f in filenames:
    d = np.loadtxt(f)
    ax.plot(d[:,0],d[:,1],color='k')
#d = np.loadtxt('Rg0')
#ax.plot(d[:,0],d[:,1],color='k')
ax.set_ylabel('velocity,km/s')
ax.set_xlabel('period,s') 
ax.set_title('Rayleigh group')
#ax.set_ylim([2.0,3.0])   

'''
#draw Love phase
fig = plt.figure(4)
ax = fig.add_subplot(111)
d0 = np.loadtxt('Lc.true')
ax.plot(d0[:,0],d0[:,1])
filenames = glob('Models/Lc*')
#d = np.loadtxt('Lc.inv')
#ax.plot(d[:,0],d[:,1],color='k')
for f in filenames:
    d = np.loadtxt(f)
    ax.plot(d[:,0],d[:,1],color='k')
ax.set_ylabel('velocity,km/s')
ax.set_xlabel('period,s') 
ax.set_title('Love phase')
#ax.set_ylim([2.0,3.0]) 

#draw Love group
fig = plt.figure(5)
ax = fig.add_subplot(111)
d0 = np.loadtxt('Lg.true')
ax.plot(d0[:,0],d0[:,1])
filenames = glob('Models/Lg*')
#d = np.loadtxt('Lg.inv')
#ax.plot(d[:,0],d[:,1],color='k')
for f in filenames:
    d = np.loadtxt(f)
    ax.plot(d[:,0],d[:,1],color='k')
ax.set_ylabel('velocity,km/s')
ax.set_xlabel('period,s')
ax.set_title('Love group')
#ax.set_ylim([2.0,3.0]) 
'''
plt.savefig('group.pdf')
plt.show()
