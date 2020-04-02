import numpy as np 
import matplotlib.pyplot as plt 

d = np.loadtxt('RC')
d1 = np.loadtxt('Rc.true')

plt.plot(d[:,0],d[:,1])
plt.plot(d1[:,0],d1[:,1])
print(d[:,1],d1[:,1])
plt.show()