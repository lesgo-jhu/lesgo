import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('out.dat')
fb = A[0,:]
fa = A[1,:]
fdb = A[2,:]
fda = A[3,:]

fig = plt.figure()
ax = fig.add_subplot(2,2,1)
ax.plot(fb,'k--')
ax.plot(fdb,'k')
plt.title('Beta')

ax = fig.add_subplot(2,2,2)
ax.plot(fa,'k--')
ax.plot(fda,'k')
plt.title('Alpha')

ax = fig.add_subplot(2,2,3)
ax.plot(abs(fb-fdb),'k')

ax = fig.add_subplot(2,2,4)
ax.plot(abs(fa-fda),'k')

plt.show()
