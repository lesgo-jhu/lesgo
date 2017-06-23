import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('std-srf')

dt = 0.33906250000000004

# x = np.loadtxt('output-startup/x.dat')
# y = np.loadtxt('output-startup/y.dat')
# u = np.loadtxt('output-startup/u.dat')
# u = np.reshape(u,[np.size(y),np.size(x)], order='C')
# X,Y = np.meshgrid(x,y)
# plt.pcolor(X,Y,u,cmap="seismic")
# plt.show()

omega = np.loadtxt('output-startup/omega.dat')
beta = np.loadtxt('output-startup/beta.dat')
gen_torque = np.loadtxt('output-startup/gen_torque.dat')
uhat = np.loadtxt('output-startup/uhat.dat')
Ctp = np.loadtxt('output-startup/Ctp.dat')
Cpp = np.loadtxt('output-startup/Cpp.dat')
P = 0.5*np.pi*63**2*Cpp*uhat**3
aero_torque = P / omega
thrust = 0.5*np.pi*63**2*Ctp*uhat**2

t = np.arange(0.0, np.size(omega,0)*dt, dt)

fig = plt.figure(figsize=(6.5,4))
ax = fig.add_subplot(3,3,1)
ax.plot(t/60,gen_torque/1E6)
plt.ylabel(r'$T$ (MN)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,2)
ax.plot(t/60,gen_torque/1E6)
plt.ylabel(r'$Q$ (MN-m)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,3)
ax.plot(t/60,omega)
plt.ylabel(r'$\omega$ (Hz)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,4)
ax.plot(t/60,uhat)
plt.ylabel(r'$\hat{u}$ (m/s)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,5)
ax.plot(t/60,Ctp)
plt.ylabel(r'$C_T^\prime$')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,6)
ax.plot(t/60,Cpp)
plt.ylabel(r'$C_P^\prime$')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,7)
ax.plot(t/60,omega/uhat*126/2)
plt.ylabel(r'$\lambda^\prime$')
plt.xlabel(r'time (min)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,8)
ax.plot(t/60,aero_torque/1E6)
plt.ylabel(r'$P/\omega$ (MN-m)')
plt.xlabel(r'time (min)')
plt.xlim([0,7])

ax = fig.add_subplot(3,3,9)
plt.plot(t/60,P/1E6)
plt.xlabel('time (min)')
plt.ylabel('Power (MW)')
plt.xlim([0,7])

plt.tight_layout()
plt.savefig('startup.pdf')