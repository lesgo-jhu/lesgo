import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('std-srf')

dt = 6.06375

omega = np.loadtxt('output-mpc/omega.dat')
beta = np.loadtxt('output-mpc/beta.dat')
gen_torque = np.loadtxt('output-mpc/gen_torque.dat')
uhat = np.loadtxt('output-mpc/uhat.dat')
Ctp = np.loadtxt('output-mpc/Ctp.dat')
Cpp = np.loadtxt('output-mpc/Cpp.dat')
Pref = np.loadtxt('output-mpc/Pref.dat')
Pfarm = np.loadtxt('output-mpc/Pfarm.dat')
alpha = np.loadtxt('output-mpc/alpha.dat')

t = np.arange(0.0, np.size(Pref,0)*dt, dt)

fig = plt.figure(figsize=(6.5,4))
ax = fig.add_subplot(3,3,1)
ax.plot(t/60,beta)
plt.ylabel(r'$\beta$ ($^\circ$)')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,2)
ax.plot(t/60,gen_torque/1E6)
plt.ylabel(r'$T_g$ (MN-m)')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,3)
ax.plot(t/60,omega)
plt.ylabel(r'$\omega$ (Hz)')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,4)
ax.plot(t/60,uhat)
plt.ylabel(r'$\hat{u}$ (m/s)')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,5)
ax.plot(t/60,Ctp)
plt.ylabel(r'$C_T^\prime$')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,6)
ax.plot(t/60,Cpp)
plt.ylabel(r'$C_P^\prime$')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,7)
ax.plot(t/60,omega/uhat*126/2)
plt.ylabel(r'$\lambda^\prime$')
plt.xlabel(r'time (min)')
plt.xlim([0,40])

ax = fig.add_subplot(3,3,8)
plt.plot(t/60,alpha)
plt.xlabel('time (min)')
plt.ylabel(r'$\alpha$')
plt.xlim([0,40])


ax = fig.add_subplot(3,3,9)
plt.plot(t/60,Pfarm/1E6,'k', label='Farm Power')
plt.plot(t/60,Pref/1E6,'r--', label='Reference')
plt.xlabel('time (min)')
plt.ylabel('Power (MW)')
plt.xlim([0,40])
plt.legend()

plt.tight_layout()
# plt.show(
plt.savefig('mpc.pdf')