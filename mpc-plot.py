import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.style.use('std-srf')

dt = 6.7134375000000004

omega = np.loadtxt('output-mpc/omega.dat')
beta = np.loadtxt('output-mpc/beta.dat')
torque_gain = np.loadtxt('output-mpc/torque_gain.dat')
uhat = np.loadtxt('output-mpc/uhat.dat')
Ctp = np.loadtxt('output-mpc/Ctp.dat')
Cpp = np.loadtxt('output-mpc/Cpp.dat')
Pref = np.loadtxt('output-mpc/Pref.dat')
Pfarm = np.loadtxt('output-mpc/Pfarm.dat')

xl = [0,50]

print(np.size(Pref))

t = np.arange(0.0, (np.size(Pref,0)-0.0001)*dt, dt)

fig = plt.figure(figsize=(6.5,4))
ax = fig.add_subplot(3,3,1)
ax.plot(t/60,beta)
plt.ylabel(r'$\beta$ ($^\circ$)')
plt.xlim(xl)

ax = fig.add_subplot(3,3,2)
ax.plot(t/60,torque_gain*omega**2/1E6)
plt.ylabel(r'$T_g$ (MN-m)')
plt.xlim(xl)

ax = fig.add_subplot(3,3,3)
ax.plot(t/60,omega)
plt.ylabel(r'$\omega$ (rad/s)')
plt.xlim(xl)

ax = fig.add_subplot(3,3,4)
ax.plot(t/60,uhat)
plt.ylabel(r'$\hat{u}$ (m/s)')
plt.xlim(xl)

ax = fig.add_subplot(3,3,5)
ax.plot(t/60,Ctp)
plt.ylabel(r'$C_T^\prime$')
plt.xlim(xl)

ax = fig.add_subplot(3,3,6)
ax.plot(t/60,Cpp)
plt.ylabel(r'$C_P^\prime$')
plt.xlim(xl)

ax = fig.add_subplot(3,3,7)
ax.plot(t/60,omega/uhat*126/2)
plt.ylabel(r'$\lambda^\prime$')
plt.xlabel(r'time (min)')
plt.xlim(xl)

# ax = fig.add_subplot(3,3,8)
# plt.plot(t/60,alpha)
# ax.legend(['1','2','3','4'])
# plt.xlabel('time (min)')
# plt.ylabel(r'$\alpha$')
# plt.xlim(xl)

ax = fig.add_subplot(3,3,8)
lh1, = plt.plot(t/60,Pfarm/1E6,'k', label='Farm Power')
lh2, = plt.plot(t/60,Pref/1E6,'r--', label='Reference')
plt.xlabel('time (min)')
plt.ylabel('Power (MW)')
plt.xlim(xl)
ax = fig.add_subplot(3,3,9)
plt.legend(handles=[lh1, lh2])

plt.tight_layout()
# plt.show(
plt.savefig('mpc.pdf')
