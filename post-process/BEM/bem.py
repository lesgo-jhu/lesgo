#This script will do BEM calculations
import matplotlib.pyplot as plt
import numpy as np 

# Function to save data file
def save(fname, **kwargs):
  f = open(fname, "wt")
  for k, v in kwargs.items():
    print >>f, "%s=%s" % (k, repr(v))
  f.close()

N= 100                                        # Number of blade segments
tolerance=1.                                 # Initialize tolerance

#w=1.1                                       # Rotor speed (rad/s)
w=9.155                                      # Rotor speed (rev/min)
w*=np.pi/30.
db=(63.-1.5)/N                         # Length of each section
r=np.linspace(1.5+db/2,63.-db/2.,N)       # Blade elements
#r=np.linspace(1.5,63.,N)       # Blade elements
V0=8.                                       # Reference velocity
B=3.                                        # Number of blades
rho= 1.                                     # Air density

# NREL 5MW Reference Turbine Properties
# radius(m) chord(m) twist(deg) airfoil type
rcta = np.array(  [  [2.8667,     3.542 ,   13.308 ,   0],
                     [5.6     ,   3.854 ,  13.308  ,   0],
                     [8.3333 ,    4.167  ,  13.308   ,  1],
                     [11.75  ,    4.557  ,  13.308  ,   2],
                     [15.85  ,    4.652  ,  11.48   ,   3],
                     [19.95  ,    4.458  ,  10.162  ,   3],
                     [24.05  ,    4.249  ,  9.011   ,   4],
                     [28.15  ,    4.007  ,  7.795  ,    5],
                     [32.25  ,    3.748  ,  6.544  ,    5],
                     [36.35  ,    3.502  ,  5.361  ,    6],
                     [40.45  ,    3.256  ,  4.188  ,    6],
                     [44.55  ,    3.01   ,  3.125   ,   7],
                     [48.65  ,    2.764  ,  2.319  ,    7],
                     [52.75  ,    2.518  ,  1.526   ,   7],
                     [56.1667 ,   2.313  ,  0.863   ,   7],
                     [58.9   ,    2.086  ,  0.37    ,   7],
                     [61.6333 ,   1.419  ,  0.106   ,   7] ] )

# Chord for each element
chord=np.interp(r,rcta[:,0],rcta[:,1])

# Twist angle at each element
twist=np.interp(r,rcta[:,0],rcta[:,2])

# Type of airfoil identifier. 
#airfoilType=np.rint(np.interp(r,rcta[:,0],rcta[:,3]));   
airfoilType = np.rint(np.interp(r,rcta[:,0],rcta[:,3])) 

sigma=chord*B/(2.*np.pi*r)                              # Solidity

# Airfoil types in NREL 5MW Reference
Airfoils= [ "Cylinder1", 
    "Cylinder2" ,
    "DU40_A17" ,
    "DU35_A17" , 
    "DU30_A17" , 
    "DU25_A17" , 
    "DU21_A17" , 
    "NACA64_A17" ]

# Step 1
# Induction factor
a=np.zeros( (N,1) )
ap=a.copy()

phi=a.copy()
aoa=a.copy()
cl=a.copy()
cd=a.copy()
cn=a.copy()
ct=a.copy()
lift=a.copy()
drag=a.copy()
Vrel=a.copy()
Ft=a.copy()
iter=0

# Variables for calculating power
Ai = np.zeros( (N-1,1) )  
Bi = Ai.copy()  


#for airfoil in Airfoils:

while (tolerance > .0000001): 
  if iter>100:
    break
  aold=a.copy()
  apold=ap.copy()
#  for airfoil in Airfoils:
  for i, radius in enumerate(r):

    # Step 2 Compute flow angle phi
    phi[i]=np.rad2deg(np.arctan( ((1.-a[i])*V0) / ((1.+ap[i])*w*r[i])))
    # Step 3 Compute local angle of attack alpha=phi-theta
    aoa[i]=(phi[i]-twist[i])
    # Step 4 Read Lift and drag coefficients
    all = np.loadtxt("./airfoilProperties/"+Airfoils[ int(airfoilType[i]) ] ) 
    alpha = all[:,0]   # Twist angle
    clList= all[:,1]   # Lift coefficient
    cdList= all[:,2]   # Drag coefficient
    cd[i]= np.interp(aoa[i], alpha, cdList)
    cl[i]= np.interp(aoa[i], alpha, clList)

#    if int(airfoilType[i])>=max(airfoilType):
#        alpha_2 = all[:,0]   # Twist angle
#        clList_2= all[:,1]   # Lift coefficient
#        cdList_2= all[:,2]   # Drag coefficient
#        cd[i]= np.interp(aoa[i], alpha, cdList)
#        cl[i]= np.interp(aoa[i], alpha, clList)
#    else:
#        all = np.loadtxt("./airfoilProperties/"+Airfoils[ int(airfoilType[i+1]) ] ) 
#        alpha_2 = all[:,0]   # Twist angle
#        clList_2= all[:,1]   # Lift coefficient
#        cdList_2= all[:,2]   # Drag coefficient
#        cd[i]=np.interp( aoa[i], [r[i], r[i+1]],
#                     [ np.interp(aoa[i], alpha, cdList)[0],
#                     np.interp(aoa[i], alpha_2, cdList_2)[0]] )
#        cl[i]=np.interp( aoa[i], [r[i], r[i+1]],
#                     [ np.interp(aoa[i], alpha, clList)[0],
#                     np.interp(aoa[i], alpha_2, clList_2)[0]] )
    # Step 5 Compute Cn and Ct
    cn[i]=cl[i]*np.cos(np.deg2rad(phi[i]))+cd[i]*np.sin(np.deg2rad(phi[i]))
    ct[i]=cl[i]*np.sin(np.deg2rad(phi[i]))-cd[i]*np.cos(np.deg2rad(phi[i]))

    # Step 6 with tip-loss correction
    ftip = B/2 * ( 63.-r[i] ) / ( r[i] * np.sin( np.deg2rad ( phi[i] ) )) 
    Ftip = 2./np.pi * np.arccos( np.exp( -ftip ) )
    fhub = B/2 *  (63.-r[i] ) / ( 63. * np.sin( np.deg2rad ( phi[i] ) )) 
    Fhub = 2./np.pi * np.arccos( np.exp( -fhub ) )

    F = Ftip * Fhub
#    F=1.
    a[i]=1./(4.*F*np.sin(np.deg2rad(phi[i]))**2./(sigma[i]*cn[i])+1)
    ap[i]=1./(4.*F*np.sin(np.deg2rad(phi[i]))*np.cos(np.deg2rad(phi[i]))/
          (sigma[i]*ct[i])-1.)

    # Calculate tangential force
    Vrel[i]=w*r[i]*(1.+ap[i])/np.cos(np.deg2rad(phi[i]));
#    Ft[i]=ct[i]*1./2.*rho*chord[i]*Vrel[i]**2
    lift[i]=0.5*Vrel[i]**2 * chord[i] *cl[i]
    drag[i]=0.5*Vrel[i]**2 * chord[i] *cd[i]
    Ft[i] = ( lift[i] * np.sin( np.deg2rad ( phi[i] ) ) -
              drag[i] * np.cos( np.deg2rad ( phi[i] ) ) )
  # Calculate the tolerance and store previous value of (a) for new iteration
  tolerance=np.sum(np.abs(aold-a))/N
  # Underrelax Solution
  a=aold+0.5*(a-aold)
  ap=apold+0.5*(ap-apold)
  iter+=1
Vaxial=V0*(1.-a)
Vtangential=((w*r).T*(1.+ap).T).T

#M=0  # Moment
#for i, ai in enumerate(Ai):
#    Ai[i]= ( Ft[i+1] - Ft[i] ) / ( r[i+1] - r[i] ) 
#    Bi[i]= ( Ft[i] * r[i+1] - Ft[i+1] * r[i] ) / ( r[i+1] - r[i] ) 
#    M += ( 1./3 * Ai[i] * (r[i+1]**3 - r[i]**3) + 
#           1./2 * Bi[i] * (r[i+1]**2 - r[i]**2) )
# Calculate power output

power=np.sum(Ft.T*r.T)*w*B*db
#power=M*w*B

print 'Number of iterations was ', iter
print 'Power is ', power

####################################################
# Save data
np.savetxt('power', np.array([power, power]).T)
np.savetxt('RotSpeed', np.array([w, w]).T)

np.savetxt('Cl', np.array([r ,cl]).T)
np.savetxt('Cd', np.array([r ,cd]).T)
np.savetxt('alpha', np.array([r ,aoa]).T)
np.savetxt('Ft', np.array([r ,Ft]).T)
np.savetxt('lift', np.array([r ,lift]).T)
np.savetxt('drag', np.array([r ,drag]).T)
np.savetxt('Vrel', np.array([r ,Vrel]).T)
np.savetxt('Vaxial', np.array([r ,Vaxial]).T)
np.savetxt('Vtangential', np.array([r ,Vtangential]).T)


####################################################
# now, plot the data:

##### Cl
#~ plt.plot(r, cl,'-',color='black')
#~ plt.xlabel(r'$r$')
#~ plt.ylabel(r'$C_L$')
#~ plt.savefig('Cl.eps')

##### Ft
#~ plt.clf()
#~ plt.plot(r, Ft,'-',color='black')
#~ plt.xlabel(r'$r$')
#~ plt.ylabel(r'$F_T$')
#~ plt.savefig('Ft.eps')

















