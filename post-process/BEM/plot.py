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

w=1.1                                       # Rotor speed (rad/s)
db=(61.6333-2.8667)/N                         # Length of each section
r=np.linspace(2.8667+db/2,61.6333-db/2.,N)       # Blade elements
V0=10.                                       # Reference velocity
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

# Type of airfoil identifier. Add one because arrays in matlab start at 1
airfoilType=np.rint(np.interp(r,rcta[:,0],rcta[:,3]));   

sigma=chord*B/(2*np.pi*r)                              # Solidity

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
Vrel=a.copy()
Ft=a.copy()
iter=0

#for airfoil in Airfoils:

while (tolerance > .00001): 
  if iter>100:
    break
  aold=a.copy()
  apold=ap.copy()
#  for airfoil in Airfoils:
  for i, radius in enumerate(r):
    all = np.loadtxt("./airfoilProperties/"+Airfoils[ int(airfoilType[i]) ] ) 
    alpha = all[:,0]   # Twist angle
    clList= all[:,1]   # Lift coefficient
    cdList= all[:,2]   # Drag coefficient
    # Step 2 Compute flow angle phi
    phi[i]=np.rad2deg(np.arctan( ((1.-a[i])*V0) / ((1.+ap[i])*w*r[i])))  # Incidence angle
    # Step 3 Compute local angle of attack alpha=phi-theta
    aoa[i]=(phi[i]-twist[i])
    # Step 4 Read Lift and drag coefficients
    cd[i]=np.interp(aoa[i],alpha, cdList)
    cl[i]=np.interp(aoa[i],alpha, clList)
    # Step 5 Compute Cn and Ct
    cn[i]=cl[i]*np.cos(np.deg2rad(phi[i]))+cd[i]*np.sin(np.deg2rad(phi[i]))
    ct[i]=cl[i]*np.sin(np.deg2rad(phi[i]))-cd[i]*np.cos(np.deg2rad(phi[i]))
    # Step 6 Calculate a and a'
    a[i]=1./(4.*np.sin(np.deg2rad(phi[i]))**2/(sigma[i]*cn[i])+1)
    ap[i]=1./(4.*np.sin(np.deg2rad(phi[i]))*np.cos(np.deg2rad(phi[i]))/(sigma[i]*cn[i])-1.);
    # Calculate tangential force
    Vrel[i]=w*r[i]*(1.+ap[i])/np.cos(np.deg2rad(phi[i]));
    Ft[i]=ct[i]*1./2.*rho*chord[i]*Vrel[i]**2;

  # Calculate the tolerance and store previous value of (a) for new iteration
  tolerance=np.sum(np.abs(aold-a))/N
  # Underrelax Solution
  a=aold+0.5*(a-aold)
  ap=apold+0.5*(ap-apold)
  iter+=1
print 'Number of iterations was ', iter

####################################################
# Save data
np.savetxt('Cl', np.array([r ,cl]).T)
np.savetxt('Cd', np.array([r ,cd]).T)
np.savetxt('alpha', np.array([r ,aoa]).T)
np.savetxt('Ft', np.array([r ,Ft]).T)


####################################################
# now, plot the data:

##### Cl
plt.plot(r, cl,'-',color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$C_L$')
plt.savefig('Cl.eps')

##### Ft
plt.clf()
plt.plot(r, Ft,'-',color='black')
plt.xlabel(r'$r$')
plt.ylabel(r'$F_T$')
plt.savefig('Ft.eps')

















