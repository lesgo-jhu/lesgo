#!/usr/bin/python
# Read lesgo.data and create the necessary data and plots from it
# The format read is CGNS which is read using the HDF5 library
# MATLAB >> h5read('output_10.cgns','/Base/Zone/Solution/VelocityX/ data')
# data - this stores the data for all plots (its index is data_file
# data_file - name of data file (read from lesgo.data)
# point1 point2 - points that start the line
# N - number of points in line, arclength- the array of sampling
# xyz - all the points in x,y, and z access as: [0] for x, [1] for y, etc

import os
from array import *
import h5py
import numpy as np

def main():

    # Read the input file 'lego.dat'
    read_input_file()

    # Post process the data and write to file
    post_process()






def read_input_file():
    f=open('../lesgo.data','r')
    global data
    global fields
    fields=['VelocityX','VelocityY','VelocityZ']
    data={}
    # Variables to be read from file
    # point1 - coordinates (x1,y1,z1)
    # point2 - coordinates (x2,y2,z2)
    # N - number of points in line
    for line in f:
        if not line.startswith('#') and line.strip():
            line=line.strip()
            #~ print line
            name_data=line.split(';')[0]
            data[ name_data ]= {}
            data[ name_data ] ['point1']  = np.array(eval(line.split(';')[1]))
            data[ name_data ] ['point2']  = np.array(eval(line.split(';')[2]))
            data[ name_data ] ['N']  = eval(line.split(';')[3])
            data[ name_data ]['xyz']= {}

def post_process():

    file_name = '../output/output_10.cgns'
    
    f= h5py.File( file_name , "r")
    
    x_f = np.asarray(f['/Base/Zone/GridCoordinates/CoordinateX/ data'])[0,0,:]
    y_f = np.asarray(f['/Base/Zone/GridCoordinates/CoordinateY/ data'])[0,:,0]
    z_f = np.asarray(f['/Base/Zone/GridCoordinates/CoordinateZ/ data'])[:,0,0]
        
    u = np.asarray(f['/Base/Zone/Solution/VelocityX/ data'])
    v = np.asarray(f['/Base/Zone/Solution/VelocityY/ data'])
    w = np.asarray(f['/Base/Zone/Solution/VelocityZ/ data'])

    for name_data in data:
        # Create variables to be used
        p1=data[name_data]['point1']
        p2=data[name_data]['point2']
        N=data[name_data]['N']
        data[name_data]['xyz'][0]=  np.linspace( p1[0], p2[0],N )      
        data[name_data]['xyz'][1]=  np.linspace( p1[1], p2[1],N ) 
        data[name_data]['xyz'][2]=  np.linspace( p1[2], p2[2],N )
        line_len=np.linalg.norm(p1-p2)  # line length
        x=data[name_data]['arclength']=np.linspace(0,line_len ,N)

        for field in fields:
            data[name_data][field]=[]

            for i in range(N):
                data[name_data][field].append( interp3(x_f, y_f, z_f,  u,
                data[name_data]['xyz'][0][i], data[name_data]['xyz'][1][i],
                data[name_data]['xyz'][2][i]) )
        
        writeData(name_data,'y U, V, W',x,data[name_data]['VelocityX'])

def interp3(x, y, z, u, xi, yi, zi):
    " Interpolate 3D function"
    "x,y,z,u are (nx,ny,nz) arrays"
    " xi yi zi are single numbers onto which to interpolate "

    nz,ny,nx = np.shape(u)
    for i in range(nx-1):
        if x[i] <= xi and xi <= x[i+1] or x[i] >= xi and xi >= x[i+1]:
            xd=(xi-x[i])/(x[i+1]-x[i])  
            i1=i
            i2=i+1          
    for j in range(ny-1):
        if y[j] <= yi and yi <= y[j+1] or y[j] >= yi and yi >= y[j+1] :
            yd=(yi-y[j])/(y[j+1]-y[j])
            j1=j
            j2=j+1      
    for k in range(nz-1):                
        if z[k] <= zi and zi <= z[k+1] or z[k] >= zi and zi >= z[k+1] :
            zd=(zi-z[k])/(z[k+1]-z[k])
            k1=k
            k2=k+1

    c00=u[k1,j1,i1]*(1-xd)+u[k2,j1,i1]*xd
    c10=u[k1,j2,i1]*(1-xd)+u[k2,j2,i1]*xd
    c01=u[k1,j1,i2]*(1-xd)+u[k2,j1,i2]*xd
    c11=u[k1,j2,i2]*(1-xd)+u[k2,j2,i2]*xd
    c0=c00*(1-yd)+c10*yd
    c1=c01*(1-yd)+c11*yd
    c=c0*(1-zd)+c1*zd
    return c


def writeData(file_in,header,x,y):
    # Create directories to write mean data
    if not os.path.exists('./Data'):
        os.makedirs('./Data')
    f = open( './Data/'+file_in+'.dat' , "w")
    f.write(header+"\n")
    for i,j in zip(x,y):
        f.write(str(i)+" "+str(j) + "\n")
    f.close()


def plotFile(file_in,x,y, parameters):
##########################################################
# Function for plotting    
    plt.clf()
    x, y = np.loadtxt(file_in, unpack=True)
    plt.plot( x[1:], y[1:],'-o',color='black')
    plt.xlabel(x[0])
    plt.ylabel(y[0])
    plt.legend(loc='best')
    plt.savefig('./plots/'+file_in+'.eps')
    plt.clf()

# Run the main function
if __name__ == '__main__':
	main()










