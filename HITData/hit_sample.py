#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hit.py
#  
#  Copyright 2015 tony <tony@mini-wind>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
import colormaps as cmaps

# Database Python import
import pyJHTDB
import pyJHTDB.dbinfo
import pyJHTDB.interpolator

def main():
    '''
    This code will obtain a time history over a plane of filtered 
    Homogeneous isotropic turbulence data from the JHU data base
    '''
    ## From the database ##
    pi = np.float32(np.pi)

    # total time
    tot_t = 2.048
    
    # The time at which to sample from the database
    t = tot_t / 2
    
    # Size of the plain in grid points
    nx, ny, nz = 92, 32, 32
    #~ nx, ny, nz = 64, 64, 64
    #~ nx, ny, nz = 15, 12, 16
    #~ nx, ny, nz = 1024, 192, 384
    #~ nx, ny, nz = 64, 64, 128

    # Domain Length from the HIT database
    Lx = 2.1
    Ly = 1.0
    Lz = 1.0

    # LES Filter width (depends on ny spectral direction)
    fw = Ly / ny

    # Define the size of the domain to be used
    x1 = np.linspace(0, Lx, nx)
    y1 = np.linspace(0, Ly, ny)
    z1 = np.linspace(0, Lz, nz)

    # Generate mesh
    x, y, z = np.meshgrid(x1, y1, z1, indexing='ij')

    # Generate the points
    #~ points = np.array((x, y, z), dtype = 'float32').reshape(-1, 3, order='F')
    points = np.zeros((nx, ny, nz, 3), dtype = 'float32', order='F')

    # Generate the velocity fields
    u = np.zeros((nx, ny, nz), dtype = 'float32', order='F')
    v = np.zeros((nx, ny, nz), dtype = 'float32', order='F')
    w = np.zeros((nx, ny, nz), dtype = 'float32', order='F')

    print('shape points', np.shape(points))
    print('shape x', np.shape(x))
    print('shape y', np.shape(y))
    print('shape z', np.shape(z))

    # Asign values to points
    for i in range(len(x1)):
        for j in range(len(y1)):
            for k in range(len(z1)):
                points[i,j,k,0]=x1[i]
                points[i,j,k,1]=y1[j]
                points[i,j,k,2]=z1[k]

    # Name of the filtered file
    fname = 'Filtered' + '_nx_' + str(nx) + '_ny_' + str(ny) + '_nz_' + str(nz)

    # Read the files from current directory if available
    if os.path.isfile('./u' + fname + '.npy'):
        #~ u = np.fromfile('./u'+fname, dtype=np.float32).swapaxes(0,2)
        #~ v = np.fromfile('./v'+fname, dtype=np.float32).swapaxes(0,2)
        #~ w = np.fromfile('./w'+fname, dtype=np.float32).swapaxes(0,2)
        # Save data into numpy format

        # Load the numpy data
        u = np.load('./u' + fname + '.npy')
        v = np.load('./v' + fname + '.npy')
        w = np.load('./w' + fname + '.npy')

    # If not available extract the data from HIT database (JHU)
    else:

        # load shared library
        lTDB = pyJHTDB.libJHTDB(auth_token='com.gmail.tonyinme-7a6d4581')
    
        #initialize webservices
        lTDB.initialize()

        print("Error NOT Here 1")
    
        for k in range(len(z1)):

            # Get filtered Data
            uvw = lTDB.getBoxFilter(t, points[:,:,k,:].reshape(-1, 3),
                field = 'velocity',
                filter_width = fw)

            #~ uvw = lTDB.getData(
                   #~ t,
                   #~ points[:,:,k,:].reshape(-1, 3),
                   #~ sinterp = 4,
                   #~ getFunction='getVelocity')

            # Asign the components of velocity
            u[:,:,k] = uvw[:,0].reshape(nx, ny)
            v[:,:,k] = uvw[:,1].reshape(nx, ny)
            w[:,:,k] = uvw[:,2].reshape(nx, ny)

            print('Filtered Velocity obtained for k=', k)

        # Save data into numpy format
        np.save('./u' + fname, u)
        np.save('./v' + fname, v)
        np.save('./w' + fname, w)
        
    # Write the data to a file fortran binary
    write_binary('./binary_u' + fname, u.swapaxes(0,2))
    write_binary('./binary_v' + fname, v.swapaxes(0,2))
    write_binary('./binary_w' + fname, w.swapaxes(0,2))

    # Save figures
    plt.pcolormesh(y[0,:,:], z[0,:,:], u[0,:,:], shading='gouraud', 
                    cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    plt.xlim([0., 2*np.pi])
    plt.ylim([0., 2*np.pi])
    plt.savefig('u.jpg', dpi=500)
    plt.clf()

    plt.pcolormesh(x[:,3,:], z[:,3,:], u[:,3,:], shading='gouraud', 
                    cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    plt.xlim([0., 2*np.pi])
    plt.ylim([0., 2*np.pi])
    plt.savefig('v.jpg', dpi=500)
    plt.clf()

    plt.pcolormesh(x[:,:,1], y[:,:,1], u[:,:,1], shading='gouraud', 
                    cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    plt.xlim([0., 2*np.pi])
    plt.ylim([0., 2*np.pi])
    plt.savefig('w.jpg', dpi=500)

    # Write Vapor Data
    #~ write_vapor('vaporfile', Lx, Ly, Lz, nx, ny, nz, ['u', 'v', 'w'])

def write_binary(filename, u):
    '''
    Write the data in the proper binary format for vapor
    '''
    #~ u.astype('float64').tofile(filename)
    np.savetxt(filename, u.reshape(-1, order='F'))
    
def write_vapor(filename, L_x, L_y, L_z, nx, ny, nz, variables):
    '''
    Write data in vapor format
    '''
    # minL used to scale dimension so the start at 1.0
    minL=min([L_x, L_y, L_z])
    cmd1=['vdfcreate -dimension '+str(nx)+'x'+str(ny) +'x'+str(nz)+
        ' -level 0 -extents 0:0:0:'+str(L_x/minL)+
        ':'+str(L_y/minL)+':'+str(L_z/minL)+
        ' -vars3d '+':'.join(variables)+' '+ filename+'.vdf']
    os.system(cmd1[0])
    print(cmd1[0])

    fname = 'Filtered' + '_nx_' + str(nx) + '_ny_' + str(ny) + '_nz_' + str(nz)

    for var in variables:
        cmd2=['raw2vdf -varname '+ var +' '+
               filename+'.vdf ' + var+fname]
        os.system(cmd2[0])
        print(cmd2[0])


if __name__ == '__main__':
	main()

