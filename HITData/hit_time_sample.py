#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  hit_sample.py
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
    
    You need certain input parameters to be able to reach a desired
    turbulence intensity inflow.
    
    First, you need to specify a turbulence intensity
    TI = this is the desired turbulence intensity as u'/U

    Then, you need to specify a time series. The code will sweep through the
    database at a sweeping velocity corresponding to the turbulence intensity
    You need to specity the initial and final times for sweeping
    t0 = initial time in the database
    tf = final time on the database
    
    The widht of the domain is specified as:
    Ly = spanwise direction (database width=2pi)
    Lz = spanwise direction (database width=2pi)
    Lx = Streamwise direction (computed according to the time sweeping)

    Grid points
    ny, nz = this is the number of grid points
    
    Filter width
    fw = the size of the filter width, by default it is the grid spacing
    '''

    # Desired turbulence intensity
    TI = 0.1

    # Make working directory and switch to it
    dir = 'DataTI10'
    if not os.path.exists(dir):
        os.makedirs(dir)
    os.chdir(dir)

    ## From the database ##
    pi = np.float32(np.pi)

    # Fluctuations in the database (u'_DB)
    up = 0.686

    # The sweeping velocity (u'_DB / [u'_LES/U_LES])
    Usweeping = up / TI

    # Total time sample
    #~ t0 = 0.
    #~ tf = 10.056
    t0 = 0.
    tf = 10.056

    # Domain Length from the HIT database
    Lx = np.float32((tf - t0) * Usweeping)
    Ly = Lz = 1.0

    # Size of the plain in grid points
    ny =  nz = 32
    nx = int(Lx * ny / Ly)
    #~ nx, ny, nz = 5, 5, 5

    # Time sampling
    time = np.linspace(t0, tf, nx)

    # LES Filter width (depends on ny spectral direction)
    fw = Ly / ny

    # Open file simulation details (fsd)
    fsd = open('simulationdetails.txt','w')
    fsd.write('Simulation Details' + '\n')
    fsd.write('up=' + str(up) + '\n')
    fsd.write('TI=' + str(TI) + '\n')
    fsd.write('Usweeping=' + str(Usweeping) + '\n')
    fsd.write('Lx=' + str(Lx) + '\n')
    fsd.write('Ly=' + str(Ly) + '\n')
    fsd.write('Lz=' + str(Lz) + '\n')
    fsd.write('Nx=' + str(nx) + '\n')
    fsd.write('Ny=' + str(ny) + '\n')
    fsd.write('Nz=' + str(nz) + '\n')
    fsd.write('Time start=' + str(t0) + '\n')
    fsd.write('Time final=' + str(tf) + '\n')
    fsd.close()

    # Define the size of the domain to be used
    x1 = Usweeping * time
    y1 = np.linspace(0, Ly, ny)
    z1 = np.linspace(0, Lz, nz)

    # Generate mesh
    x, y, z = np.meshgrid(x1, y1, z1, indexing='ij')

    # Generate the points
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

        print('Reading Data')

        # Load the numpy data
        u = np.load('./u' + fname + '.npy')
        v = np.load('./v' + fname + '.npy')
        w = np.load('./w' + fname + '.npy')

    # If not available extract the data from HIT database (JHU)
    else:

        print('Initializing Database')

        # load shared library
        lTDB = pyJHTDB.libJHTDB(auth_token='com.gmail.tonyinme-7a6d4581')
    
        #initialize webservices
        lTDB.initialize()

        for i in range(len(x1)):

            print('Calling Database for i=', i)

            # Get filtered Data
            uvw = lTDB.getBoxFilter(time[i], 
                points[i,:,:,:].reshape(-1, 3),
                field = 'velocity',
                filter_width = fw)

            #~ uvw = lTDB.getData(
                   #~ t,
                   #~ points[:,:,k,:].reshape(-1, 3),
                   #~ sinterp = 4,
                   #~ getFunction='getVelocity')

            # Asign the components of velocity
            u[i,:,:] = uvw[:,0].reshape(ny, nz)
            v[i,:,:] = uvw[:,1].reshape(ny, nz)
            w[i,:,:] = uvw[:,2].reshape(ny, nz)

            print('Filtered Velocity obtained for i=', i)

        # Save data into numpy format
        np.save('./u' + fname, u)
        np.save('./v' + fname, v)
        np.save('./w' + fname, w)
        
    # Write the data to a file fortran binary
    write_binary('./binary_u' + fname, u.swapaxes(0,2))
    write_binary('./binary_v' + fname, v.swapaxes(0,2))
    write_binary('./binary_w' + fname, w.swapaxes(0,2))

    write_lesgo_data('./lesgo_u' + fname, u.swapaxes(0,2))
    write_lesgo_data('./lesgo_v' + fname, v.swapaxes(0,2))
    write_lesgo_data('./lesgo_w' + fname, w.swapaxes(0,2))

    # Save figures
    plt.pcolormesh(y[0,:,:], z[0,:,:], u[0,:,:], shading='gouraud', cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    #~ plt.xlim([0., 2*np.pi])
    #~ plt.ylim([0., 2*np.pi])
    plt.savefig('u.jpg', dpi=500)
    plt.clf()

    plt.pcolormesh(x[:,3,:], z[:,3,:], u[:,3,:], shading='gouraud', cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    #~ plt.xlim([0., 2*np.pi])
    #~ plt.ylim([0., 2*np.pi])
    plt.savefig('v.jpg', dpi=500)
    plt.clf()

    plt.pcolormesh(x[:,:,1], y[:,:,1], u[:,:,1], shading='gouraud', cmap=cmaps.inferno)
    plt.colorbar()
    # Set the colorscale
    plt.gca().set_aspect('equal', adjustable='box')

    # Axis limits
    #~ plt.xlim([0., 2*np.pi])
    #~ plt.ylim([0., 2*np.pi])
    plt.savefig('w.jpg', dpi=500)

    # Write Vapor Data
    #~ write_vapor('vaporfile', Lx, Ly, Lz, nx, ny, nz, ['u', 'v', 'w'])

def write_binary(filename, u):
    '''
    Write the data in the proper binary format for vapor
    '''
    u.astype('float32').tofile(filename)

def write_lesgo_data(filename, u):
    '''
    Write the data in the proper binary format for lesgo
    '''
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
               filename+'.vdf ' + 'binary_'+var+fname]
        os.system(cmd2[0])
        print(cmd2[0])

if __name__ == '__main__':
	main()

