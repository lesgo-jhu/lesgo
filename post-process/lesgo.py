#!/usr/bin/python
################################################################################
## Written by: 
##
##   Luis 'Tony' Martinez <tony.mtos@gmail.com> (Johns Hopkins University)
##
##   Copyright (C) 2012-2013, Johns Hopkins University
##
##   This file is part of Lesgo
##
##   Lesgo is free software: you can redistribute it 
##   and/or modify it under the terms of the GNU General Public License as 
##   published by the Free Software Foundation, either version 3 of the 
##   License, or (at your option) any later version.
##
##   Lesgo is distributed in the hope that it will be 
##   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
# Read lesgo.data and create the necessary data and plots from it
# The format read is CGNS which is read using the HDF5 library
# Many classes are contained which store and post-process the data
# main(): Run the program
# planeClass: a class for storing planes of data
# read_input_file(): a function for reading the input file lesgo.data
# fieldClass: a class which stores the fields to be plotted
# linesClass: a class for storing and operating on 1D data
# interp3(): function for interpolating linearly in 3D grid

import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    '''
    Main function 
    reads in the inpute file lesgo.data and
    creates the data files and plots specified in lesgo.data
    '''
    # Read the input file 'lego.dat'
    [extract_1D_data, extract_2D_data, plot_1D_data, plot_2D_data,
                      lines1D, planes, fieldList] =  read_input_file()

    # Extract all the data
    if extract_1D_data:
        print 'Extracting and writing data'
        # Loop through all the line plots
        for line in lines1D:
            # Extract the field data
            for field in fieldList:
                line.addData(File=field.File, field=field.name)
            line.writeData()
        print 'Done extracting and writing data'    

    if plot_1D_data:
        print 'Plotting Data'
        # Loop through all the line plots
        for line in lines1D:
            line.readData()  # Extract the field data
            for field in fieldList: # Plot each field
                line.plot(field)
        print 'Done plotting Data'
   
    if extract_2D_data:
        print 'Extracting 2D data'
        for plane in planes: # Extract each plane
            for field in fieldList: # Extract each field
                plane.extractData( fieldObj=field )
                plane.contour( fieldObj=field )
        print 'Done Extracting Data'

    if plot_2D_data:
        print 'Plotting 2D Data'
        for plane in planes: # Extract each plane
            for field in fieldList: # Extract each field
                plane.readData( fieldObj=field )
                plane.contour( fieldObj=field )
        print 'Done plotting 2D data'




class planeClass(object):
    '''
    Create a plane of data which is to be plotted
    
    Arguments:
    cutPlaneVector= points the normal vector to the plane
    pointPlane= A point in the plane
    maxVal= The maximum value for the legend
    minVal= The minimum value for the legend
    label= The text that goes in the legend
    upVector= The vector pointing up in the visualization

    functions:

    contour(self, file_name=):
    Creates the contour and stores it as an image

    '''

    version='0.5'

    def __init__(self, name=None, normal=None, dis=None):
        '''
        The initialization of the class object
        '''
        # Name of the plane
        self.name=name

        # Plane coordinate distance
        self.dis=dis

        # Normal plane
        self.normal=normal

        # Clipping of the domain
        self.clip=None
        # clip should be the limits in x, y, z
        # clip=[x0, x1, y0, y1, z0, z1]

        # Data coordinates
        self.x=None
        self.y=None
        self.z=None
        self.dx=None
        self.dy=None

        # The data from every field as a dictionary
        self.data={}

        # Create file to write Data if not existent
        if not os.path.exists('./Data'):
            os.makedirs('./Data')
        
    def extractData(self, fieldObj=None):
        '''
        This will provide the data field
        fieldObj is an object of the fieldClass
        '''
        field=fieldObj.name
        File=fieldObj.File
        # Point to the right location of the file
        File = '/home/tony/Dropbox/' + File
        f= h5py.File( File , "r")

        # Load coordinates and field
        path='/Base/Zone/GridCoordinates/'
        x_f = np.asarray(f[path+'CoordinateX/ data'])[0,0,:]
        y_f = np.asarray(f[path+'CoordinateY/ data'])[0,:,0]
        z_f = np.asarray(f[path+'CoordinateZ/ data'])[:,0,0]
        u = np.asarray(f['/Base/Zone/Solution/'+field+'/ data'])

        # For the plane
        self.x=np.unique(np.clip(x_f,self.clip[0],self.clip[1]))
        self.y=np.unique(np.clip(y_f,self.clip[2],self.clip[3]))
        self.z=np.unique(np.clip(z_f,self.clip[4],self.clip[5]))
        # Create the plane points
        if self.normal=='x': 
            self.x=[self.dis]; self.dx, self.dy=np.meshgrid(self.y,self.z)
            # Length of 2D array
            n1, n2=len(self.y), len(self.z)
        if self.normal=='y': 
            self.y=[self.dis]; self.dx, self.dy=np.meshgrid(self.x,self.z) 
            # Length of 2D array
            n1, n2=len(self.x), len(self.z)
        if self.normal=='z': 
            self.z=[self.dis]; self.dx, self.dy=np.meshgrid(self.x,self.y)
            # Length of 2D array
            n1, n2=len(self.x), len(self.y)

        self.data[field]=np.empty([n1,n2])

        if self.normal=='x': 
            x=self.x[0]
            for j, y in enumerate(self.y):
                for k, z in enumerate(self.z):
                    self.data[field][j,k]=interp3(x_f, y_f, z_f, u, x, y, z)
        if self.normal=='y': 
            y=self.y[0]
            for i, x in enumerate(self.x):
                for k, z in enumerate(self.z):
                    self.data[field][i,k]=interp3(x_f, y_f, z_f, u, x, y, z)
        if self.normal=='z': 
            z=self.z[0]
            for i, x in enumerate(self.x):
                for j, y in enumerate(self.y):
                    self.data[field][i,j]=interp3(x_f, y_f, z_f, u, x, y, z)

        np.savez('./Data/'+self.name+field, self.dx, self.dy,
                          self.data[field] )

    def readData(self, fieldObj=None):
        '''
        This will provide the data field
        fieldObj is an object of the fieldClass
        The data is read from numpy file (loadtxt)
        '''
        field=fieldObj.name
        npzfile= np.load('./Data/'+self.name+field+'.npz')
        self.dx = npzfile['arr_0']
        self.dy = npzfile['arr_1']
        self.data[field] = npzfile['arr_2']

    #################################
    def contour(self, fieldObj=None):
        '''
        Create contour
        '''
        # Obtain information form field object
        field=fieldObj.name
        label=fieldObj.label
        fMax=fieldObj.Max
        fMin=fieldObj.Min

        if not os.path.exists('./Contours'):
            os.makedirs('./Contours')
            
        for color in ['gray', 'hsv', 'spectral','gist_ncar','jet']:
            #~ matplotlib.rcdefaults()
            plt.clf()
            #~ print self.data[field].shape
            #~ print self.dx.shape, self.dy.shape
    
            # Adjust size of the figure
            #~ a,b=matplotlib.rcParams['figure.figsize']
            #~ mx=np.amax(self.dx)-np.amin(self.dx)
            #~ my=np.amax(self.dy)-np.amin(self.dy)
            #~ a*=max(mx,my)/mx
            #~ b*=max(mx,my)/my
            #~ matplotlib.rcParams['figure.figsize']=[a,b]
            #~ print a,b,matplotlib.rcParams['figure.figsize']
    
            plt.pcolormesh(self.dx, self.dy, self.data[field].T,
                           cmap = plt.get_cmap(color), vmin=fMin, vmax=fMax )
    
            # Set the colorscale
            plt.gca().set_aspect('equal', adjustable='box')
    
            # Axis limits
            plt.xlim([np.amin(self.dx), np.amax(self.dx)])
            plt.ylim([np.amin(self.dy), np.amax(self.dy)])
    
            # Tick [ara,eters
            plt.tick_params(\
                axis='both',       # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                left='off',        # ticks along the left edge are off
                right='off',       # ticks along the right edge are off
                labelbottom='off', # labels along the bottom edge are off
                labelleft='off') # labels along the bottom edge are off
    
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(plt.gca())
            cax = divider.append_axes("right", "2%", pad="3%")
    
            # Colorbar
            cb=plt.colorbar(cax=cax)
    
            # Set the numbe rof tick labels
            if fMax or fMin:
                cb.set_ticks(np.linspace( fMin, fMax, 3 ))
    
            #~ cb.set_label(label, rotation=0)
            
            # Figure attached to limits
            plt.tight_layout()
    
            plt.savefig('./Contours/'+self.name+field+'-'+color+'.jpg')
            plt.close()

######################        
def read_input_file():
    '''
    This function reads in the lesgo.data file
    The file contains all the information with respect to
    which data to extract and how to plot it
    '''
    # Variables to be read from file
    # point1 - coordinates (x1,y1,z1)
    # point2 - coordinates (x2,y2,z2)
    # N - number of points in line

    f=open('./lesgo.data','r')

    # List of line objects
    lines1D=[]

    # List of fields in domain
    fieldList=[]

    # List of planes
    planes=[]

    # Flags to specify which section of the input file we are on
    fields=None
    data_1D=None
    data_2D=None

    # Used to clip the domain
    clip=None
    
    # Loop for all the lines in the file
    for line in f:
        # Ignore commented lines
        if not line.startswith('#') and line.strip():

            # Store the quantities given in the file to the dictionary    
            line=line.strip()

            if line.split('=')[0]=='extract_1D_data':
                extract_1D_data=eval(line.split('=')[1])
            elif line.split('=')[0]=='extract_2D_data':
                extract_2D_data=eval(line.split('=')[1])
            elif line.split('=')[0]=='plot_1D_data':
                plot_1D_data=eval(line.split('=')[1])
            elif line.split('=')[0]=='plot_2D_data':
                plot_2D_data=eval(line.split('=')[1])

            elif line == 'fields:':
                fields='yes'
                data_1D=None
                data_2D=None

            elif line == '1D_data:':
                fields=None
                data_1D='yes'
                data_2D=None

            elif line == '2D_data:':
                fields=None
                data_1D=None
                data_2D='yes'
            
            # Read in the 1D data into lines class list
            elif fields and not data_1D and not data_2D:
                name= line.split(';')[0].strip()
                File= line.split(';')[1].strip()
                label= line.split(';')[2].strip()
                # Declare field class object in list               
                fieldList.append( fieldClass(name=name, label=label, File=File))

                # Maximum and minimum values, if given
                if len(line.split(';'))>3: 
                    fieldList[-1].Min= float(line.split(';')[3])
                if len(line.split(';'))>4: 
                    fieldList[-1].Max= float(line.split(';')[4])

            # Read in the 1D data into lines class list
            elif data_1D and not data_2D:                
                name=line.split(';')[0]
                p1=np.array(eval(line.split(';')[1]))
                p2=np.array(eval(line.split(';')[2]))
                N=eval(line.split(';')[3])
                # Declare the class object as part of list
                lines1D.append( linesClass( p1=p1, p2=p2, N=N, name=name ) )

            # Read in the 1D data into lines class list
            elif data_2D and not data_1D:      
                if line.split('=')[0].strip()=='clip':
                    clip=eval(line.split('=')[1])
                else:
                    # Store the quantities given in the file to the dictionary    
                    line=line.strip()
                    #~ print line
                    name_data=line.split(';')[0]
                    # The arguments passed as input to create class
                    args=line.split(';')[1]
                    # Create class object for plane
                    planes.append(eval('planeClass( ' + args + ' )') )
                    planes[-1].name=name_data
                    if clip: planes[-1].clip=clip

    result=[extract_1D_data, extract_2D_data, plot_1D_data, plot_2D_data,
            lines1D, planes, fieldList]
    return result

class fieldClass(object):
    '''
    Stores the field with its information
    '''

    version='0.1'

    def __init__(self, name=None, label=None, File=None ):
        # Name of the field
        self.name=name

        # String with label field
        self.label=label

        # File containing the field
        self.File=File

        # Maximum value of field
        self.Max=None

        # Minimum value of field
        self.Min=None
    
class linesClass(object):
    '''
    Holds the information for all the fields in the line plot
    '''

    version='0.1'

    def __init__(self, p1=None, p2=None, N=None, name=None):
        '''
        The initialization of the class object
        p1 point 1 of the line
        p2 point 2 of the line
        N total number of line points
        name name of the line
        '''
        self.p1=p1
        self.p2=p2
        self.N=N
        self.name=name
        self.data={} # The data list is a dictionary
        self.data['x']=  list( np.linspace( p1[0], p2[0],N ) )
        self.data['y']=  list( np.linspace( p1[1], p2[1],N ) )
        self.data['z']=  list( np.linspace( p1[2], p2[2],N ) )
        line_len=np.linalg.norm(p1-p2)  # line length
        self.data['arclength']= list( np.linspace(0,line_len ,N) )

    def addData(self, File=None, field=None):
        '''
        This will provide the data field
        '''
        # Point to the right location of the file
        File = '/home/tony/Dropbox/'+File
        f= h5py.File( File , "r")
        path='/Base/Zone/GridCoordinates/'
        x_f = np.asarray(f[path+'CoordinateX/ data'])[0,0,:]
        y_f = np.asarray(f[path+'CoordinateY/ data'])[0,:,0]
        z_f = np.asarray(f[path+'CoordinateZ/ data'])[:,0,0]

        u = np.asarray(f['/Base/Zone/Solution/'+field+'/ data'])

        # Intetpolate the  value of the point from field        
        dummy_array=[]
        for i in xrange(0,self.N):
            u_i= [interp3(x_f, y_f, z_f, u,
                self.data['x'][i], self.data['y'][i],self.data['z'][i])][0]
            dummy_array.append(u_i) 
        self.data[field]=dummy_array  # Save to object
        
    def writeData(self):
        '''
        Write the data to a file
        '''
        # Name of the data file
        name=self.name

        # Create directories to write mean data
        if not os.path.exists('./Data'):
            os.makedirs('./Data')
        f = open( './Data/'+name+'.dat' , "w")
    
        N=self.N

        # The first fields to be written
        field_list=['arclength','x', 'y', 'z']
    
        # The fields not wanted to be written
        unwanted=[ 'N', 'point1', 'point2']; unwanted.extend(field_list)

        # The list of sorted arrays
        field_list.extend([x for x in sorted(self.data) if x not in unwanted])

        # Write data to file
        for val in (field_list):
            f.write('{0: >12}'.format(str(val)))
            f.write(' ')
        f.write('\n')
        for i in range(0,N):
            for val in (field_list):
                num='{:3.8f}'.format(self.data[val][i])
                f.write( '{0: >12}'.format(str(num)) )
                f.write(' ')
            f.write('\n')
        f.close()
        
    def readData(self):
        '''
        Read the file stored and plot all of its quantities 
        The output is a dictionary plot_data
        which contains all the fields
        '''
        # Open the file
        f=open('./Data/'+self.name+'.dat','r')
        
        # Store the name of the fields for the dictionary
        fields=f.readline().split()
        
        # Create the dictionary
        plot_data={}

        # Create the fields data as empyt lists        
        for field in fields:
            self.data[field]=[]
        
        # Loop for all the lines in the file and store 
        # the arrays for each field
        for line in f:
            for i, field in enumerate(fields):        
                self.data[field].append(line.split()[i])
            

        
    def plot(self, field, **parameters):
        '''
        Function for plotting 
        field object from fieldClass   
        '''
        if not os.path.exists('./plots/'+self.name):
            os.makedirs('./plots/'+self.name)
            
        plt.clf()
        unwanted=['arclength', 'x', 'y', 'z']
        x=self.data['arclength']
        plt.clf()
        y=self.data[field.name]
        plt.plot( x,y,'-o',color='black')
        plt.xlabel(r'$r/D$')
        plt.ylabel(field.label)
        plt.savefig('./plots/'+self.name+'/'+field.name+'.eps')
        plt.clf()

def interp3(x, y, z, u, xi, yi, zi):
    ''' 
    Interpolate 3D function
    x(nx),y(ny),z(nz),u(nx,ny,nz) are arrays
    xi yi zi are single numbers onto which to interpolate 
    '''
    nz,ny,nx = np.shape(u)
    for i in range(nx-1):
        if x[i] <= xi and xi < x[i+1]:
            i0=i
            i1=i+1          
    for j in range(ny-1):
        if y[j] <= yi and yi < y[j+1]:
            j0=j
            j1=j+1      
    for k in range(nz-1):                
        if z[k] <= zi and zi < z[k+1]:
            k0=k
            k1=k+1
    xd=(xi-x[i0])/(x[i1]-x[i0])  
    yd=(yi-y[j0])/(y[j1]-y[j0])
    zd=(zi-z[k0])/(z[k1]-z[k0])
    c00=u[k0,j0,i0]*(1.-xd)+u[k0,j0,i1]*xd
    c10=u[k0,j1,i0]*(1.-xd)+u[k0,j1,i1]*xd
    c01=u[k1,j0,i0]*(1.-xd)+u[k1,j0,i1]*xd
    c11=u[k1,j1,i0]*(1.-xd)+u[k1,j1,i1]*xd
    c0=c00*(1.-yd)+c10*yd
    c1=c01*(1.-yd)+c11*yd
    c=c0*(1.-zd)+c1*zd
    return c


########################################
def combined_plot(compare_plot, params):
    '''
    Plot the list of combined files
    '''
    # Loop through all the files specified
    for field, items in labels.iteritems():
        plt.clf()
        for i, file_in in enumerate(params):    
            plot_data=read_plot_file(file_in)
            x=plot_data['arclength']
            y=plot_data[field]
            plt.plot( x,y,dashes[i],color='black',label=file_in)
        plt.xlabel('arclength')
        plt.ylabel(field)            
        #~ plt.legend(loc='best')
        # Place legend to the right
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Save figure with legend
        plt.savefig('./plots/'+compare_plot+'-'+field+'.eps',bbox_inches='tight')

# Run the main function
if __name__ == '__main__':
	main()










