#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  atm_postprocess.py
#
#  Copyright 2013 tony <tony@latrobe-300>
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
#  This code will read blade output from the ATM Library
#  and make the required plots

import os
import matplotlib.pyplot as plt
import numpy as np



################################################################################
##                      The main function of this file                        ##
################################################################################
def main():

    ###############################################
    # Plot blade output
    ###############################################
    # Create directories to save plots
    if not os.path.exists('./plots'):
        os.makedirs('./plots')

    # Create directories to write mean data
    if not os.path.exists('./meanData'):
        os.makedirs('./meanData')

    # The files to be plotted (ylabel,dummy)
    files={'alpha':(r"$\alpha$",1.), 'Cl':(r"$C_{l}$",1.),
           'Cd' :(r"$C_{d}$",1.),'lift':(r"Lift",1.),'drag':(r"Drag",1.),
           'Vrel':(r"$V_{rel}$",1.), 'Vaxial':(r"$V_{axial}$",8.),
           'Vtangential':(r"$V_{tangential}$",1.) }

    for file_in, parameters in files.iteritems():
        plotFileBlade(file_in,'r',parameters)

    ###############################################
    # Plot turbine output
    ###############################################
    # The files to be plotted (ylabel, divideby)
    files={'power':(r'P (W)',1000000.)}

    for file_in, parameters in files.iteritems():
      plotFileTurbine(file_in,parameters,file_in,'t (s)',parameters[0])

##########################################################
# Function for writing mean data output to file
def writeMeanData(file_in,x,y):
    f = open( './meanData/'+file_in + '.dat', "w")
    for i,j in zip(x,y):
        f.write(str(i)+" "+str(j) + "\n")
    f.close()

##########################################################
# Function for plotting (blade properties)
def plotFileBlade(file_in,labelX,parameters):
# Read file location
 read_location='../turbineOutput/turbine1/'
 myFile=open(read_location+file_in)
 x,y=readAndAverage(myFile)

 # Write the plotted data to file
 writeMeanData(file_in,x,y)

 plt.clf()
 plt.plot([r/63.0 for r in x],  [r/parameters[1] for r in y],
          '-o',label='LES',color='black')
 del x,y
 x, y = np.loadtxt('./BEM/'+file_in, unpack=True)
 plt.plot([r/63.0 for r in x],  [r/parameters[1] for r in y] ,
          '--',label='BEM',color='black')
 plt.xlabel(labelX)
 plt.ylabel(parameters[0])
 plt.legend(loc='best')
 plt.savefig('./plots/'+file_in+'.eps')
 plt.clf()

##########################################################
# Function for plotting (power)
def plotFileTurbine(file_in,parameters,run_label,labelX,labelY):
# Read file location
 read_location='../turbineOutput/turbine1/'
 myFile=open(read_location+file_in)
 x,y=readFile(myFile)

 plt.clf()
 plt.plot(x,[float(k)/parameters[1] for k in y],'-o', color='black',label='LES')
 x, y = np.loadtxt('./BEM/'+file_in, unpack=True)
 plt.axhline(y/parameters[1],color='black',label='BEM')
 plt.xlabel(labelX)
 plt.ylabel(labelY)
 plt.legend(loc='best')
 plt.savefig('./plots/'+file_in+'.eps')
 plt.clf()

##########################################################
# Read the file and average it (for blade quantities)
def readAndAverage(file_in):
 j=500
 for i, line in enumerate(file_in):
   if i==j:
     y=[float(s) for s in line.split()] # Initial value of y
   if i>j:
     y = [sum(pair) for pair in zip (y, [float(s) for s in line.split()])]
 print i-j+1, 'Number of lines for', file_in
 y=[k/(i-j+1) for k in y] # Take the average
 myfile=open('../turbineOutput/turbine1/blade')
 line=myfile.readline()
 line=myfile.readline()
# Convert the elements in line to numbers
 x=[float(s) for s in line.split()] 
 del y[0:2],x[0:2]

 return x,y

##########################################################
# Load BEM data
def load(fname):
 ret = {}
 for line in open(fname, "rt"):
   k, v = line.strip().split("=", 1)
   ret[k] = eval(v)
 return ret

##########################################################
# Read the file and average it (For single number turbine output (power, etc))
def readFile(file_in):
 x,y=[],[]
 for i, line in enumerate(file_in):
   if i>0:
     x.append(i)
     y.append(line.split()[1]) # Initial value of y
 return x,y


##########################################################
# Run the main function
if __name__ == '__main__':
	main()


















