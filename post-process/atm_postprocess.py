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
##########################################################
# The main function of this file
##########################################################
def main():
 # Create directories to save plots
 if not os.path.exists('./blades'):
   os.makedirs('./blades')
	
 # The files to be plotted
 files={'alpha':(r"$\alpha$",0), 'Cl':(r"$C_{l}$",0), 
        'Cd' :(r"$C_{d}$",0)}

 for file_in, parameters in files.iteritems():
   plotFile(file_in,'r',parameters[0])


##########################################################
# Function for plotting (single line)
def plotFile(file_in,labelX,labelY):
# Read file location
 read_location='../turbineOutput/'
 myFile=open(read_location+file_in)
 x,y=readAndAverage(myFile)
 plt.clf()
 plt.plot(x,y,'-o',label='LES',color='black')
 x, y = np.loadtxt('./BEM/'+file_in, unpack=True)
 plt.plot(x,y,'-d',label='BEM',color='black')
 plt.xlabel(labelX)
 plt.ylabel(labelY)
 plt.legend()
 plt.savefig(file_in+'.png')
 plt.clf()

##########################################################
# Read the file and average it
def readAndAverage(file_in):
 for i, line in enumerate(file_in):
   if i==1:
     y=[float(s) for s in line.split()] # Initial value of y
     # Not sure what this does haha
   if i>1:
     y = [sum(pair) for pair in zip (y, [float(s) for s in line.split()])]
 y=[k/i for k in y] # Take the average 
 myfile=open('../turbineOutput/blade')
 line=myfile.readline()
 line=myfile.readline()
 x=[float(s) for s in line.split()] # Convert the elements in line to numbers
 del y[0],x[0]
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
# Run the main function
if __name__ == '__main__':
	main()


















