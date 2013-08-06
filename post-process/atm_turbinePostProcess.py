#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#  
#  Copyright 2013 tony <tony@wind>
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
#  This code will plot turbine output as a function of time

import matplotlib.pyplot as plt
import numpy as np
import os

##########################################################
# The main function of this file
# This is the code that will be run
def main():

 # The files to be plotted (ylabel, divideby)
 files={'power':(r'P (W)',1000000.)}

 for file_in, parameters in files.iteritems():
   plotFile(file_in,parameters,file_in,'r',parameters[0])


##########################################################
# Function for plotting (single line)
def plotFile(file_in,parameters,run_label,labelX,labelY):
# Read file location
 read_location='../turbineOutput/'
 myFile=open(read_location+file_in)
 x,y=readFile(myFile)
 plt.clf()
 plt.plot(x,[float(k)/parameters[1] for k in y],'-o', color='black')
 plt.xlabel(labelX)
 plt.ylabel(labelY)
 plt.savefig(file_in+'.png')
 plt.clf()


##########################################################
# Read the file and average it
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
