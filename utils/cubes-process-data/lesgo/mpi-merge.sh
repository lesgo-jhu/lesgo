##
##  Copyright (C) 2011-2013  Johns Hopkins University
##
##  This file is part of lesgo.
##
##  lesgo is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  lesgo is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
##

#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` {nproc}"
  exit 1
fi

NPROC=$1

MCR=vel_avg.mcr

# Create string for data files
for (( i=0 ; i < $NPROC ; i++ ))
do 

 DATASET="$DATASET \"|MFBD|/vel_avg.dat.c$i\""

done
rpl -- '<DATASET>' "$DATASET" $MCR

# Create interpolation for each zone
for (( i=0 ; i < $NPROC ; i++ ))
do 
  echo "\$!LINEARINTERPOLATE" >> $MCR
  echo "  SOURCEZONES =  [$(($i+1))]" >> $MCR
  echo "  DESTINATIONZONE = $(($NPROC+1))" >> $MCR
  echo "  VARLIST =  [4-7]" >> $MCR
  echo "  LINEARINTERPCONST = 0" >> $MCR
  echo "  LINEARINTERPMODE = DONTCHANGE" >> $MCR
done

# Write file
echo '$!WRITEDATASET "|MFBD|/vel_avg.dat"' >> $MCR
echo "  INCLUDETEXT = NO" >> $MCR
echo "  INCLUDEGEOM = NO" >> $MCR
echo "  INCLUDECUSTOMLABELS = NO" >> $MCR
echo "  INCLUDEAUTOGENFACENEIGHBORS = YES" >> $MCR
echo "  ASSOCIATELAYOUTWITHDATAFILE = NO" >> $MCR
echo "  ZONELIST =  [$(($NPROC+1))]" >> $MCR
echo "  BINARY = NO" >> $MCR
echo "  USEPOINTFORMAT = NO" >> $MCR
echo "  PRECISION = 9" >> $MCR
echo "  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT" >> $MCR
echo "\$!RemoveVar |MFBD|" >> $MCR

tec360 -b -p $MCR

