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

if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` {nproc} {Lx} {Nx}"
  exit 1
fi

NPROC=$1
Lx=$2
Nx=$3
Dx=`echo "$Lx/$Nx" | bc`
Dx3=`echo "$Dx*$Dx*$Dx" | bc`

MCR=force_tot.mcr

# Create string for data files
for (( i=0 ; i < $NPROC ; i++ ))
do 

  if [ "$NPROC" -eq "1" ]; then
    DAT=force_avg.dat
  else
    DAT=force_avg.dat.c$i
  fi

   DATASET="$DATASET \"|MFBD|/$DAT\""

done
rpl -- '<DATASET>' "$DATASET" $MCR
rpl -- '<NPROC>' "$NPROC" $MCR

tec360 -b -p $MCR

FTOT_RAW=`tail -n 1 Out.txt | awk '{ print $2 }'`
#FTOT=`echo "$FTOT_RAW*1" | bc`
#echo $FTOT
rm -f force_tot.dat
echo "$FTOT_RAW" >> force_tot.dat
