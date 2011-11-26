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

