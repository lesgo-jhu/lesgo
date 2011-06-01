#!/bin/bash

if [ $# -ne 3 ]
then
  echo "Usage: `basename $0` {nproc}"
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
