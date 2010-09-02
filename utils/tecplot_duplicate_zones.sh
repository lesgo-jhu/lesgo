#!/bin/bash
echo "Enter output macro string:"
read OUTMCR
#echo "Enter TECPATH string:"
#read TECPATH
TECPATH=$(pwd)
echo "Enter number of zones to duplicate:"
read NZONE
echo "Enter Tecplot equation to move duplicated data:"
read EQUATION

echo '#!MC 1200' > $OUTMCR
echo '# Created by Tecplot 360 build 12.2.0.9077' >> $OUTMCR
echo "\$!VarSet |MFBD| = '$TECPATH'" >> $OUTMCR

for (( ZONEID=1; ZONEID<=$NZONE; ZONEID++ ))
do
DESTINATIONZONE=$(($NZONE+$ZONEID))
echo '$!DUPLICATEZONE' >> $OUTMCR
echo "  SOURCEZONE = $ZONEID" >> $OUTMCR
echo "  DESTINATIONZONE = $DESTINATIONZONE" >> $OUTMCR
done

ZONESTART=$(($NZONE+1))
ZONEEND=$(($NZONE+$NZONE))
echo "\$!ALTERDATA  [$ZONESTART-$ZONEEND]" >> $OUTMCR
echo "  EQUATION = $EQUATION" >> $OUTMCR

echo '$!RemoveVar |MFBD|' >> $OUTMCR
exit 0


