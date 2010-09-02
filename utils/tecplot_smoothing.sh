#!/bin/bash
echo "Enter output macro filename:"
read OUTMCR
#echo "Enter TECPATH string:"
#read TECPATH
TECPATH=$(pwd)
echo "Enter number of zones to smooth:"
read NZONE
echo "Enter variable ID to smooth:"
read VARID
echo "Enter number of smoothing passes (enter 3 for default):"
read SPASS

echo '#!MC 1200' > $OUTMCR
echo '# Created by Tecplot 360 build 12.2.0.9077' >> $OUTMCR
echo "\$!VarSet |MFBD| = '$TECPATH'" >> $OUTMCR

for (( ZONEID=1; ZONEID<=$NZONE; ZONEID++ ))
do
echo '$!SMOOTH' >> $OUTMCR
echo "  ZONE = $ZONEID" >> $OUTMCR
echo "  VAR = $VARID" >> $OUTMCR
echo "  NUMSMOOTHPASSES = $SPASS" >> $OUTMCR
echo '  SMOOTHWEIGHT = 0.5' >> $OUTMCR
echo '  SMOOTHBNDRYCOND = FIXED' >> $OUTMCR	
done



echo '$!RemoveVar |MFBD|' >> $OUTMCR
exit 0
