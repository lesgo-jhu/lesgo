##
##  Copyright (C) 2010-2013  Johns Hopkins University
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
