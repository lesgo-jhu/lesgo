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
echo "Enter zone to start duplication:"
read ZONESTART
echo "Enter number of zones to duplicate:"
read NZONE
echo "Enter Tecplot equation to move duplicated data:"
read EQUATION

echo '#!MC 1200' > $OUTMCR
echo '# Created by Tecplot 360 build 12.2.0.9077' >> $OUTMCR
echo "\$!VarSet |MFBD| = '$TECPATH'" >> $OUTMCR


ZONEEND=$(($ZONESTART+$NZONE-1))
for (( ZONEID=$ZONESTART; ZONEID<=$ZONEEND; ZONEID++ ))
do
DESTINATIONZONE=$(($NZONE+$ZONEID))
echo '$!DUPLICATEZONE' >> $OUTMCR
echo "  SOURCEZONE = $ZONEID" >> $OUTMCR
echo "  DESTINATIONZONE = $DESTINATIONZONE" >> $OUTMCR
done

ZONESTART=$(($ZONEEND+1))
ZONEEND=$(($ZONESTART+$NZONE-1))
echo "\$!ALTERDATA  [$ZONESTART-$ZONEEND]" >> $OUTMCR
echo "  EQUATION = '$EQUATION'" >> $OUTMCR

echo '$!RemoveVar |MFBD|' >> $OUTMCR
exit 0


