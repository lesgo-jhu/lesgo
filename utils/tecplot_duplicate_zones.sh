##
##  Copyright 2010,2011 Johns Hopkins University
##
##  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
##  use this file except in compliance with the License. You may obtain a copy of
##  the License at:
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
##  Unless required by applicable law or agreed to in writing, software 
##  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
##  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
##  License for the specific language governing permissions and limitations under
##  the License.
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


