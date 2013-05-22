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
