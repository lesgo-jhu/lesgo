##
##  Copyright 2011 Johns Hopkins University
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
######################################################################
# 
# This script generates Tecplot macros used for creating RNS reference
# regions. This script currently works with reference regions that are
# defined by a x-plane defined by 3 points. The 3 points are read from
# the output of lesgo (typically lesgo.out) which are contained in
# tables with the headers labeled as:
# 
# | ID | NPOINT | AREA | P1.X | P1.Y | P1.Z | P2.X | P2.Y | P2.Z | P3.X | P3.Y | P3.Z |
#
# There is a table for each of the r, \beta, and b-elements. The data
# in the table defines a plane. To get the bounding box P1 is used as
# the starting point and P3+(dX,0,0) is used as the second. The value
# for dX is computed using "alpha_dist" from lesgo_param.out, and P2.Z
# and P3.Z as:
#
#   dX = 2 * alph_dist * ( P3.Z - P2.Z )
# 
# Also note that the points in the x-plane are arranged as:
#
#                   P3
#                   |
#                   |
#                   |
#      P1-----------P2
#
######################################################################
ELEMENT=( "R_ELEM" "BETA_ELEM" "B_ELEM" )
#SECTION=( "Filling R_ELEM Reference Plane Arrays" \
#          "Filling BETA_ELEM Reference Plane Arrays" \
#          "Filling B_ELEM Reference Plane Arrays" )

declare -i istart
declare -i istop
declare -i i
declare -i data_read

# Grab the factor of the box height that the first plane is
# upstream of the center of the box
alpha_dist=$(grep alpha_dist lesgo_param.out | cut -f2 -d:)

echo "Global values:"
echo "  alpha_dist = $alpha_dist"
echo ""

function print_macro_header {
    echo "#!MC 1300"
    echo "# Created by Tecplot 360 build 13.1.0.15185"
    echo "\$!VarSet |MFBD| = '.'"
}

function print_tecplot_box {

    local _p1=( "$1" "$2" "$3" )
    local _p2=( "$4" "$5" "$6" )

    echo "\$!CREATERECTANGULARZONE"
    echo "  IMAX = 2"
    echo "  JMAX = 2"
    echo "  KMAX = 2"
    echo "  X1 = ${_p1[0]}"
    echo "  Y1 = ${_p1[1]}"
    echo "  Z1 = ${_p1[2]}"
    echo "  X2 = ${_p2[0]}"
    echo "  Y2 = ${_p2[1]}"
    echo "  Z2 = ${_p2[2]}"
    echo "  XVAR = 1"
    echo "  YVAR = 2"
    echo "  ZVAR = 3"
}

function print_macro_tailer {
    echo "\$!RemoveVar |MFBD|"
}

for elem in "${ELEMENT[@]}"; 
do

    outfile=$(echo "${elem}-ref-regions.mcr" | tr '[:upper:]' '[:lower:]')
    echo "Creating $outfile:"    

    print_macro_header > "$outfile"

    HEADER="Filling $elem Reference Plane Arrays"

    istart=0
    # Get the line that the section header exists
    istart=$(grep -n "$HEADER" lesgo.out | cut -f1 -d:)

    istop=0
    i=$istart
    read_data=0

    while [ $istop -eq 0 ]; 
    do
        i=$i+1

        # Read the line into an array
        line=(`sed -n "${i}p" lesgo.out`)
        # Check if the first line element is an integer
        if [[ "${line[0]}" =~ ^([0-9]+)$ ]]; then
         
            echo "  creating region ${line[0]}"
          
            [[ $read_data -eq 0 ]] && read_data=1

            P1=( ${line[3]} ${line[4]} ${line[5]} )
            P2=( ${line[6]} ${line[7]} ${line[8]} )
            P3=( ${line[9]} ${line[10]} ${line[11]} )

            # Compute offset for second bounding point
            P3[0]=$(echo "${P3[0]} + 2*${alpha_dist}*(${P3[2]}-${P2[2]})" | bc)
            
            # Create reference region with specified bounding points
            print_tecplot_box "${P1[@]}" "${P3[@]}" >> "$outfile"

        else

            [[ $read_data -eq 1 ]] && istop=1
                
        fi
    done
    
    print_macro_tailer >> "$outfile"

done