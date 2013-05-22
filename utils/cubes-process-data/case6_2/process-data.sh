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

MACROS=("u_scale_y-0.mcr" "u_scale_z-0.5.mcr" "v_scale_z-0.5.mcr");

for MCR in "${MACROS[@]}" 
do 
  tec360 -b -p $MCR
done
