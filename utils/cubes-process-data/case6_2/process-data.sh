##
##  Copyright (C) 2011-2013  Johns Hopkins University
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

MACROS=("u_scale_y-0.mcr" "u_scale_z-0.5.mcr" "v_scale_z-0.5.mcr");

for MCR in "${MACROS[@]}" 
do 
  tec360 -b -p $MCR
done
