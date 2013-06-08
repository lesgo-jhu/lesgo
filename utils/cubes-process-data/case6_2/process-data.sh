#!/bin/bash

MACROS=("u_scale_y-0.mcr" "u_scale_z-0.5.mcr" "v_scale_z-0.5.mcr");

for MCR in "${MACROS[@]}" 
do 
  tec360 -b -p $MCR
done
