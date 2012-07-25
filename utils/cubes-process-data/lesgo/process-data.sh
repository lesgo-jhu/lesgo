#!/bin/bash

# Set flag to merge MPI vel_avg.dat.c* files to single vel_avg.dat file
MPIMERGE=1
NPROC=4

# Distance to subtract to zero cube
XOFFSET=1.5625; YOFFSET=2.0625

# Spatial dimensions as set in the computational domain
Nx=64; XMIN=0; XMAX=8.0
Ny=64; YMIN=0; YMAX=8.0
Nz=29; ZMIN=0.0625; ZMAX=3.5625

# Velocity spatial scale (DON'T TOUCH)
Us=0.4
Vs=1.2

# Perform MPI merge
if [ "$MPIMERGE" == "1" ]; then
  MCR=vel_avg.mcr
  cp $MCR.rpl $MCR
  rpl -- '<NX>' $Nx $MCR
  rpl -- '<XMIN>' $XMIN $MCR
  rpl -- '<XMAX>' $XMAX $MCR
  rpl -- '<NY>' $Ny $MCR
  rpl -- '<YMIN>' $YMIN $MCR
  rpl -- '<YMAX>' $YMAX $MCR
  rpl -- '<NZ>' $Nz $MCR
  rpl -- '<ZMIN>' $ZMIN $MCR
  rpl -- '<ZMAX>' $ZMAX $MCR
  
  ./mpi-merge.sh $NPROC

fi

# Compute reference velocity
cp -v vel_ref.mcr.rpl vel_ref.mcr
rpl -- '<XOFFSET>' $XOFFSET vel_ref.mcr
rpl -- '<YOFFSET>' $YOFFSET vel_ref.mcr
tec360 -b -p vel_ref.mcr

# Extract reference velocity
Ur=`tail -1 vel_ref.dat`
echo
echo Reference velocity : $Ur
echo

# Replace reference velocity
cp -v vel_extract.mcr.rpl vel_extract.mcr
rpl -- '<XOFFSET>' $XOFFSET vel_extract.mcr
rpl -- '<YOFFSET>' $YOFFSET vel_extract.mcr
rpl -- '<NY>' $Ny vel_extract.mcr
rpl -- '<YMIN>' `echo "$YMIN-$YOFFSET" | bc` vel_extract.mcr
rpl -- '<YMAX>' `echo "$YMAX-$YOFFSET" | bc` vel_extract.mcr
rpl -- '<NZ>' $Nz vel_extract.mcr
rpl -- '<ZMIN>' $ZMIN vel_extract.mcr
rpl -- '<ZMAX>' $ZMAX vel_extract.mcr

# Extract velocity profiles
tec360 -b -p vel_extract.mcr

# Replace velocity scales
cp -v vel_scale.mcr.rpl vel_scale.mcr
rpl -- '<US>' $Us vel_scale.mcr
rpl -- '<VS>' $Vs vel_scale.mcr
rpl -- '<UR>' $Ur vel_scale.mcr

# Scale and displace velocity profiles
tec360 -b -p vel_scale.mcr

# Compute total forces
cp -v force_tot.mcr.rpl force_tot.mcr
if [ "$MPIMERGE" == "1" ]; then
  ./force-tot.sh $NPROC `echo "$XMAX-$XMIN" | bc` $Nx
else
  ./force-tot.sh 1 `echo "$XMAX-$XMIN" | bc` $Nx
fi
