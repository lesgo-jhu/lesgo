#!/bin/bash
# This script creates a Tecplot file with all the zplane_avg variables together in one zone

# Parameters to set:
    # Number of points in the z-direction (entire domain)
        nz=32  
    # Location of zplane_avg files
        file_path='../output'
    # Location of sh/temp_ files
        temp_path='.'


# Create the Tecplot header for zplane_avg_combined.dat
    cp ${temp_path}/temp_zplane_avg_combined.dat ${file_path}/zplane_avg_combined.dat
    rpl -w temp_nz ${nz} ${file_path}/zplane_avg_combined.dat

# Remove the headers from the data files -- store as temp_*
    # All but tau needs to remove z variable as well -- and order matters
    line_start=$((nz+4))
    more +4 ${file_path}/tau_zplane_avg.dat > ${temp_path}/temp_tau_zplane_avg.dat
    more +${line_start} ${file_path}/ddz_zplane_avg.dat > ${temp_path}/temp_ddz_zplane_avg.dat
    more +${line_start} ${file_path}/force_zplane_avg.dat > ${temp_path}/temp_force_zplane_avg.dat    
    more +${line_start} ${file_path}/vel2_zplane_avg.dat > ${temp_path}/temp_vel2_zplane_avg.dat
    more +${line_start} ${file_path}/vel_zplane_avg.dat > ${temp_path}/temp_vel_zplane_avg.dat
    more +${line_start} ${file_path}/cs_opt2_zplane.dat > ${temp_path}/temp_cs_opt2_zplane.dat
    more +${line_start} ${file_path}/rs_zplane.dat > ${temp_path}/temp_rs_zplane.dat    

# Combine the files
    cat ${file_path}/zplane_avg_combined.dat ${temp_path}/temp_{tau,ddz,force,vel2,vel}_zplane_avg.dat ${temp_path}/temp_{cs_opt2,rs}_zplane.dat > ${file_path}/combined.dat

# Clean up
    cp ${file_path}/combined.dat ${file_path}/zplane_avg_combined.dat
    rm ${file_path}/combined.dat
    rm ${temp_path}/temp_{tau,ddz,force,vel2,vel}_zplane_avg.dat 
    rm ${temp_path}/temp_{cs_opt2,rs}_zplane.dat
    
