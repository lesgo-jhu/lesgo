#!/bin/bash

# For launching lesgo jobs on njord cluster.
# Searches lesgo configuration file for number of processors to be used then launches mpi job.
# Author: Joel Bretheim, jbretheim@gmail.com

conf=lesgo.conf

nproc=`grep "nproc" $conf | awk '{print $3}'`

printf "Found nproc = $nproc in $conf \n"

mpirun -np $nproc ./lesgo-mpi-* > lesgo.out &
