#!/bin/bash

# For launching lesgo jobs on njord cluster.
# Author: Joel Bretheim, jbretheim@gmail.com

# specify the configuration file
conf=lesgo.conf

# grab number of processors from conf
nproc=`grep "nproc" $conf | awk '{print $3}'`
printf "Found nproc = $nproc in $conf \n"

# launch job
mpirun -np $nproc ./lesgo-mpi* > lesgo.out &
