#!/bin/bash
#SBATCH --job-name=lesgo
#SBATCH --output=lesgo.out
#SBATCH --ntasks=16
#SBATCH --time=4:00:00

module load mpich/gcc-4.7.2 fftw3/mpich/gcc-4.7.2 hdf5-1.8.15/mpich/gcc-4.7.2 cgns-3.2.1/mpich/gcc-4.7.2

mpirun lesgo-mpi
