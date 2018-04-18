cmake_minimum_required(VERSION 2.8.9)

################################################################################
#                               User Inputs
################################################################################

# Specify Compiler
# Options: Intel and GNU
set(CMAKE_Fortran_COMPILER_ID "GNU")

# Read and write binary files as big_endian or little_endian
# Options: DEFAULT, LITTLE, or BIG
set(WRITE_ENDIAN "DEFAULT")
set(READ_ENDIAN "DEFAULT")

# Options to enable various features in the code
option(USE_MPI "Enable parallelization support." ON)
option(USE_CPS "Enable concurrent precursor simulations." OFF)
option(USE_HIT "Enable homogeneous isotropic turbulence input." OFF)
option(USE_STREAKS "Shift the domain in the y (spanwise) direction to eliminate streaks effects on averaging." OFF)
option(USE_LVLSET "Enable level set." OFF)
option(USE_TURBINES "Enable actuator disk model without rotation." OFF)
option(USE_ATM "Enable actuator line/section model" OFF)
option(USE_DYN_TN "Dynamically update T, the timescale for Lagrangian averaging, using Taylor timescale." OFF)
option(USE_SAFETYMODE "Set extra safety measures." ON)
option(USE_CGNS "Enable CGNS data output." OFF)

################################################################################
#                         Set up build (do not change)
################################################################################
# Set compiler versions for Cmake
#  Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(compiler_name "ifort")
#  GNU
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(compiler_name "gfortran")
endif ()
find_program(WHICH NAMES which)
exec_program(${WHICH} ARGS ${compiler_name} OUTPUT_VARIABLE compiler_path)
set(CMAKE_Fortran_COMPILER ${compiler_path})

# Project details
enable_language (Fortran)
project (lesgo)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

# Search paths
set(include_search_paths /cm/shared/apps/Intel/fftw/3.3.4/include
    /opt/cgns-3.2.1-hdf5-mpi/include /opt/hdf5-1.8.12-mpi/include
    /usr/local/hdf5/include /usr/local/cgns/include
    /act/fftw3/mpich/gcc-4.7.2/include /act/hdf5-1.8.15/mpich/gcc-4.7.2/include
    /act/cgns-3.3.0/mpich/gcc-4.7.2/include $ENV{HOME}/cgns/include/)
set(library_search_paths /cm/shared/apps/Intel/fftw/3.3.4/lib
    /opt/hdf5-1.8.12-mpi/lib  /opt/cgns-3.2.1-hdf5-mpi/lib /usr/local/hdf5/lib
    /usr/local/cgns/lib /act/fftw3/mpich/gcc-4.7.2/lib
    /act/hdf5-1.8.15/mpich/gcc-4.7.2/lib /act/cgns-3.3.0/mpich/gcc-4.7.2/lib
    $ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)

# FFTW
find_path(fftw_include fftw3.f ${include_search_paths})
include_directories(${fftw_include})
find_library(fftw3_library fftw3 ${library_search_paths})
link_libraries(${fftw3_library})

# CGNS
if (USE_CGNS)
    find_path(cgns_include cgns.mod ${include_search_paths})
    include_directories(${cgns_include})
    find_library(cgns_library cgns ${library_search_paths})
    link_libraries(${cgns_library})
    find_library(hdf5_library hdf5 ${library_search_paths})
    link_libraries(${hdf5_library})
endif (USE_CGNS)

# MPI
if (USE_MPI)
    set (CMAKE_Fortran_COMPILER mpif90)
endif (USE_MPI)

# Compiler flags
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # Intel
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -funroll-loops -no-wrap-margin -assume byterecl")
else ()
    # GNU
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -ffast-math -funroll-loops -fall-intrinsics -Wall -fcheck=all")
endif ()

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
add_subdirectory(src)
