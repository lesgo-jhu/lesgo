CGNS and HDF5 may be used for IO in LESGO, but they must be installed 
separately to take advantage of these features.
Below are instructions for installing these libraries 
(based on CGNS 3.3.0)

* HDF5 library [www.hdfgroup.org/HDF5/release/obtain5.html](http://www.hdfgroup.org/HDF5/release/obtain5.html)
* CGNS library [cgns.sourceforge.net](http://cgns.sourceforge.net/)

It is required to build with Fortran and parallel IO support.

The following bash commands should install HDF5 from the SVN 
repository and CGNS from the git repository. The install locations 
need to be specified in the Makefile. To build these libraries you 
will need gfortran, gcc, zlib, an mpi version(tested with openmpi 
and mvapich2)

## Bash script
The following bash script has been used to successfully install these libraries. However, it may have to be changed
for different systems of library versions.

    #!/bin/bash

    # Location to install hdf5
    export HDF5INSTALL=$HOME/hdf5
    # Location to install CGNS
    export CGNSINSTALL=$HOME/cgns
    # MPI compiler mpicc
    export MPICC=mpicc

    # HDF5
    svn co -q https://svn.hdfgroup.uiuc.edu/hdf5/branches/hdf5_1_8
    cd hdf5_1_8 
    CC=$MPICC ./configure --enable-fortran --enable-parallel --enable-shared --prefix=$HDF5INSTALL
    make
    make install

    # CGNS
    git clone -b master https://github.com/CGNS/CGNS.git CGNS
    cd CGNS/src    
    export FC=gfortran
    export F77=gfortran
    export CC=mpicc
    export FLIBS="-Wl,--no-as-needed -ldl -lz"
    export LIBS="-Wl,--no-as-needed -ldl -lz"
    ./configure --verbose \
    --enable-parallel \
    --with-fortran \
    --with-hdf5=$HDF5INSTALL \
    --with-mpi \
    --enable-64bit \
    --enable-lfs \
    --disable-cgnstools \
    --disable-x \
    --prefix=$CGNSINSTALL \
    --libdir=$CGNSINSTALL/lib \
    --includedir=$CGNSINSTALL/include    
    make
    make install
