# Installing on MacOS with Homebrew

LESGO comes with support for MacOS with libraries and dependencies installed 
using [Homebrew](https://brew.sh/). If some of these libaries or dependencies
are already installed and you do not want to overwrite these using the Homebrew
versions, some modification of the CMakeLists.txt file will be 
necessary.

1. Install homebrew from [https://brew.sh/](https://brew.sh/)
2. Install gcc, open-mpi, fftw, and cmake in the terminal:
    
        brew install gcc
        brew install open-mpi
        brew install fftw
        brew install cmake

3. Change the environment flag for lesgo on line 14 of CMakeLists.txt to:

        set(hostname "darwin-brew")

4. Compile in terminal:

        ./build-lesgo

