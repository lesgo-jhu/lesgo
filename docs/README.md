LESGO solves the filtered Navier-Stokes equations in the high-Reynolds number
limit on a Cartesian mesh. Originally designed to simulate flow in the
atmospheric boundary layer, LESGO has been extended and used to simulate flow
over tree canopies, wall-mounted cubes, and wind turbine arrays, among other things.
At its core is the LES flow solver. Built on top of the solver are modules that provide
additional functionality such as immersed boundary methods, wind farm modeling, and so on.

LESGO was originally based on the code presented in John D. Albertson's PhD Thesis "Large
eddy simulation of land-atmosphere interaction" (University of California, Davis
1996). In the intervening years, many researchers have [contributed](contributors.html)
to LESGO's code base and used LESGO in dozens of scientific [publications](publications.html).

If you use LESGO for any purpose, please [cite](citing.html) appropriately.

LESGO is distributed without any warranty or technical support.

## Features
* [Flow solver](solver.html)
* [Subgrid scale models](subgrid.html)
* [Wall models](wall-model.html)
* [Level set immersed boundary method](levelset.html)
* [Wind turbine modeling](wind.html)
* [Concurrent precursor simulations](precursor.html)

## Getting started
To [get started](start.html), download the code using the buttons on the left.
You'll need a modern Fortran compiler that supports C preprocessor directives,
[CMake](https://cmake.org/), and [FFTW3](http://www.fftw.org/).

## Licensing
LESGO is a free, open-source tool published under the
[GNU General Public License Version 3](http://www.gnu.org/licenses/)

## Acknowledgements
Development of LESGO has been supported in part by the National Science Foundation.