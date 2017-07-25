# Flow solver
The LES flow solver uses a pseudospectral discretization in the longitudinal and 
spanwise directions (x and y, respectively) and second-order centered finite differencing
in the wall-normal (z) direction. The equations are the filtered N-S equations in the high-Re limit in which molecular
viscosity effects are generally neglected. An eddy viscosity closure model is
applied for the sub-grid scale stresses. Boundary conditions along the x-y
perimeter are periodic. The wall-normal boundary condition options include various wall models which prescribe the surface stresses (for LES), no-slip velocity (for DNS), and stress-free. The default wall-normal configuration is a half-channel (wall model or no-slip at bottom wall, stress-free at the top to mimic a full channel's center plane), but a full-channel can readily be simulated for many (but not all) combinations of wall model and SGS model and DNS. Solid objects can be represented using a level set immersed boundary method. In addition to the bottom wall, solid surfaces represented by the immersed boundary method
also apply a log-law stress boundary condition. For temporal integration, the explicit second-order Adams-Bashforth method is used. Time stepping occurs using either a static time step or CFL-based dynamic time step. Continuity is
preserved by solving the pressure Poisson equation using a direct TDMA solver.

## MPI Domain decomposition

The LES flow solver is parallelized with MPI. The flow domain is evenly divided
in the vertical (z) direction between the MPI processes. As such, the number of vertical gridpoints (nz) should be evenly divisible by the number of processors (nproc). Each MPI-process has
its "own" z-levels indexed by k = 1 to nz-1. The z-levels k = 0 and k = nz are
overlap nodes and are only for holding data that has been copied from the
process below or above the current one. For the main flow solver, the only
overlap data that is required is for the vertical derivatives. Since the
spectral discretization occurs along the slab planes, no additional overlap data
is required. In the future, it may be advantageous to also decompose the slabs
along the y-direction and utilize the parallel FFTW3 MPI library in order to
increase scalability for large domains. For the pressure solver, a pipelining 
technique was chosen to parallelize the TDMA.

## Variable locations

The code employs a staggered grid along the vertical direction. This means that
not all variables are stored on the same grid. We call the two grids the 'uv'
grid and the 'w' grid, where the grids are distinguished based on the location
of the velocity components. The only difference is that the 'uv' grid is shift
+dz/2 in the vertical direction, where dz is the vertical grid spacing. Also,
the 'w' grid conforms to the physical domain such that the first grid point lies
at z=0 and the last grid point at z=L_z where L_z is the domain height.

The following is a list of the variables and the grid they exist on. The names
are the same as found in lesgo.


|  Variables                 | 'uv' grid | 'w' grid  |
| -------------------------- |:---------:|:---------:|
|  u, v, p                   |     X     |           |
|  w                         |           |     X     |
|  dudx, dudy, dvdx, dvdy    |     X     |           |
|  dudz, dvdz                |           |     X     |
|  dwdx, dwdy                |           |     X     |
|  dwdz                      |     X     |           |
|  dpdx, dpdy                |     X     |           |
|  dpdz                      |           |     X     |
|  txx, txy, tyy, tzz        |     X     |           |
|  txz, tyz                  |           |     X     |
|  RHSx, RHSy                |     X     |           |
|  RHSz                      |           |     X     |
|  divtx, divty              |     X     |           |
|  divtz                     |           |     X     |

When writing data to file the variables are interpolated to the same grid.
The statistics may not be on the same grid as the corresponding instantaneous 
variables above. See section below (" A) Time Averaged Data ") for notes on
which grid each statistic is recorded on.

## Data output
The output data from lesgo contains two types: 1) restart data and 2)
visualization data. Both types are described here. It should be noted that the
file names for serial runs will only be listed. In the case of MPI runs some
files will be appended with '.c&lt;id&gt;' where &lt;id&gt; is the z-ordered MPI rank used
in the simulation. For files where this is applicable, it will be noted by the
term 'MPI appendix'

Below are listed the output files along with a description
of their contents. These file use the 'MPI appendix'

  1) vel.out    : Contains the core part of the restart data for lesgo. In this
                  file the velocity, right hand side of the momentum equation, and
                  several SGS variables are store. Essentially all the data
                  required from the previous time step or steps is located in these
                  files.

  2) tavg.out   : Contains the time averaged data from the previous simulation. In
                  this file the running averages of the velocity, Reynolds stress,
                  etc are store here and are used to maintain the averaging
                  between simulations.

The visualization data is located in the directory "output". The "output"
directory is created by making a 'mkdir' system call and should work fine on all
Unix based systems. By default these file are written in binary formatted files
using direct fortran write calls. The output to this directory depends on the
settings in the OUTPUT block in the input file "lesgo.conf". In this section we
will discuss the output from the core LES solver of lesgo. All data output from
various modules will be discussed in their respective sections.

For large data files, two binary formats can be used: 1) lesgo binary format and 2) CGNS.
These files are denoted below with the dummy file extension '.ext'. For lesgo binary, the
extension will actually be '.bin' and for CGNS the actual extension with '.cgns'. For 
smaller data files, ASCII is used and end with '.dat'.

A) Time averaged data on entire domain (these file use the 'MPI appendix'). 
Note the comments in brackets [] which denote the grid that each statistic is 
recorded on.

  * veluv_avg.ext : Mean u, v, w velocity field [on uv grid]
  * velw_avg.ext  : Mean w velocity [on w grid]
  * vel2_avg.ext  : Mean of square products of the velocity field
                  : [all on uv grid except for ww which is on w grid]
  * force_avg.ext : Mean field of modeled forces (IBM, turbines, RNS, etc)
  * tau_avg.ext   : Mean of the sub-grid scale stress tensor
                  : [all on uv grid except txz, tyz which are on w grid]
  * rs.ext        : Reynolds stresses [same grids as corresponding vel2]
  * cs_opt2.ext   : Square of Smagorinsky coefficient from SGS modeling

B) Domain Data (these files use the 'MPI appendix')

  * vel.&lt;timestep&gt;.ext      : Instantaneous velocity field at time step &lt;timestep&gt;

C) Sampled Point Data

  * vel.x-&lt;xloc&gt;.y-&lt;yloc&gt;.z-&lt;zloc&gt;.dat :
    Instantaneous velocity sampled at point (&lt;xloc&gt;,&lt;yloc&gt;,&lt;zloc&gt;)

D) Data along x-planes (these files use the 'MPI appendix')

  * vel.x-&lt;xloc&gt;.&lt;timestep&gt;.ext :
    Instantaneous velocity sampled at x location &lt;xloc&gt; for the time step &lt;timestep&gt;

E) Data along y-planes (these files use the 'MPI appendix')

  * vel.y-&lt;yloc&gt;.&lt;timestep&gt;.ext :
    Instantaneous velocity sampled at y location &lt;yloc&gt; for the time step &lt;timestep&gt;

F) Data along z-planes

  * vel.z-&lt;zloc&gt;.&lt;timestep&gt;.ext : 
    Instantaneous velocity sampled at z location &lt;zloc&gt; for the time step &lt;timestep&gt;
