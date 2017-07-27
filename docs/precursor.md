# Concurrent precursor method
LESO is able to use an inflow condition turbulent inflow generated from a
precursor simulation instead of the standard periodic boundary condition. This
allows LESGO uses to simulation developing flow over wind turbine arrays or
arrays of cubes, etc. The concurrent precursor simulation (CPS) module provides
the framework for implementing this inflow condition.

## Overview
Traditionally, when using inflow data from a precursor simulation, it is
required that a complete simulation be conducted before the target simulation
can be performed. During the precursor simulation the inflow data would be
sampled and periodically written to file. Subsequently, this inflow data would
then be read in by the target simulation. While this approach is conceptually
simple, it does have several drawbacks

* Only one simulation is performed at a given time, i.e., must have the data
     from the precursor simulation before the target simulation can be executed.

* May require significant disk space for large simulations.

* Requires significant I/O which may be a hindrance for good computational
     efficiency.

* Must coordinate the precursor simulation with the target simulation to
     ensure enough inflow data for the target simulation.

To alleviate these issues, the CPS module performs the precursor simulation
concurrently with the target simulation. Conceptually the approach is the same,
except of writing the sampled inflow data from the precursor simulation to file
it is copied in memory directly to the target simulation using MPI. Since both
simulations are executed simultaneously, there is no waiting for the precursor
simulation to complete and direct memory copies remove the I/O overhead both in
speed and storage space.

## Implementation

In the CPS module the precursor simulation is conducted in the 'producer' domain
and is called the 'red' domain. The target simulation occurs in the 'consumer'
domain which is labeled the 'blue' domain. MPI communication between these
domains and within themselves is controlled by defining appropriate MPI
communicators.

At the start of the simulation the "MPI_COMM_WORLD" communicator is split into
two local communicators "localComm" where the 'red' and 'blue' each have their
own. lesgo then takes the local communicator and create a communicator called
"comm" using the MPI Cartesian topology functions. The communicator "comm" is
used for all MPI communication within the 'red' or 'blue' domains and is the
intracommuncator for these domains. The "comm" communicator is also used for
standard MPI with CPS turned off, which results in no special treatment for
point-to-point communication when using CPS.

For communication between the 'red' and 'blue' domain an intercommunicator
"interComm" is created. It essentially builds a communication bridge for each
process in the 'red' domain to communicate with the corresponding process in the
'blue' domain.

When using the CPS module, the simulation is executed as multiple program,
multiple data (MPMD) paradigm. Therefore, the simulation is launched with "N"
process and "N/2" get assigned to the 'red' domain and the other "N/2" are
assigned to the 'blue' domain. Within the global "MPI_COMM_WORLD" communicator,
the domains have global ranks assigned to them such that 0 to N/2-1 is given to
the 'red' domain and N/2 to N-1 is assigned to the 'blue' domain. Once the local
communicator "comm" is created each domain has the local ranks 0 to N/2-1
assigned to each of the processes. The intercommunicator then takes these local
ranks and maps rank "n" from 'red' to rank "n" in 'blue' creating the bridge for
copying the inflow data.

There is a capability to shift the domain in order to eliminate the streaks.
This option will shift the domain to the side of the precursor domain.

## Usage

The first step is to build in the CPS support into lesgo by setting

    USE_MPI true
    USE_CPS true

in "CMakeLists.txt". You will need to build in CPS support for both the 'red' and
'blue' domain executables. Other support may be built for required functionality
but at a minimum CPS is needed. An example for this is the setup for developing
flow over a [wind farm](wind.html) using the
[actuator disk model](actuator-disk.html). In the 'red' domain we'd have
standard boundary layer flow so we'd set

    USE_MPI true
    USE_CPS true

as above and build the executable. Then we'd also need to include the wind
turbines module for the 'blue' domain so we'd set

    USE_MPI true
    USE_CPS true
    USE_TURBINES true

and build lesgo.

Once the executables are built you can then setup your cases. For now we'll call
the executable for the 'red' domain "lesgo-red" and the one for the 'blue'
domain "lesgo-blue". Each domain will need it's own run directory so you'll have
to create these; we'll call these directories "red" and "blue". The executables
should then be placed in their respective run directory. A copy of the input
file "lesgo.conf" will have to be place in each of the run directories.

Now the input files have to be configured. The number of processors should be
set what will be used for each domain. So, for example, if each domain will use
4 processes, then

    nproc = 4

in "lesgo.conf" for both cases. When submitting the job, you have to request the
total number of processes being used. Continuing with the example you have to
request 8 process if the 'red' and 'blue' domains use 4 each. The next important
setting is the inflow flag "inflow". For the red domain it must be set to

    inflow = .false.

where the 'blue' domain will use

    inflow = .true.

To provide the inflow condition while numerically still using periodic boundary
conditions, a fringe method is applied. The fringe method is a well established
technique for forcing the velocity field to the desired, sampled field over a
small region called the fringe region. There are several settings which control
the location and size of this fringe region.

One constraint of the CPS module is that the grid spacing of the 'red' and
'blue' domains must match. The domain lengths can be different, but the grid
spacing must be the same to ensure accurate results. Another constraint is that
because these simulation are synchronized, the same value for "dt" must be
specified unless dynamic time stepping is used, then the same CFL value should
be used for both domains.

The specifics on launching the simulation depends on which MPI implementation is
used. This is discussed below, assuming we are using executables named
"lesgo-red" and "lesgo-blue" and run directories "red" and "blue" for the "red"
and "blue" domains, respectively, with a total of "N" processes.

MPICH2 launch command

    mpiexec -f <nodefile> -wdir red -n <N/2> ./lesgo-red : -wdir blue -n <N/2> ./lesgo-blue

When running, all diagnostic information is written to standard out with no
specific order.

## References
Stevens RJAM, Graham J, Meneveau C. "[A concurrent precursor inflow method for Large Eddy Simulations and applications to finite length wind farms](https:/doi.org/10.1016/j.renene.2014.01.024)." *Renewable Energy* **68** (2014). 46-50.