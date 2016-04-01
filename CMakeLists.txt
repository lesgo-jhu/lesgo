cmake_minimum_required(VERSION 2.8.9)

################################################################################
#                               User Inputs
################################################################################

# Specify Compiler
# Options: Intel and GNU
set(CMAKE_Fortran_COMPILER_ID "GNU")

# Set the name of the computer you are compiling on.
# Supported computers: panoramix, njord, marcc
set(hostname "njord")

#--Set global VERBOSE flag;
set(VERBOSE false)
#--Sets inclusion of development sections
set(DEVEL false)
# Set extra safety measures 
set(USE_SAFETYMODE true)

# Output additional quantities to file. May be useful
# for diagnosing issues
set(OUTPUT_EXTRA false)

# Use FFTW3 libraries (if no FFTW2 will be used)
set(USE_FFTW3 true)

# Enable parallization support (openmp currently not supported)
set(USE_MPI true)
set(USE_OPENMP false)

# Read and write binary files as big_endian or little_endian
# Options: DEFAULT, LITTLE, or BIG
set(WRITE_ENDIAN "DEFAULT")
set(READ_ENDIAN "DEFAULT")

# Select a turbine model
set(USE_TURBINES false)
set(USE_ATM false)

# Select the data output
set(USE_CGNS false)
set(USE_BINARY true)

# Select level set
set(USE_LVLSET false)
set(USE_RNSE_LS false)
set(USE_CYC_SKEW_LS false)

# Concurrent precursor simulation 
# Requires configuration of two builds: upstream and downstream
set(USE_CPS false)

# Shift the domain in the y (spanwise) direction to eliminate
# streaks effects on averaging
set(USE_STREAKS false)

################################################################################
#                         Set up build (do not change)
################################################################################
# Set compilter versions for Cmake
#  Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(compiler_name "ifort")
#  GNU
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_definitions(-DPPGFOFTRAN)
    set(compiler_name "gfortran")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
find_program(WHICH NAMES which)
exec_program(${WHICH} ARGS ${compiler_name} OUTPUT_VARIABLE compiler_path)
set(CMAKE_Fortran_COMPILER ${compiler_path})

# Project details
enable_language (Fortran)
project (lesgo)
set(exec_name "lesgo")

# Common inlclude and link directories
include_directories($ENV{HOME}/include /usr/local/include /usr/include)
link_directories($ENV{HOME}/lib /usr/local/lib /usr/lib64 $ENV{HOME}/lib /usr/local/lib  /usr/lib64)

# Include host-specific paths
#  panoramix
if (${hostname} STREQUAL "panoramix")
    include_directories(/opt/cgns-3.2.1-hdf5-mpi/include /opt/hdf5-1.8.12-mpi/include /usr/local/hdf5/include /usr/local/cgns/include /act/fftw3/mpich/gcc-4.7.2/include /act/hdf5-1.8.15/mpich/gcc-4.7.2/include /act/cgns-3.3.0/mpich/gcc-4.7.2/include)
    link_directories(/opt/hdf5-1.8.12-mpi/lib  /opt/cgns-3.2.1-hdf5-mpi/lib /usr/local/hdf5/lib  /usr/local/cgns/lib /act/fftw3/mpich/gcc-4.7.2/lib /act/hdf5-1.8.15/mpich/gcc-4.7.2/lib /act/cgns-3.3.0/mpich/gcc-4.7.2/lib)
#  marcc
elseif (${hostname} STREQUAL "marcc")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/act/fftw3/openmpi-1.6/gcc-4.7.2/include /act/hdf5-1.8.14/openmpi-1.6/gcc-4.7.2/include /act/cgns-3.2.1/openmpi-1.6/gcc-4.7.2/include /cm/shared/apps/Intel/fftw/3.3.4/include)
        link_directories(/act/fftw3/openmpi-1.6/gcc-4.7.2/lib /act/hdf5-1.8.14/openmpi-1.6/gcc-4.7.2/lib /act/cgns-3.2.1/openmpi-1.6/gcc-4.7.2/lib /cm/shared/apps/Intel/fftw/3.3.4/lib)
    # GNU -- These have to be installed locally for now.
    else (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif (USE_CGNS)
    endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
#  njord
elseif (${hostname} STREQUAL "njord")
    # nothing yet, but can add if needed
#  all other systems
else (${hostname} STREQUAL "panoramix")
    # nothing yet, but can add if needed
endif (${hostname} STREQUAL "panoramix")

# List common sources
set (Sources cfl_util.f90 clocks.f90 convec.f90 derivatives.f90 divstress_uv.f90 divstress_w.f90 dns_stress.f90 emul_complex.f90 fft.f90 finalize.f90 forcing.f90 fringe_util.f90 functions.f90 grid.f90 ic.f90 ic_dns.f90 initial.f90 initialize.f90 input_util.f90 interpolag_Sdep.f90 interpolag_Ssim.f90 io.f90 lagrange_Sdep.f90 lagrange_Ssim.f90 main.f90 messages.f90 open_file.f90 padd.f90 param.f90 param_output.f90 press_stag_array.f90 ran3.f90 rmsdiv.f90 scaledep_dynamic.f90 sgs_param.f90 sgs_stag_util.f90 sim_param.f90 stat_defs.f90 std_dynamic.f90 string_util.f90 test_filtermodule.f90 tridag_array.f90 types.f90 unpadd.f90 wallstress.f90 wallstress_dns.f90)

if (USE_MPI)
    add_definitions(-DPPMPI)
    set (CMAKE_Fortran_COMPILER mpif90)
    set(Sources ${Sources} mpi_transpose_mod.f90 tridag_array_pipelined.f90 mpi_defs.f90)
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    set(Sources ${Sources} concurrent_precursor.f90)
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_STREAKS)
    add_definitions(-DPPSTREAKS)
    set(exec_name "${exec_name}-streaks")
endif (USE_STREAKS)

if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    set(Sources ${Sources} level_set_base.f90 level_set.f90 linear_simple.f90)
    set(exec_name "${exec_name}-ls")
endif (USE_LVLSET)

if (USE_RNS_LS)
    add_definitions(-DPPRNS_LS)
    set(Sources ${Sources} rns_base_ls.f90 rns_ls.f90)
    set(exec_name "${exec_name}-rns")
    if (USE_CYL_SKEW_LS)
        set(Sources ${Sources} rns_cyl_skew_ls.f90)
    endif (USE_CYL_SKEW_LS)
endif (USE_RNS_LS)

if (USE_CYL_SKEW_LS)
    add_definitions(-DPPCYL_SKEW_LS)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lcgaml")
    set(Sources ${Sources} cyl_skew_base_ls.f90 cyl_skew_ls.f90)
    set(exec_name "${exec_name}-cs")
endif (USE_CYL_SKEW_LS)

if (USE_TURBINES)
    add_definitions(-DPPTURBINES)
    set(Sources ${Sources} turbines.f90 turbines_base.f90)
    set(exec_name "${exec_name}-turbines")
endif (USE_TURBINES)

if (USE_ATM)
    add_definitions(-DPPATM)
    set(Sources ${Sources} atm_base.f90 atm_input_util.f90 actuator_turbine_model.f90 atm_lesgo_interface.f90)
    set(exec_name "${exec_name}-ATM")
endif (USE_ATM)

if (OUTPUT_EXTRA)
    add_definitions(-DPPOUTPUT_EXTRA)
    set(exec_name "${exec_name}-exout")
endif (OUTPUT_EXTRA)

if (USE_DYN_TN)
    add_definitions(-DPPDYN_TN)
    set(exec_name "${exec_name}-dyntn")
endif (USE_DYN_TN)

if (USE_BINARY)
    add_definitions(-DPPBINARY)
    set(exec_name "${exec_name}-binary")
endif (USE_BINARY)

if (USE_SAFETYMODE)
    add_definitions(-DPPSAFETYMODE)
else (USE_SAFETYMODE)
    set(exec_name "${exec_name}-safety_off")
endif (USE_SAFETYMODE)

if (USE_CGNS)
    add_definitions(-DPPCGNS)
    link_libraries(cgns hdf5)    
endif (USE_CGNS)

if (USE_FFTW3)
    add_definitions(-DPPFFTW3)
    link_libraries(fftw3)    
else (USE_FFTW3)
    link_libraries(rfftw fftw)    
endif (USE_FFTW3)

if (VERBOSE)
    add_definitions(-DPPVERBOSE)
endif (VERBOSE)

if (DEVEL)
    add_definitions(-DPPDEVEL)
endif (DEVEL)

if (WRITE_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPWRITE_LITTLE_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "LITTLE")

if (WRITE_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPWRITE_BIG_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "BIG")

if (READ_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPREAD_LITTLE_ENDIAN)
endif (READ_ENDIAN STREQUAL "LITTLE")

if (READ_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPREAD_BIG_ENDIAN)
endif (READ_ENDIAN STREQUAL "BIG")

# Intel
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -funroll-loops -no-wrap-margin -assume byterecl")
# GNU
else (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "-cpp -O2 -ffast-math -funroll-loops")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

add_executable(${exec_name} ${Sources})
