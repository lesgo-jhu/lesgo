# Executable name
set(exec_name "lesgo")
set(libs core io)

# Common sources
set (sources clocks.f90 convec.f90 derivatives.f90 divstress_uv.f90
    divstress_w.f90 finalize.f90 forcing.f90 fringe_util.f90 functions.f90
    initialize.f90 initial.f90 main.f90 param_output.f90 press_stag_array.f90
    rmsdiv.f90 tridag_array.f90 wallstress.f90 init_random_seed.f90)

set (core_sources param.f90 sim_param.f90 messages.f90 string_util.f90
    test_filtermodule.f90 fft.f90 emul_complex.f90 functions.f90
    grid.f90 iwmles.f90 cfl_util.f90 input_util.f90 )

if (USE_MPI)
    add_definitions(-DPPMPI)
    set(core_sources ${core_sources} mpi_defs.f90)
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    set(sources ${sources} concurrent_precursor.f90)
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_HIT)
    add_definitions(-DPPHIT)
    set(sources ${sources} hit_inflow.f90)
    set(exec_name "${exec_name}-HIT")
endif (USE_HIT)

if (USE_STREAKS)
    add_definitions(-DPPSTREAKS)
    set(exec_name "${exec_name}-streaks")
endif (USE_STREAKS)

if (OUTPUT_EXTRA)
    add_definitions(-DPPOUTPUT_EXTRA)
    set(exec_name "${exec_name}-exout")
endif (OUTPUT_EXTRA)

if (USE_DYN_TN)
    add_definitions(-DPPDYN_TN)
    set(exec_name "${exec_name}-dyntn")
endif (USE_DYN_TN)

if (USE_SAFETYMODE)
    add_definitions(-DPPSAFETYMODE)
else (USE_SAFETYMODE)
    set(exec_name "${exec_name}-safety_off")
endif (USE_SAFETYMODE)

if (USE_CGNS)
    add_definitions(-DPPCGNS)
    set(exec_name "${exec_name}-cgns")
endif (USE_CGNS)

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

# Basic libraries
add_subdirectory(sgs)
add_subdirectory(io)

# Actuator disk model
if (USE_TURBINES)
    add_definitions(-DPPTURBINES)
    add_subdirectory(turbines)
    set(exec_name "${exec_name}-turbines")
    set(libs ${libs} turbines)
endif (USE_TURBINES)

# Actuator line model
if (USE_ATM)
    add_definitions(-DPPATM)
    add_subdirectory(atm)
    set(exec_name "${exec_name}-ATM")
    set(libs ${libs} atm)
endif (USE_ATM)

# Level set
if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    add_subdirectory(level_set)
    set(exec_name "${exec_name}-ls")
    set(libs ${libs} level_set)
endif (USE_LVLSET)

# mmin lesgo executable
add_library(core ${core_sources})
add_executable(${exec_name} ${sources})
target_link_libraries(${exec_name} ${libs})
