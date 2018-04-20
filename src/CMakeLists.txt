# Executable name
set(exec_name "lesgo")

# Common sources
set (Sources cfl_util.f90 clocks.f90 convec.f90 derivatives.f90 data_writer.f90
     divstress_uv.f90 divstress_w.f90 emul_complex.f90 fft.f90 finalize.f90
     forcing.f90 fringe_util.f90 functions.f90 grid.f90
     initial.f90 initialize.f90 input_util.f90 interpolag_Sdep.f90
     interpolag_Ssim.f90 io.f90 lagrange_Sdep.f90 lagrange_Ssim.f90 main.f90
     messages.f90 param.f90 param_output.f90
     press_stag_array.f90 rmsdiv.f90 scaledep_dynamic.f90
     sgs_param.f90 sgs_stag_util.f90 sim_param.f90 stat_defs.f90
     std_dynamic.f90 string_util.f90 test_filtermodule.f90 tridag_array.f90
     time_average.f90 types.f90 wallstress.f90 init_random_seed.f90
     iwmles.f90)

if (USE_MPI)
    add_definitions(-DPPMPI)
    set(Sources ${Sources} mpi_transpose_mod.f90 mpi_defs.f90)
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    set(Sources ${Sources} concurrent_precursor.f90)
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_HIT)
    add_definitions(-DPPHIT)
    set(Sources ${Sources} hit_inflow.f90)
    set(exec_name "${exec_name}-HIT")
endif (USE_HIT)

if (USE_STREAKS)
    add_definitions(-DPPSTREAKS)
    set(exec_name "${exec_name}-streaks")
endif (USE_STREAKS)

if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    set(Sources ${Sources} level_set_base.f90 level_set.f90 linear_simple.f90
        trees_pre_ls.f90 trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90
        trees_global_fmask_ls.f90)
    set(exec_name "${exec_name}-ls")
endif (USE_LVLSET)

if (USE_TURBINES)
    add_definitions(-DPPTURBINES)
    set(Sources ${Sources} turbines.f90 turbine_indicator.f90)
    set(exec_name "${exec_name}-turbines")
endif (USE_TURBINES)

if (USE_ATM)
    add_definitions(-DPPATM)
    set(Sources ${Sources} atm_base.f90 atm_input_util.f90
        actuator_turbine_model.f90 atm_lesgo_interface.f90)
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

# main lesgo executable
add_executable(${exec_name} ${Sources})

