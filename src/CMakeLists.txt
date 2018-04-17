# Common sources
set (Sources cfl_util.f90 clocks.f90 convec.f90 derivatives.f90 divstress_uv.f90
     divstress_w.f90 emul_complex.f90 fft.f90 finalize.f90
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
    set(Sources ${Sources} mpi_transpose_mod.f90 mpi_defs.f90)
endif (USE_MPI)

if (USE_CPS)
    set(Sources ${Sources} concurrent_precursor.f90)
endif (USE_CPS)

if (USE_HIT)
    set(Sources ${Sources} hit_inflow.f90)
endif (USE_HIT)

if (USE_LVLSET)
    set(Sources ${Sources} level_set_base.f90 level_set.f90 linear_simple.f90
        trees_pre_ls.f90 trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90
        trees_global_fmask_ls.f90)
endif (USE_LVLSET)

if (USE_TURBINES)
    set(Sources ${Sources} turbines.f90 turbine_indicator.f90)
endif (USE_TURBINES)

if (USE_ATM)
    set(Sources ${Sources} atm_base.f90 atm_input_util.f90
        actuator_turbine_model.f90 atm_lesgo_interface.f90)
endif (USE_ATM)

# main lesgo executable
add_executable(${exec_name} ${Sources})
