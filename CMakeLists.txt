cmake_minimum_required(VERSION 2.8.9)

# Project details
enable_language (Fortran)
project (wake_model)

# Library links
link_libraries(fftw3 blas lapack)

# Common include and link directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
include_directories(/usr/include)
link_directories(/usr/lib)

# List common sourcesi
set (Sources types.f90 wake_model.f90 wake_model_base.f90 param.f90 util.f90 messages.f90 string_util.f90 open_file.f90 init_random_seed.f90 functions.f90 grid.f90 sim_param.f90 cubic_spline.f90 bicubic_spline.f90 tridiagonal.f90 turbines.f90 stat_defs.f90 minimize.f90 conjugate_gradient.f90 line_search.f90 lbfgsb.f90 pchip.f90 bi_pchip.f90 wake_model_estimator.f90 turbines_mpc.f90 wake_model_adjoint.f90 input_util.f90)
#set (CMAKE_Fortran_FLAGS "-cpp -O3 -DPPTURBINES -ffast-math -funroll-loops -fall-intrinsics -Wall")

set (CMAKE_Fortran_FLAGS "-cpp -O0 -g -fbacktrace -ffpe-trap=invalid,zero -DPPTURBINES -fall-intrinsics -Wall")

#set (CMAKE_Fortran_FLAGS "-cpp -O3 -DPPTURBINES -ffast-math -funroll-loops -fall-intrinsics -fexternal-blas -Wall")

# main lesgo executable
add_executable(mpc mpc.f90 ${Sources})
add_executable(startup startup.f90  ${Sources})
