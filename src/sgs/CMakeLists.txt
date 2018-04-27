add_library(sgs interpolag_Ssim.f90 lagrange_Ssim.f90 sgs_param.f90
    std_dynamic.f90 interpolag_Sdep.f90  lagrange_Sdep.f90 scaledep_dynamic.f90
    sgs_stag_util.f90)
target_link_libraries(sgs core)