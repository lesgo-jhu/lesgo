add_library(level_set level_set_base.f90 level_set.f90 linear_simple.f90
        trees_pre_ls.f90 trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90
        trees_global_fmask_ls.f90)
target_link_libraries(level_set core)