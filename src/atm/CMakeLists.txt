add_library(atm atm_base.f90 atm_input_util.f90 actuator_turbine_model.f90
    atm_lesgo_interface.f90)
target_link_libraries(atm core)