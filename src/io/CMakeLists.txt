add_library(io io.f90 time_average.f90 data_writer.f90)
target_link_libraries(io core sgs)