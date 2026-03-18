cmake_minimum_required(VERSION 3.21)

project(fdlbm
    VERSION 0.1
    LANGUAGES Fortran
)

option(DUGKS "Enable the Discrete Unified Gas Kinetic Scheme" OFF)

add_library(fdlbm
    src/precision.F90
    src/lattice.F90
    src/periodic_lw.F90
    src/periodic_lw_fem.f90
    src/periodic_dugks.F90
    src/taylor_green.f90
    src/gnuplot_io.F90
    )

if(DUGKS)
    target_compile_definitions(fdlbm PUBLIC DUGKS)
endif()

# Create the executable target
add_executable(main_tgv src/main_taylor_green.F90)
target_link_libraries(main_tgv PRIVATE fdlbm)