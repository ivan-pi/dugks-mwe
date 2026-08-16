cmake_minimum_required(VERSION 3.21)

project(fdlbm
    VERSION 0.1
    LANGUAGES Fortran
)

option(DUGKS "Enable the Discrete Unified Gas Kinetic Scheme" OFF)
option(WITH_FFTW "Enable periodic mass solver" OFF)

#find_package(OpenMP COMPONENTS Fortran)

add_library(fdlbm
    src/precision.F90
    src/lattice.F90
    src/periodic_lw.F90
    src/periodic_lw_fem.F90
    src/periodic_dugks.F90
    src/taylor_green.f90
    src/gnuplot_io.F90
)

if(WITH_FFTW)
    target_sources(fdlbm PRIVATE src/periodic_solver.F90)
    target_link_libraries(fdlbm PUBLIC armpl)
    target_link_directories(fdlbm PUBLIC
        /opt/arm/armpl_26.01_flang-21/lib/)
    #    /opt/homebrew/Cellar/fftw/3.3.10_3/lib
    #    /opt/homebrew/Cellar/libomp/22.1.1/lib)
    target_include_directories(fdlbm PRIVATE
        /opt/arm/armpl_26.01_flang-21/include/)
        #/opt/homebrew/Cellar/fftw/3.3.10_3/include)
endif()

if(DUGKS)
    target_compile_definitions(fdlbm PUBLIC DUGKS)
endif()

# Create the executable target
add_executable(main_tgv src/main_taylor_green.F90)
target_link_libraries(main_tgv PRIVATE fdlbm)

if(WITH_FFTW)
    target_compile_definitions(main_tgv PRIVATE WITH_FFTW)
endif()

if(WITH_FFTW)
    add_executable(test_periodic_solver src/test_periodic_solver.f90 src/periodic_solver.F90)
    target_link_libraries(test_periodic_solver PRIVATE fftw3)
    target_link_directories(test_periodic_solver PRIVATE /opt/homebrew/Cellar/fftw/3.3.10_3/lib)
    target_include_directories(test_periodic_solver PRIVATE /opt/homebrew/Cellar/fftw/3.3.10_3/include)
endif()