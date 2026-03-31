================================================================                
FFTW many_dft Test Suite
2026-04-01
================================================================

OVERVIEW:
This MWE (Minimum Working Example) is designed to validate the 
"Advanced" (many_dft) interface of FFTW-compatible libraries. 
It specifically tests the handling of:
  - 2D Strided/Embedded memory layouts (nembed)
  - Batch processing (howmany)
  - Spectral differentiation accuracy (Laplacian operator)

SYSTEM REQUIREMENTS:
  - C++17 compliant compiler (g++)
  - GFortran (for interface parity tests)
  - FFTW3 or ArmPL 26.0+

COMPILATION:
The provided Makefile supports standard Fortran, C++ and linker parameters.
By default it the GNU compilers are expected.

  Use the variable FFTW3_DIR to configure the location, i.e.
    make FFTW3_DIR=/path/to/fftw3 all

  To run the tests:
    make run

  To create a source distribution tarball:
    make dist

FILES:
  - Makefile:                 Multi-vendor build system
  - test_suite.cpp:           Core C++ testing logic
  - test_fftw_interface.f90:  Fortran interface check
  - run.sh:                   Build and run tests.
================================================================
