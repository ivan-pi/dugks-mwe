# dugks-mwe

Preprocessor flags:
* `-DSPLIT` - use split collision kernel (no offloading available)
* `-DDUGKS` - use DUGKS mode (default is an FVM scheme)
* `-DBGK_OFFLOAD` - use target offload in the streaming kernel
* `-DWITH_SP` - use single precision (32-bit)

Environment variables:
* `N` - grid size
* OpenMP variables (`OPENMP_NUM_THREADS`)

To-Do:
* Implement target offload in collision step (`collision_bgk.f90`)
* OpenMP data regions (`lattice.F90`)

Build steps:
```
$ make FC=gfortran FCFLAGS="-fopenmp -O3 -march=native"
```
The flags `LDFLAGS` and `LDLIBS` can be modified if needed.

Running the program:
```
$ N=120 ./main_taylor_green <CFL>
```
`CFL` is the Courant number (a number less than 1, e.g. 0.5)

