# dugks-mwe

Build instructions:
```
cmake -B build --preset gnu-release [-DDUGKS=ON|OFF]
cmake --build build -j 4
```

Preprocessor flags:
* `-DSPLIT` - use split collision kernel (no offloading available)
* `-DDUGKS` - use DUGKS mode (default is an FVM scheme)
* `-DBGK_OFFLOAD` - use target offload in the streaming kernel
* `-DWITH_SP` - use single precision (32-bit)

Environment variables:
* `N` - grid size
* OpenMP variables (`OPENMP_NUM_THREADS`)

To-Do:
* Better setting parsing
* Implement target offload in collision step
* OpenMP data regions (`lattice.F90`)

To run the progam use,
```
N=120 ./main_taylor_green <CFL>
```
where `CFL` is the Courant number (< 1)

---

Older Make-based build:
```
$ make FC=gfortran FCFLAGS="-fopenmp -O3 -march=native"
```
The flags `LDFLAGS` and `LDLIBS` can be modified if needed.
