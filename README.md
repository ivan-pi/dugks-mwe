# dugks-mwe

Preprocessor flags:
* `-DSPLIT` - use split collision kernel (no offloading available)
* `-DDUGKS` - use DUGKS mode (default is an FVM scheme)
* `-DBGK_OFFLOAD` - use target offload in the streaming kernel

Environment variables:
* `N` - grid size
* OpenMP variables (`OPENMP_NUM_THREADS`)

```
$ make FC=gfortran FCFLAGS="-fopenmp -O3 -DBGK_OFFLOAD"
$ N=120 ./main_taylor_green <CFL>
```

`CFL` is the Courant number (a number less than 1, e.g. 0.5)

