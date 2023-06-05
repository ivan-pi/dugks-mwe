# dugks-mwe

Compile time defintions available:
* `-DSPLIT` - use split collision kernel (no offloading available)
* `-DDUGKS` - use DUGKS mode (default is an FVM scheme)
* `-DBGK_OFFLOAD` - use target offload in the streaming kernel

```
$ make FC=gfortran FCFLAGS="-fopenmp -O3 -DBGK_OFFLOAD"
$ ./main_taylor_green <CFL>
```

`CFL` is the Courant number (a number less than 1, e.g. 0.5)

