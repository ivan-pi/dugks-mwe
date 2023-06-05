# dugks-mwe

```
$ make FC=gfortran FCFLAGS="-fopenmp -O3"
$ ./main_taylor_green <CFL>
```

`CFL` is the Courant number (a number less than 1, e.g. 0.5)