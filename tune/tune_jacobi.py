# tune_jacobi.py
import numpy as np
from kernel_tuner import tune_kernel

N = 1024
size = N * N

fsrc = np.random.random_sample(size).astype(np.float64)
fdst = np.zeros_like(fsrc)

args = [fsrc,fdst]

tune_params = {
    "NCOLLAPSE" : [1, 2],
    "TILE1" : [1, 2, 4, 8, 16, 32],
    "TILE2" : [64, 128, 256, 512],
}
#    "TILE3" : [128, 256, 512]

print("compiling")
result, _ = tune_kernel(
    kernel_name="time_jacobi",
    kernel_source="jacobi_kernels.F90",
    problem_size=(N, N),
    arguments=args,
    tune_params=tune_params,
    lang="C",
    compiler="gfortran-15",
    compiler_options=["-cpp","-fopenmp", "-O3", "-mcpu=native", "-ffast-math"])
