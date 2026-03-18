


import numpy as np
from kernel_tuner import tune_kernel



size = 512 * 512 * 9

fsrc = np.ones(size,dtype=np.float64)
fdst = np.zeros_like(fsrc)

dt = np.float64(0.01)

args = [fsrc,fdst,dt]

tune_params = dict()
tune_params["KERN"] = [1,2]
tune_params["NCOLLAPSE"] = [2,3]
tune_params["TILE_1"] = [1,2,4,8,16]
tune_params["TILE_2"] = [16,32,64,128,256]

print("compiling")
result, _ = tune_kernel(
    kernel_name="time",
    kernel_source="periodic_lw_mwe.F90",
    problem_size=(512,512),
    arguments=args,
    tune_params=tune_params,
    lang="C",
    compiler="gfortran-15",
    compiler_options=["-cpp","-fopenmp", "-O3", "-mcpu=native", "-ffast-math"])

