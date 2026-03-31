
## Taylor-Green Vortex

## FIG1

- High Mach number leads to saturation due to compressibility error
- Before saturation, convergence at second order
- With bigger grids, the CFL reaches 1 soon, leading to instability
- Some methods diplay a mininum, others don't. Why?
- (The Bardow-FEM implementation could still be hiding a bug.)

## FIG2

Comparison with results of Zhu et al. (2017).
Serves as implementation check for DUGKS and Bardow's scheme.

- We used equilibrium initialization only; Zhu et al. (2017) used the non-equilibrium initilization by Skordos.
- The FEM method follows Bardow, Karlin, Gusev; we use row-lumping just like they did.
- Results closely match the values reported by Zhu et al. (2017), n.b. the overlapping markers.
- Convergence with second order in space.
- Bardow-FEM discretization shows early sign of saturation.

- In the future it would be interesting to compare further schemes schemes:
  * Bardow FEM with full mass matrix
  * Sofonea & Sekerka (2003), JCP 184, 422-434
  * Shrestha et al. (2016), PRE 93, 023306
  * Chen et al. (2020), CaMwA 80, 3066-3081, https://doi.org/10.1016/j.camwa.2020.10.022


TODO:

* run the FEM case at a bigger resolution
* add the mass matrix, to see difference in the spectral behavior
* check single-precision convergence instead
* instead of velocity norm, look at the fitted viscosity instead; also Shrestha did something like this. Perhaps this would make it easier to interpret the roots, where the errors drop to near zero.

Comment:

It is important to understand the error as a function of CFL number and dt/tau ratio.
For instance in the adaptive method of Fakhari, grid blocks run at various CFL numbers,
depending on the level of the grid, there will be CFL = 1.0 for the leaf blocks that run standard LBM.
Coarser blocks will have CFL values of 1/2, 1/4, ... depending on the level.

The hyper-convergent behavior also matters for making comparison between methods.
If unaware of these values, it may lead to false conclusions.


## FIG3

Comparison with results of Wu et al. (2018).
Also serves as an implementation check for DUGKS

u_0 = 1, nu = 0.01, dt = 1E-5, dt = 4 tau,

meaning that cs^2 = 4000, cs = 200.
Error is measured at t = 1s

This gives the settings

Re = 100
Ma = 1 / (20 sqrt(10)) =~ 0.0158

The total time simulates is t / dt = 100000 steps.
Alternatively, the time can be determined from the vortex decay time.

I must say here, that DUGKS is harder to parallelize. Due to the
intermediate collision at the cell edges, it favors a AoS memory layout.
With higher-order velocity lattices, there may be other challenges too.

## FIG4

Here we scrutinize the case picked in S&M (2024).

The settings are:
- Re = 1
- Mach = sqrt(3) * 0.01 (on the coarsest mesh)
- CFL (unclear)

In the original paper, CFL is varied with grid spacing, which
is a bit weird. I assume they wanted to keep the CFL fixed, as in
regular LBM where it is implicitly fixed to 1.

The high viscosity, means that the fluid dissipates very quickly, so
the number of time-steps needed is low.

On the other hand, the choice of parameters is such that dt / tau < 1,
meaning that the implicit collision technique is not really necessary here.
The whole benefit of the implicit collision, is being able to push
into high dt/tau regimes, which give good bandwidth.

## FIG5

The goal here was to repeat the study of Krämer et al. (2020), using the
Natrium library.

In Fig. 4 (page 44) of their article they show convergence of the
time-averaged errors in a simulation of a two-dimensional TGV.

The settings are
- Re = 10^6
- Ma = 10^(-4)

The time-step size is arbitrarily fixed at 10^(-7) and the simulation proceeds
for 100000 steps, that is until tmax = 0.01.

The time error was measured every tenth time step.

It is not specified exactly, how they compute the 2-norm; I'm guessing
it is a integral in the finite element sense.

Also for the time-averaging, they may have actually used the physical time-step.
When reporting the error, the area is assumed to be [0,2*pi]^2.

## FIG6

Comparison with results of Krämer.

Settings:
* Ma = 0.008
* CFL = 0.1

* nu = pi / 5

* tmax = ln(10) t_c, i.e. exp(-2 nu t_max) = 0.1

Domain [0,2 pi]^2, containing four vortices

This example is done in Excel

## FIG7

Comparison with results of Chen (2021)

## FIG8

Comparison with results of Guo & Zhao (2003)

## FIG9

Shear wave case of Sofonea & Sekerka (2003)



## Higher-order accuracy

Here the purpose is showing 4th and 6th order convergence are possible, and the error can be driven
all the way down to machine limits by use of diffusiving scaling (expensive!).

With the Lax-Wendroff discretization, there is no verdict on what is the "best" stencil.



