# MHD2D Solver

PDE solver for 2D ideal MHD problems, using the ENO (Essentially Non-Oscillatory) scheme.

This is a research project, it has some rough edges and is not very user friendly, but it should be functional.

The source is an unholy amalgamation of C-like old code and partially rewritten object-oriented code. My apologies.

It has also not been optimized to avoid unnecessary bugs (unit tests? what is that?), and parallelization consists of "run several copies with different input files". That being said, as there is a "Publications" section at the end of this document, it did its job pretty well.

(Note: this document can be converted to HTML with, for example, `pandoc -f markdown_github -t html README.md -o README.html`.)

## Table of contents

* [Requirements](#requirements)
* [Compilation](#compilation)
  * [IDE](#ide)
  * [Command line](#command-line)
* [Usage](#usage)
  * [Section `[problem]`](#section-problem)
    * [Test problem: shock tube](#test-problem-shock-tube)
    * [Test problem: 2D explosion](#test-problem-2d-explosion)
    * [Test problem: Orszag-Tang vortex](#test-problem-orszag-tang-vortex)
    * [Test problem: Linear waves](#test-problem-linear-waves)
    * [Plasma sheet](#plasma-sheet)
  * [Section `[divb corrector]`](#section-divb-corrector)
  * [Section `[time]`](#section-time)
  * [Section `[output grid]`](#section-output-grid)
  * [Section `[output non grid]`](#section-output-non-grid)
  * [Section `[logging]`](#section-logging)
* [Subprojects](#subprojects)
  * [Euler 1D solver](#euler-1d-solver)
  * [RK3 order](#rk3-order)
* [Publications](#publications)

## Requirements

* Boost ( https://www.boost.org/ ), tested with version 1.55.0. It doesn't have to be compiled, only header files are needed.

* GCC, any version supporting at least C++11. Other compilers should work, but there are no prepared Makefiles or project files for them.

* Optional: The Code::Blocks IDE ( https://www.codeblocks.org/ ).

## Compilation

### IDE

The solver was developed in the
[Code::Blocks](https://www.codeblocks.org/) IDE.

To compile in C::B, first open the project file, `mhd_2d.cbp`.
Then, go to `Project > Build options...`, and in the left panel select the root node of the project (`MHD2D Solver`).
On the right, in the `Compiler settings` tab, open the `Other compiler options` sub-tab and update the `-isystem` flag to point to the directory with the Boost header files.

You should now be able to compile the project with the C::B's `Build` command.

### Command line

An Makefile is included.
Add the Boost headers to the compiler's search path and run `make` in the source directory.
Alternatively, update the Makefile to include the `-isystem path/to/boost_1_55_0` flag in the compilation step.

(Aside: Why `-isystem` and not `-I`? Because including it as a system library suppresses the compilation warnings from Boost headers.)

## Usage

```
./mhd2d input-file.ini
```

The simulation is configured with an input file, given in INI format.
The configuration is split into several sections:

* `[problem]`: Problem to solve and its parameters.
* `[divb corrector]`: Parameters for the subroutine for cleaning the magnetic field divergence.
* `[time]`: Time stepping parameters.
* `[output grid]`: Output parameters for full simulation grid.
* `[output non grid]`: Output parameters for aggregate data.
* `[logging]`: Turn on or off several logging categories, mainly used for debugging.

If a configuration option is missing, the program will print a warning and use a hardcoded default value.

The `example-inputs` directory contains a handful of typical input files.

The `papers-inputs` directory contains (most of) the input files and scripts used to analyse results presented in [publications](#publications).
In subdirectory `2d` are the input files for the main project, while in subdirectories `1d` and `rk3-order` are the input files for the two [subprojects](#subprojects).

### Section `[problem]`

There are several predefined test problems, and the core "plasma sheet" configuration.
The problem type only defines the initial state of the simulation grid, the same simulation code is used for all problems.

#### Test problem: shock tube

A 1D problem where two fluids in a tube are separated by a membrane, which is removed at time `t = 0`.

```
[problem]
type = shock tube

Nx =   32
Ny =    8
Lx =  2.0
Ly =  0.5
start_x = -1.0
start_y = -0.25

rho_l = 1.0
u_l   = 0.0
v_l   = 0.0
w_l   = 0.0
Bx_l  = 0.75
By_l  = 1.0
Bz_l  = 0.0
p_l   = 1.0

rho_r = 0.125
u_r   = 0.0
v_r   = 0.0
w_r   = 0.0
Bx_r  = 0.75
By_r  = -1.0
Bz_r  = 0.0
p_r   = 0.1

gamma = 2.0
time method = rk3
space method = eno-lf
boundary left   = neumann
boundary right  = neumann
boundary top    = neumann
boundary bottom = neumann

halt on negative pressure = true
```

`Lx`, `Ly` (doubles) determine the physical size of the simulation box, while `Nx`, `Ny` (integers) are the total amount of grid points in each direction
(`dx = Lx / Nx` and `dy = Ly / Ny` can be different).
`start_x`, `start_y` determine the coordinates of the lower left corner of the simulation grid.

The `_l` and `_r` variables describe, respectively, the plasma state on the left and right side of the shock tube, and the `interface angle` is the angle between the shock tube and the _x_ axis, given in degrees.

`gamma` is the ratio of specific heats.

`time method` is the method used for time-stepping: `euler` (Euler method) or `rk3` (third order total variation diminishing Runge-Kutta method).

`space method` is the method used for spatial integration: `central fd` (second order central finite difference), `eno-roe` (ENO method, Roe version), or `eno-lf` (ENO method with Lax-Friedrichs flux splitting).

The four `boundary` parameters define the conditions at their respective boundaries: `periodic` (periodic BC), `dirichlet` (Dirichlet BC), `neumann` (Neumann BC), or `open` (no BC is imposed).

If the `halt on negative pressure` parameter is `true`, the simulation will abort if negative pressure is encountered during calculation. If it is `false`, negative pressure is forced to zero and calculation resumes.

#### Test problem: 2D explosion

A 2D version of the shock tube problem, where two fluids are separated by a circular membrane, which is removed at time `t = 0`.

```
[problem]
type = explosion

Nx =  64
Ny =  64
Lx = 1.0
Ly = 1.0
start_x = -0.5
start_y = -0.5

transition type   = linear
transition points = 1

r = 0.2

rho_in = 1.0
u_in   = 0.0
v_in   = 0.0
w_in   = 0.0
Bx_in  = 0.2
By_in  = 0.0
Bz_in  = 0.0
p_in   = 1.0

rho_out = 0.125
u_out   = 0.0
v_out   = 0.0
w_out   = 0.0
Bx_out  = 0.2
By_out  = 0.0
Bz_out  = 0.0
p_out   = 0.1

gamma = 1.66666666667
time method = rk3
space method = eno-lf
boundary left   = neumann
boundary right  = neumann
boundary top    = neumann
boundary bottom = neumann

halt on negative pressure = true
```

`transition type` describes the interface between the two fluids: `jump` (a sharp transition) or `linear` (linearly interpolated transition).
If `transition type = linear`, an additional parameter `transition points` gives the width of the transition.

`r` is the radius of the membrane (i.e., the inner fluid).

The `_in` variables describe the state of the inner fluid, while the `_out` variables describe the outer fluid.

The rest of the parameters are the same as in the shock tube problem.

#### Test problem: Orszag-Tang vortex

The Orszag-Tang vortex is a standard test problem used for testing numerical simulations of plasma.

```
[problem]
type = orszag-tang vortex

Nx = 64
Ny = 64

x0 percent = 0.0
y0 percent = 0.0

gamma = 1.66666666667
time method = rk3
space method = eno-lf

halt on negative pressure = true
```

The OT vortex is defined on a square periodic grid of size 2π×2π.

The `x0 percent` and `y0 percent` determine the offset of the initial state from the canonical one; an offset of 0.5 would shift the box by π.

The rest of the parameters are the same as in the shock tube problem.

#### Test problem: Linear waves

For small-amplitude waves, the system is nearly linear, and waves travel at a constant velocity.
This problem can be used to measure the accuracy of the simulation;
with periodic boundary conditions, an excited wave moving at velocity `U` will return to the initial position at time `T = Lx/U`,
and the difference between the initial state and state at time `T` allows us to estimate the error.

```
[problem]
type = wave test

Nx = 128
Ny = 128
Lx = 1.0
Ly = 1.0
start_x = -0.5
start_y = -0.5

rho = 1.0
u   = 0.0
v   = 0.0
w   = 0.0
Bx  = 1.0
By  = 1.4142135623730950488016887242097
Bz  = 0.5
p   = 0.6

wave         = negative fast magnetosonic
wavenumber x = 1
wavenumber y = 0
epsilon      = 1.0e-7

gamma = 1.66666666667
time method = rk3
space method = eno-lf
boundary left   = periodic
boundary right  = periodic
boundary top    = periodic
boundary bottom = periodic

halt on negative pressure = true
```

The `rho` to `p` values give the uniform background state of the plasma.

`wave` determines which of the seven MHD waves will be excited: `negative fast magnetosonic`, `negative alfven`, `negative slow magnetosonic`, `entropy`, `positive slow magnetosonic`, `positive alfven`, `positive fast magnetosonic`.

`wavenumber x` and `wavenumber y` determine how many waves will be excited in each direction (i.e., the wavelength),
and `epsilon` determines the amplitude.

The rest of the parameters are the same as in the shock tube problem.

#### Plasma sheet

The core problem this simulation code was developed for.

Initially, there are four distinct plasma states, "up", "down", "left", and "right".
The physical configuration of the simulation area is show in the following diagram:
```
+-------------------------+
|            up           |
+------------+------------+
|    left    |   right    |
+------------+------------+
|           down          |
+-------------------------+
```

```
[problem]
type = plasma sheet

Nx =  512
Ny =   96
Lx = 32.0
Ly =  6.0
start_x = -16.0
start_y = -3.0
sheet start      = 0.0
sheet start foot = 0.0
sheet thickness  = 1.0
sheet profile    = constant

transition type   = linear
transition points = 1

rho_l = 1.0
u_l   = -0.5
v_l   = 0.0
w_l   = 0.0
Bx_l  = 0.0
By_l  = 0.0
Bz_l  = 0.0
p_l   = 1.0

rho_r = 1.0
u_r   = 0.0
v_r   = 0.0
w_r   = 0.0
Bx_r  = 0.0
By_r  = 0.0
Bz_r  = 0.0
p_r   = 1.0

rho_u = 0.155
u_u   = 0.0
v_u   = 0.0
w_u   = 0.0
Bx_u  = -1.3
By_u  = 0.0
Bz_u  = 0.0
p_u   = 0.155

rho_d = 0.155
u_d   = 0.0
v_d   = 0.0
w_d   = 0.0
Bx_d  = 1.3
By_d  = 0.0
Bz_d  = 0.0
p_d   = 0.155

gamma = 1.66666666667
time method = rk3
space method = eno-lf
boundary left   = dirichlet
boundary right  = dirichlet
boundary top    = neumann
boundary bottom = neumann

halt on negative pressure = true
```

`sheet start` gives the _x_ coordinate of the boundary between the "left" and "right" area.
`sheet start foot` gives the _x_ coordinate of the top and bottom edge of the aforementioned boundary; if `sheet start foot < sheet start`, the boundary itself is a half-ellipse; if `sheet start foot = sheet start` it is a straight line (bug: the program crashes for `sheet start foot > sheet start`).
`sheet thickness` is the height of the "left" and "right" areas.

`sheet profile` determines the vertical profile of the magnetic field in the "left" and "right" areas: `constant` (magnetic field is uniform over the entire area) or `tanh` (magnetic field scaled by `tanh(y)`, where `y` is normalized so that it equals `1` at the lower edge and `1` at the upper edge; pressure is also adjusted to keep total pressure uniform).

`transition type` describes the interface between adjacent areas: `jump` (a sharp transition) or `linear` (linearly interpolated transition).
If `transition type = linear`, an additional parameter `transition points` gives the width of the transition. (Note: a linear transition does not work well if `sheet start foot != sheet start`.)

The `_l`, `_r`, `_u`, and `_d` variables describe, respectively, the "left", "right", "up", and "down" plasma states.

The rest of the parameters are the same as in the shock tube problem.

### Section `[divb corrector]`

Due to the peculiarities of the simulation method used, accumulated numerical errors may generate non-physical divergence in the magnetic field.
This program uses the SOR (Successive Over-Relaxation) method to solve a Poisson equation whose solution is used to remove the divergence.

```
[divb corrector]
stepping rate = 1
method = sor

sor rsteps = 100000
sor rmax   = 1.0e-5
sor steps    = 10
sor divbmax  = 1.0e-6
sor overrelaxation param use static = true
sor overrelaxation param            = 1.6
```

`stepping rate` determines how often (i.e., every `stepping rate` steps) the correction routine runs.

`method` may eventually allow choosing other correction methods: currently, the only allowed value is `sor`.

SOR parameters are:

* `sor rsteps`: maximum allowed iterations; if the solution doesn't converge in this many steps, the program will print a warning and resume calculation.
* `sor rmax`: target relative error between iterations; if the relative error is smaller than this value, the SOR method has converged.
* `sor steps`: maximum number of times to iterate the cleaning procedure if the result is not satisfactory.
* `sor divbmax`: if the largest divergence in the simulation area is smaller than this value, the cleaning procedure has succeeded.
* `sor overrelaxation param use static`: if `true`, read the overrelaxation parameter from configuration (see below); if `false`, let the program determine it dynamically. (Note: the dynamic method is suspicious, as we couldn't find the paper used to implement it.)
* `sor overrelaxation param`: the over-relaxation parameter used in the SOR method.


### Section `[time]`

Time stepping parameters.

```
[time]
mode = variable
cfl number = 0.2
dt_max = 0.0078125
dt_min = 0.000001
t_max = 1.0
```
or
```
[time]
mode = constant
steps = 1000
dt    = 0.001
```

`mode` is how the time step sizes are determined: `constant` (fixed time steps) or `variable` (time steps are halved or doubled until the the CFL condition is satisfied).

If `mode = constant`:

* `steps`: how many steps to calculate.
* `dt`: time step size.

If `mode = variable`:
* `cfl number`: the CFL number (obtained from velocities in the simulation) has to be below this value.
* `dt_max`: the largest allowed time step.
* `dt_min`: the smallest allowed time step; if the CFL condition is not satisfied, the simulation will terminate with an error.
* `t_max`: target end time of the simulation (note: it may go over by a single time step).

### Section `[output grid]`

How and when to output the full grid of simulation data.

```
[output grid]
natural      = true
conservation = false
single file  = false
datafile     = %m-%ri-%rt.dat
rt format = %010.7f
ri format = %04d
rs format = %08d
binary            = true
binary use double = false
binary hint       = %m.hint
skip steps = 1
skip t     = 0.0078125
mode = time
```

`natural` determines whether to output the natural plasma variables (velocities `u, v, w`, pressure `p`).
`conservation` determines whether to output their conserved counterparts (moments `mx = rho u, my = rho v, mz = rho w`, energy `e`).
Those that belong in both categories (density `rho`, magnetic field `Bx, By, Bz`) are outputted unconditionally.

`single file` determines if the results will all be in one file, or if there will be a separate file for each result.
The name of the file(s) is set with `datafile`, where the following substitutions can be used:

* `%f`: full name of the input file.
* `%m`: name of the input file without extension.
* `%t`: current wall clock time.
* `%rs`: current step (multi-file only).
* `%ri`: file index (multi-file only).
* `%rt`: current simulation time (multi-file only).

`rt format`, `ri format`, and `rs format` are `printf`-style format templates for their respective fields.

`binary` determines whether the output file(s) will be in binary or text format.
In text format, the filed values are printed in exponential format to 6 significant digits (`%+.5e`).
In binary format, if `binary use double` is `true` then the data will be outputted as doubles (8 bytes); if it is `false`, then floats will (4 bytes) be used.
When binary format is used, the program will also provide field information and some example plots in a file named with `binary hint`; the same `%f, %m, %t` substitutions as above can be used. Field information for the text format is shown in the header row of result file(s).

The timing of file output is determined by `mode`, which can be either `time` (output when at least `skip t` time has passed since last output) or `step` (output every `skip steps` steps; yes, the naming is confusing).

### Section `[output non grid]`

Program can also output some basic and aggregate data: how many full results have been printed, current step and simulation time, step size, maximum and total divergence of the magnetic field, total energy, etc.

```
[output non grid]
datafile = %m.ngdat
skip steps = 1
skip t     = 0.001
mode = step
```

Name of the output file is set with `datafile`; the same `%f, %m, %t` substitutions as above can be used.

`mode` and `skip t` / `skip steps` have the same function as in the previous section.

### Section `[logging]`

Finally, some values can be displayed to help with debugging (separate from the full debugging mode set by `DEBUG` and friends).

```
[logging]
phi correction       = false
phi correction final = false
divb correction      = false
characteristics      = false
```

`phi correction` will print to screen the value of the relative error in the SOR method for each iteration.
`phi correction final` will report whether the SOR method converged and print the final value.
`divb correction` will report the maximum value of magnetic field divergence after the cleaning procedure.

`characteristics` will output the eigenvalues (wave velocities) and local values of characteristics to a file with the same name as grid result filename, with `char-` prefixed. (Note: it turns out that this data is pretty much useless, but it has been left in the code because it was too much trouble to remove.)

## Subprojects

### Euler 1D solver

Solver for Euler equations of gas dynamics. Located in directory `euler-1d-solver`.

Compiled and run the same way as the main MHD2D project.

The `example-inputs` subdirectory has several input files covering most of the possible input configurations.

### RK3 order

A small program to test the order of the three-step Runge-Kutta method by calculating a numerical approximation of `exp(x)` from 0 to 1.

Can be compiled as Code::Blocks project, or from the command line with:
```
gcc -std=c++14 -o rk3-order main.cpp
```

Takes one positive integer parameter, the number of steps to divide the interval into, and outputs the numerical error of Euler method and three-step Runge-Kutta method.

## Publications

[1] Rudolf Tretler, "Two-dimensional MHD simulation of the piston-induced rarefaction wave in the Earth's plasma sheet", Master's thesis, University of Electro-Communications, Tokyo, March 2015. ( http://id.nii.ac.jp/1438/00005010/ )

[2] Rudolf Tretler, Tomo Tatsuno, and Keisuke Hosokawa, "Loss of the rarefaction wave during plasma sheet thinning", Plasma and Fusion Research, vol. 15, 2401053, 2020. ( http://www.jspf.or.jp/PFR/PFR_articles/pfr2020/pfr2020_15-2401053.html )

[3] Rudolf Tretler, Tomo Tatsuno, and Keisuke Hosokawa, "Plasma sheet thinning due to loss of near-Earth magnetotail plasma", Journal of Plasma Physics, vol. 87, 905870108, 2021. ( https://doi.org/10.1017/S0022377820001580 )

[4] Rudolf Tretler, "Simulation of a simplified 2D MHD model of plasma sheet thinning", PhD thesis, University of Electro-Communications, Tokyo, March 2021. ( https://doi.org/10.18952/00009893 )
