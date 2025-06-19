# Options

There are quite a few options to choose from in DuctAPE.
As DuctAPE was developed, various options were added as different approaches were considered.
The current defaults were selected as the best set of options for the studies in the publications metioned in DuctAPE's README.
However, most of the options implemented during development have been maintained in case of future need.
The `Options` object contains various general options as well as several more detailed option objects.
To make the process of setting options easier, the `set_options` method can be used to make only the desired changes from the defaults without having to define everything in the `Options` object.
Note that the `Options` object is implemented using the `@kwdef` macro, so the `set_options` function doesn't really do anything in the case of a single operating point, but the other dispatch of `set_options` is especially helpful in initializing several of the sub-option objects to defaults of the correct size.
Note that there are several fields in the options that are used for bookkeeping, especially in the case of multiple operating points.

```@docs
DuctAPE.Options
DuctAPE.set_options
```

The major sub-categories of options include general options, pre-processing options, solver options for both determining the wake sheet positions as well as the overall aerodyanmics solve, post-processing options, and bookkeeping options.

---

## Bookkeeping Options

These are options that can be changed by the user for development/debugging purposes, but at this point, it would be wise in general usage to not change them. In future revisions, these will likely no longer be accessible to the user.

---

## General Options

The verbose and silence warnings options are simply about what get's printed as the analysis runs.
Warnings are printed when some sort of automated adjustment is made to the inputs in order to ensure they conform to the format required.
The verbose option at this level is for verbose statements that are not within any solvers.  Solver verbosity is constrolled in the individual solver options.

Occasionally, something in the preprocessing will fail, likely the LU decomposition of the linear system defining the bodies' panel system.
If such a failure occurs, DuctAPE cannot continue to the main solve and will exit.
The `hard_fail` option dictates what the exit behavior is.
If true, DuctAPE will just return `nothing` immediately, which is quicker for turn-around on single runs.
If false, DuctAPE will attempt to return an output object of the correct size and type, which is convenient for some optimization frameworks for which you'll want some output to be available even if passing a failure flag for the specific analysis.


---

## Preprocess Options

### Geometry Interpolation and Generation Options

The `autoshiftduct` option may be convenient depending on how the duct coordinates are being input.
It allows the user to input the duct coordinates at an arbitrary radial location, for example if a standard airfoil is used with leading edge at (0,0).
In the preprocessing, the duct geometry will be shifted to the radial location at which the rotor tip is coincident with the duct surface at the axial location at which the first rotor is situated.
If you are already inputting the duct geometry at the correct position, this option may be turned off, but it usually doesn't hurt to be left on.

### Paneling Options

These, in general, do not need to be touched by users; thus we do not include them in the `PanelingConstants`, but we do make them available.

### Integration Options

DuctAPE uses numerical integration to determine the influence of axisymmetric vortex and source panels on the restof the system.
There are several options for numerical integration, and the methods can be mixed and matched with their own specific options for the nominal (panel on other panels) and singular (panel on itself) cases.
The Gauss-Legendre method is useful in optimization cases to avoid noise in the integration error.
The Gauss-Kronrod method is implemented via QuadGK.jl and is more accurate, especially in for the singular cases, but is less useful for optimization purposes due to noise in the error from the adaptive nature.
Similarly, the Romberg method similar to that implemented in DFDC is adaptive and can be fast, but is also not usually the best choice.
Thus we default to Gauss-Legendre methods.

```@docs
DuctAPE.IntegrationOptions
DuctAPE.GaussLegendre
DuctAPE.GaussKronrod
DuctAPE.Romberg
```

#### Example

```@example quadrature
using DuctAPE

# set nominal options using a GaussLegendre object (which is an InterationMethod type)
# note that a convenience method is used here that takes in the number of points and
#calculates the appropriate sample locations and weights.
nominal_integration_method = DuctAPE.GaussLegendre(10)

# set singular options using a GaussKronrod object (which is an InterationMethod type)
# note that like most option structs, these are defined using @kwdef allowing the fields
#to be treated as keyword arguments.
# also note that we haven't changed the evaluation limit (default 10^7)
singular_integration_method = DuctAPE.GaussKronrod(; order=7, atol=2e-16)

# put the quadrature options together
integration_options = DuctAPE.IntegrationOptions(;
    nominal=nominal_integration_method, singular=singular_integration_method
)

# example of calling the set_options function
options = DuctAPE.set_options(; integration_options=integration_options)
```

---


## Solver Options

### Elliptic Grid Solver Options

The wake geometry is obtained by solving for approximate streamlines on an elliptic grid with the bodies as boundaries.
There are two methods available, a successive line over relaxiation (SLOR) method that can be used in isolation or as a preconditioner to a Newton solve method from the NLsolve.jl package.
The default option is to run a few iterations of SLOR to precondition and smooth out the initial grid, and then finish up with the Newton solve, usually within 3-5 iterations.
Note that the SLOR method is not implemented using ImplicitAD.jl since there isn't a clean residual definition separate from the solve method.
The Newton solve is implemented using ImplicitAD.jl which helps speed up automatic differentiation in an optimization setting, but it still benefits from a few iterations of SLOR beforehand.

```@docs
DuctAPE.SLORGridSolverOptions
DuctAPE.GridSolverOptions
```

#### Example
```@example gridsolve
using DuctAPE

# define wake grid solver settings
wake_solve_options = DuctAPE.GridSolverOptions(; atol=1e-10)

# set all options
options = DuctAPE.set_options(; grid_solver_options=wake_solve_options)
```

### Aerodynamics Solver Options

Quite a few solve methods were explored in the development of DuctAPE which can be separated into three broad categories: Fixed-point iteration solvers, quasi-Newton solvers, and Newton solvers.
In general, the fixed-point solvers have been faster and more robust than other methods, but all the methods have been kept in case they are desired for future development.


There are two methods that are implemented directly in DuctAPE: the CSOR and ModCSOR methods which are the controlled successive over relaxation fixed-point approach taken in DFDC and a modified version compatible with ImplicitAD.jl that is the current default and is currently the best (fastest/most robust) for optimization.

```@docs
DuctAPE.CSORSolverOptions
DuctAPE.ModCSORSolverOptions
```

The other methods are implemented via external dependencies and some do better than others.
In our experience, the other strictly fixed-point methods work relatively well, but are middle of the road.

```@docs
DuctAPE.SpeedMappingOptions
DuctAPE.FixedPointOptions
```

The quasi-Newton methods are hit or miss, with Minpack doing well enough to make it into some of the compound solver options discussed below, but we have had very little success wth SIAMFANLEquations up to this point.

```@docs
DuctAPE.MinpackOptions
DuctAPE.SIAMFANLEOptions
```

NonlinearSolve is generally faster than NLsolve if the problem is large enough, but we find it to be significantly less robust.  NLsolve's Anderson method is perhaps the best external method in terms of speed and robustness, but is just barely edged out by the CSOR methods.

```@docs
DuctAPE.NLsolveOptions
DuctAPE.NonlinearSolveOptions
```

Finally, there are several compound solve methods implemented, the first chaining solvers together.  The solver chain can be defined as the user wishes, but the defaults start with the fixed-point Anderson method, move to the Minpack quasi-Newton method if the fixed-point method doesn't converge in the given number of iterations, and then finishes with a full Newton method if the quasi-newton method doesn't converge.
The other compound solver combines solvers in a composite manner, typically starting with a few iterations of a full Newton method to get the solver going in the right direction and then finishing with a fixed-point method.
The `ChainSolverOptions` was at one point the default method, but once the modified CSOR method was developed, the compound solvers weren't used much.
Note that due to the way the solvers are implemented and dispatched, it is currently not possible to mix and match the CSOR methods with any of the external package methods.

```@docs
DuctAPE.ChainSolverOptions
DuctAPE.CompositeSolverOptions
```

#### Example
```@example
using DuctAPE
using LineSearches

# Define settings for NLsolve's newton method
aero_solver_options = DuctAPE.NLsolveOptions(;
    algorithm=:newton,
    atol=1e-10,
    iteration_limit=30,
    linesearch_method=LineSearches.BackTracking, #don't include parentheses on method handle
    linesearch_kwargs=(; order=3, maxstep=1e6),
)

# set all the options
DuctAPE.set_options(; solver_options=aero_solver_options)
```

!!! note "Iteration Counters"
    The `iterations` field (not to be confused with the `iterations_limit` field) in the solver options should generally not be changed.  They automatically save (in-place) the number of iterations the solver performs and can be accessed after the analysis is run.



---

## Postprocess Options

Most of the postprocess options have to do with writing the outputs to files, with the default behavior being to not write anything.  These options can be useful for debugging purposes or saving outputs, though it is usually more efficient to manually save a select few outputs rather than all the outputs.
The one option that is slightly more involved is the boundary layer options.
Currently, only Head's method is fully implemented, but there has also been some development started on Green's method.
The boundary layer options include choices regarding type of solver, with a simple 2nd-order Runge-Kutta method being the default (and appears to be best for optimization). There is also a 4th-order Runge-Kutta method implemented as well as the RadauIIA5 method from DifferentialEquations.jl which may be more accurate for single runs.
Note that the default setting is to not run the boundary layer method.
The `model_drag` option in the boundary layer options needs to be set to true if it is desired to include the drag model.

```@docs
DuctAPE.HeadsBoundaryLayerOptions
```


#### Example
```@example boundarylayer
using DuctAPE

# Define Boundary Layer Settings
boundary_layer_options = DuctAPE.HeadsBoundaryLayerOptions(;
    model_drag=true,
    separation_penalty_upper=0.1,
    separation_penalty_lower=0.1,
    separation_allowance_upper=3,
    separation_allowance_lower=25,
)

# set all the options
DuctAPE.set_options(; boundary_layer_options=boundary_layer_options)
```


---

# Advanced Options for Multi-point analyses

For using advanced options in multi-point analyses, there are various changes that need to be made to avoid run-time errors.
Here is an example for setting options with the CSOR solver.


```@example multipoint
using DuctAPE

# number of operating points to analyze
nop = 3

options = DuctAPE.set_options(;
    solver_options=DuctAPE.ModCSORSolverOptions(;
        converged=fill(false, (1, nop)), # need a convergence flag for each operating point
        iterations=zeros(Int, (1, nop)), # need a iteration count for each operating point
    ),
    write_outputs=fill(false, nop), # we need to know which of the operating point outputs to write
    outfile=fill("", nop), # we need to include names, even if they won't be used.
    output_tuple_name=fill("outs", nop), # we need to include names, even if they won't be used.
)
```

If using a compound algorithm with a multi-point solve, then each of the solvers needs to have the multiple `converged` and `iterations` fields for each operating point, and the overall solve type needs to have a `converged` and `iterations` field for each solver and each operating point.

```@example multipoint
options = DuctAPE.set_options(;
    solver_options=DuctAPE.ChainSolverOptions(;
        solvers=[ # vector of solvers to use in poly-algorithm
            DuctAPE.NLsolveOptions(;
                algorithm=:anderson,
                atol=1e-12,
                iteration_limit=200,
                converged=fill(false, (1, nop)), # flags for each operating point
                iterations=zeros(Int, (1, nop)), # counters for each operating point
            ),
            DuctAPE.MinpackOptions(;
                atol=1e-12,
                iteration_limit=100,
                converged=fill(false, (1, nop)),
                iterations=zeros(Int, (1, nop)),
            ),
        ],
        converged=fill(false, (2, nop)), # flags for each solver and each operating point
        iterations=zeros(Int, (2, nop)), # counts for each solver and each operating point
    ),
)
```
