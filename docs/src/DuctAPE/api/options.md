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

## Bookkeeping Options

These are options that can be changed by the user for development/debugging purposes, but at this point, it would be wise in general usage to not change them. In future revisions, these will likely no longer be accessible to the user.

## General Options

The verbose and silence warnings options are simply about what get's printed as the analysis runs.
Warnings are printed when some sort of automated adjustment is made to the inputs in order to ensure they conform to the format required.
The verbose option at this level is for verbose statements that are not within any solvers.  Solver verbosity is constrolled in the individual solver options.

Occasionally, something in the preprocessing will fail, likely the LU decomposition of the linear system defining the bodies' panel system.
If such a failure occurs, DuctAPE cannot continue to the main solve and will exit.
The `hard_fail` option dictates what the exit behavior is.
If true, DuctAPE will just return `nothing` immediately, which is quicker for turn-around on single runs.
If false, DuctAPE will attempt to return an output object of the correct size and type, which is convenient for some optimization frameworks for which you'll want some output to be available even if passing a failure flag for the specific analysis.

## Preprocess Options



### Geometry Interpolation and Generation Options

The `autoshiftduct` option may be convenient depending on how the duct coordinates are being input.
It allows the user to input the duct coordinates at an arbitrary radial location, for example if a standard airfoil is used with leading edge at (0,0).
In the preprocessing, the duct geometry will be shifted to the radial location at which the rotor tip is coincident with the duct surface at the axial location at which the first rotor is situated.
If you are already inputting the duct geometry at the correct position, this option may be turned off, but it usually doesn't hurt to be left on.

### Paneling Options

These, in general, do not need to be touched by users; thus we do not include them in the `PanelingConstants`, but we do make them available.

### Integration Options



```@docs
DuctAPE.IntegrationOptions
```

```@docs
DuctAPE.GaussLegendre
DuctAPE.GaussKronrod
DuctAPE.Romberg
```

## Solver Options



### Elliptic Grid Solve
```@docs
DuctAPE.SLORGridSolverOptions
DuctAPE.GridSolverOptions
```

### Aerodynamics Solve
```@docs
DuctAPE.ChainSolverOptions
DuctAPE.CompositeSolverOptions
DuctAPE.NLsolveOptions
DuctAPE.NonlinearSolveOptions
DuctAPE.MinpackOptions
DuctAPE.SIAMFANLEOptions
DuctAPE.SpeedMappingOptions
DuctAPE.FixedPointOptions
DuctAPE.CSORSolverOptions
DuctAPE.ModCSORSolverOptions
```


## Postprocess Options
```@docs
DuctAPE.BoundaryLayerOptions
DuctAPE.HeadsBoundaryLayerOptions
```

# Advanced Option Selection

DuctAPE has been written in an attempt to make as many of the available options exposed to the user as possible.  This means that there are quite a few options to select from if not using the option convenience functions.
To help the user, the majority of overarching option types are defined using the `@kwdef` macro and have default values that should be reasonable in most cases.
We will introduce some of the available options here that may be of common interest.

## General Option Selection

In general, options are all accessed through the `options` argument of the analysis function being called.
Said options are passed via an `Options` struct.

```@docs; canonical=false
DuctAPE.Options
```

Options are selected through the `set_options` function

```@docs; canonical=false
DuctAPE.set_options
```

There are three main sub-option objects for quadrature, wake geometry solver, and aerodyanmic solver; these are explained in more detail below.
In addition, there are various options for pre- and post-processing as well as miscellaneous options for things such as supressing warnings and printing verbose statements throughout the analysis, which can be seen in the docstring above.


## Quadrature

There are several implementations for different quadrature approaches depending on user desires; they include:
- [Gauss-Legendre quadature](@ref "DuctAPE.GaussLegendre") (default),
- [Gauss-Kronrod Quadrature](@ref "DuctAPE.GaussKronrod"), and
- [Romberg Quadrature](@ref "DuctAPE.Romberg") methods.

The default method is Gauss-Legendre quadrature using 8 sample points for both the nominal and singular integrals.
To modify the quadrature methods and settings, an `IntegrationOptions` struct needs to be passed to the `set_options` method.

```@docs; canonical=false
DuctAPE.IntegrationOptions
```

The `IntegraionOptions` type takes in two objects of type `IntegrationMethod`, one for the nominal integrals, and one for the singular integrals.
These methods can be mixed and matched between quadrature methods as well as settings.

For example, if one wanted to use a 10-point Gauss-Legendre method for the nominal integrals, and a order 7 Gauss-Kronrod method with an absolute tolerance of 2e-16 the following would need to be included in the `set_options` call:

```julia
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

## Elliptic Grid Solvers

As part of the pre-process, an elliptic grid defining the wake geometry is solved with a system of Poisson equations.
For this solve there currently two options:

- [SLOR](@ref "DuctAPE.SLORGridSolverOptions"): DFDC grid solver
- [Default](@ref "DuctAPE.GridSolverOptions"): Default method compatible with ImplicitAD

The SLOR (successive line over relaxation) is the method employed by DFDC.

Selection of solver and solver settings follows the same pattern as with the quadrature settings, in that the user must pass the appropriate `GridSolverOptionsType` into the `set_options` call.

For the SLOR method alone, the type is
```@docs; canonical=false
DuctAPE.SLORGridSolverOptions
```

And for the default method compatible with ImplicitAD, the type is
```@docs; canonical=false
DuctAPE.GridSolverOptions
```

As an example, this is the input that would be required to use the default method with an absolute convergence tolerance of 1e-10, and also including the quadrature settings from above:

```julia
# define wake grid solver settings
wake_solve_options = DuctAPE.GridSolverOptions(; atol=1e-10)

# set all options
options = DuctAPE.set_options(;
    integration_options=integration_options, grid_solver_options=wake_solve_options
)
```

!!! note "Convergence Flags"
    The convergence flags default to false, and in general should be left alone as they are modified in-place in the various solves by the analysis.

## Aerodynamics Solvers

There are two general types of solvers available in DuctAPE, the first is very similar to the solver in DFDC and converges a residual very similar to DFDC's.
The other type is for external solvers that converge an alternate residual that is default in DuctAPE.
The various solver options include:
- [CSOR](@ref "DuctAPE.CSORSolverOptions"): the DFDC solver
- [ModCSOR](@ref "DuctAPE.ModCSORSolverOptions"): modified DFDC solver for ImplicitAD compatibility
- [FixedPoint.jl](@ref "DuctAPE.FixedPointOptions")
- [SpeedMapping.jl](@ref "DuctAPE.SpeedMappingOptions")
- [MINPACK.jl](@ref "DuctAPE.MinpackOptions")
- [SIAMFANLEquations.jl](@ref "DuctAPE.SIAMFANLEOptions")
- [NLsolve.jl](@ref "DuctAPE.NLsolveOptions")
- [SimpleNonlinearSolve.jl](@ref "DuctAPE.NonlinearSolveOptions")

Note that the CSOR, ModCSOR, FixedPoint.jl, and SpeedMapping.jl are all different fixed-point iteration solvers, MINPACK.jl and SIAMFANLEquations.jl are primarily quasi-newton solvers, and NLsolve.jl and SimpleNonlinearSolve.jl have various solver options.

DuctAPE also has some poly-algorithm solvers that employ more than one solver.
The [Chain Solver](@ref "DuctAPE.ChainSolverOptions") option is the default which starts with a fixed-point iteration, and if it doesn't converge, moves on to a quasi-, then full Newton solver until either convergence is reached, or no convergence is found.
The other poly-algorithm that is available, but is less robust is the [Composite Solver](@ref "DuctAPE.CompositeSolverOptions") which partially converges with one solver, and finishes with another.

Each of the solve methods have a variety of different settings associated with them, detailed in their respective docstrings.
The following example should contain all the principles required to be able to adapt to the most complex use cases.

```julia
# Define settings for NLsolve's newton method
aero_solver_options = DuctAPE.NLsolveOptions(;
    algorithm=:newton,
    atol=1e-10,
    iteration_limit=30,
    linesearch_method=LineSearches.BackTracking, #don't include parentheses on method handle
    linesearch_kwargs=(; order=3, maxstep=1e6),
    additional_kwargs=(; autoscale=false),
)

# set all the options
DuctAPE.set_options(;
    integration_options=integration_options,
    grid_solver_options=wake_solve_options,
    solver_options=aero_solver_options,
)
```

!!! note "Iteration Counters"
    The `iterations` field (not to be confused with the `iterations_limit` field) in the solver options should generally not be changed.  They automatically save (in-place) the number of iterations the solver performs and can be accessed after the analysis is run.

## Boundary Layer Solvers

If desired, a one-way turbulent boundary layer can be modeled, from which an approximate viscous drag can be determined.
Currently, only the Head's method for turbulent boundary layer computation is working.
In the case of separation, a separation drag penalty is applied based on values selected in the options.

```@docs; canonical=false
DuctAPE.HeadsBoundaryLayerOptions
```

Here is an example of a possible boundary layer option setting:

```julia
# Define Boundary Layer Settings
boundary_layer_options = DuctAPE.HeadsBoundaryLayerOptions(;
    model_drag=true,
    separation_penalty_upper=0.1,
    separation_penalty_lower=0.1,
    separation_allowance_upper=3,
    separation_allowance_lower=25,
)

# set all the options
DuctAPE.set_options(;
    integration_options=integration_options,
    grid_solver_options=wake_solve_options,
    solver_options=aero_solver_options,
    boundary_layer_options=boundary_layer_options,
)
```

## Advanced Options for Multi-point analyses

For using advanced options in multi-point analyses, there are various changes that need to be made to avoid run-time errors.
Here is an example for setting options with the CSOR solver.

```julia
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

If using a poly-algorithm with a multi-point solve, then each of the solvers needs to have the multiple `converged` and `iterations` fields for each operating point, and the overall solve type needs to have a `converged` and `iterations` field for each solver and each operating point.

```julia
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
