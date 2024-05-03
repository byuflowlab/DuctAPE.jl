# Examples

```@contents
Pages = ["examples.md"]
Depth = 5
```

## Advanced Option Selection

DuctAPE has been written in an attempt to make as many of the available options exposed to the user as possible.  This means that there are quite a few options to select from if not using the option convenience functions.
To help the user, the majority of overarching option types are defined using the `@kwdef` macro and have default values that should be reasonable in most cases.
We will introduce some of the available options here that may be of common interest.

### Quadrature

There are several implementations for different quadrature approaches depending on user desires; they include:
- [Gauss-Legendre quadature](@ref "DuctAPE.GaussLegendre") (default),
- [Gauss-Kronrod Quadrature](@ref "DuctAPE.GaussKronrod"), and
- [Romberg Quadrature](@ref "DuctAPE.Romberg") methods.

The docstrings for each are linked to the names.

### Elliptic Grid Solvers

As part of the pre-process, an elliptic grid defining the wake geometry is solved with a system of Poisson equations.
For this solve there currently two options:

- [SLOR](@ref "DuctAPE.SLORGridSolverOptions"): DFDC grid solver
- [SLOR+Newton](@ref "DuctAPE.GridSolverOptions")

The SLOR (successive line over relaxation) is the method employed by DFDC, and can be used by itself, or as a preconditioner to a Newton solve (using NLsolve.jl).

### Aerodynamics Solvers

There are two general types of solvers available in DuctAPE, the first is very similar to the solver in DFDC and converges a residual very similar to DFDC's.
The other type is for external solvers that converge an alternate residual that is default in DuctAPE.
The various solver options include:
- [CSOR](@ref "DuctAPE.CSORSolverOptions"): the DFDC solver
- [FixedPoint.jl](@ref "DuctAPE.FixedPointOptions")
- [SpeedMapping.jl](@ref "DuctAPE.SpeedMappingOptions")
- [MINPACK.jl](@ref "DuctAPE.MinpackOptions")
- [SIAMFANLEquations.jl](@ref "DuctAPE.SIAMFANLEOptions")
- [NLsolve.jl](@ref "DuctAPE.NLsolveOptions")
- [SimpleNonlinearSolve.jl](@ref "DuctAPE.NonlinearSolveOptions")

The respective docstrings are linked to their names.
Note that the CSOR, FixedPoint.jl, and SpeedMapping.jl are all different fixed-point iteration solvers, MINPACK.jl and SIAMFANLEquations.jl are primarily quasi-newton solvers, and NLsolve.jl and SimpleNonlinearSolve.jl have various solver options.

DuctAPE also has some poly-algorithm solvers that employ more than one solver.
The [Chain Solver](@ref "DuctAPE.ChainSolverOptions") option is the default which starts with a fixed-point iteration, and if it doesn't converge, moves on to a quasi-, then full Newton solver until either convergence is reached, or no convergence is found.
The other poly-algorithm that is available, but is less robust is the [Composite Solver](@ref "DuctAPE.CompositeSolverOptions") which partially converges with one solver, and finishes with another.

### Other Options

The remaining option details can be found in the [Options](@ref "DuctAPE.Options") docstring.


-----


## Available Outputs

The output tuple contains many items.
The [`post_process`](@ref "DuctAPE.post_process") function docstring lists them.

### Returning the Pre-process Objects

Sometimes, it may be desireable to return the pre-process objects, including:

- `panels` which is a named tuple containing the body, rotor, and wake panel objects
- `ivb` which are the unit induced velocities on the body panels
- `solve_parameter_tuple` which contains all of the solver parameters
- `blade_elements` which contains all of the blade element geometry and airfoil information
- `linsys` which contains all the linear system objects for the panel method
- `idmaps` which contains all the index mapping used throughout the solve and post-process.

In this case, we can use the `return_inputs` keyword argument when calling the `analyze` function to return a named tuple containing those pre-process objects.

```julia
outs, ins, success_flag = dt.analyze(propulsor; return_inputs=true)
```


-----


## Multi-Point Analyses

In the case that one wants to run the same geometry at several different operating points, for example: for a range of advance ratios, there is another dispatch of the `analyze` function that takes in an input, `multipoint`, that is a vector of operating points.

```@docs; canonical=false
DuctAPE.analyze(multipoint::AbstractVector{TO},propulsor::Propulsor,options::Options) where TO<:OperatingPoint
```

If we were to continue the [tutorial](@ref "Run Analysis") and run a multi-point analysis on the example geometry given there, it might look something like this:

```julia
# - Advance Ratio Range - #
Js = range(0.0, 2.0; step=0.1)

# - Calculate Vinfs - #
D = 2.0 * rotorstator_parameters.Rtip[1] # rotor diameter
n = RPM / 60.0 # rotation rate in revolutions per second
Vinfs = Js * n * D

# - Set Operating Points - #
ops = [deepcopy(operating_point) for i in 1:length(Vinfs)]
for (iv, v) in enumerate(Vinfs)
    ops[iv].Vinf[] = v
end

# - Run Multi-point Analysis - #
outs_vec, success_flags = dt.analyze(ops, propulsor, dt.set_options(ops))
```

There are a few things to note here.
First, we want to make sure that the operating point objects we put into the input vector are unique instances.
Second, we need to use the dispatch of `set_options` that takes in the operating point vector to set up the right number of things in the background (like convergence flags for each operating point).
Third, the outputs of the analysis are vectors of the same outputs for a single analysis.


-----


## Pre-compiling the Caches

There are several available caches that can be precompiled to help speed up multiple analyses.
The first is a cache used for intermediate calculations in the pre- and post-processing phases of the analysis.
It can be preallocated using `allocate_prepost_container_cache`

```@docs; canonical=false
DuctAPE.allocate_prepost_container_cache
```

The second is a cache containing parameters used in the solver, in other words, the results of the pre-processing phase.
It can be preallocated using `allocate_solve_parameter_cache`.

```@docs; canonical=false
DuctAPE.allocate_solve_parameter_cache
```

The final precompileable cache is for intermediate calculations within the solve and can be preallocated using `allocate_solve_container_cache`

```@docs; canonical=false
DuctAPE.allocate_solve_container_cache
```


-----


## Circumventing the Automated Geometry Re-paneling

It is not advised to circument the automated geometry re-paneling, but if it must be done, the user needs to provide duct, centerbody, and wake nodes conforming to compatible geometry formatting.
The best use case for this is to use previously generated geometry or perhaps geometry exported from DFDC.

The process is not simple, but is possible.
You would have to manually run the dispatches of [`precompute_parameters`](@ref "DuctAPE.precompute_parameters") that take in the the repaneled body nodes and wake grid.
These dispatches exist for this purpose, but there is, by design, no convenience functions at this time to aid the user in easily bypassing the automated repaneling.
