# Advanced Option Selection

DuctAPE has been written in an attempt to make as many of the available options exposed to the user as possible.  This means that there are quite a few options to select from if not using the option convenience functions.
To help the user, the majority of overarching option types are defined using the `@kwdef` macro and have default values that should be reasonable in most cases.
We will introduce some of the available options here that may be of common interest.

## Quadrature

There are several implementations for different quadrature approaches depending on user desires; they include:
- [Gauss-Legendre quadature](@ref "DuctAPE.GaussLegendre") (default),
- [Gauss-Kronrod Quadrature](@ref "DuctAPE.GaussKronrod"), and
- [Romberg Quadrature](@ref "DuctAPE.Romberg") methods.

The docstrings for each are linked to the names.

## Elliptic Grid Solvers

As part of the pre-process, an elliptic grid defining the wake geometry is solved with a system of Poisson equations.
For this solve there currently two options:

- [SLOR](@ref "DuctAPE.SLORGridSolverOptions"): DFDC grid solver
- [SLOR+Newton](@ref "DuctAPE.GridSolverOptions")

The SLOR (successive line over relaxation) is the method employed by DFDC, and can be used by itself, or as a preconditioner to a Newton solve (using NLsolve.jl).

## Aerodynamics Solvers

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

## Other Options

The remaining option details can be found in the [Options](@ref "DuctAPE.Options") docstring.
