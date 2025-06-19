## Precompiled Caches

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

You may run all these simultaneously using the `initialize_all_caches` function.

```@docs; canonical=false
DuctAPE.initialize_all_caches
```

As an example of how to run this function, we'll grab [solver options](@ref "Aerodynamics Solver Options") and [paneling constants](@ref "Paneling Constants") from previous examples

```@example cache
using DuctAPE
using LineSearches

# - grab an object of SolverOptionsType defined in a previous example - #
aero_solver_options = DuctAPE.NLsolveOptions(;
    algorithm=:newton,
    atol=1e-10,
    iteration_limit=30,
    linesearch_method=LineSearches.BackTracking, #don't include parentheses on method handle
    linesearch_kwargs=(; order=3, maxstep=1e6),
)

# - grab an object of PanelingConstants type from the Getting Started tutorial - #
num_duct_inlet_panels = 30
num_center_body_inlet_panels = 30
num_panels = [30, 1, 30]
dte_minus_cbte = -1.0
num_wake_sheets = 11
wake_length = 0.8

# assemble paneling constants
paneling_constants = DuctAPE.PanelingConstants(
    num_duct_inlet_panels,
    num_center_body_inlet_panels,
    num_panels,
    dte_minus_cbte,
    num_wake_sheets,
    wake_length,
)

# - Airfoils are required for proper sizing of the caches - #
# DFDC-type airfoil object
afparams = DuctAPE.c4b.DFDCairfoil()

# specify the airfoil array
airfoils = [fill(afparams, 6)]

# - Initialize Caches - #
prepost_container_caching, solve_parameter_caching, solve_container_caching = DuctAPE.initialize_all_caches(
    aero_solver_options, paneling_constants, airfoils
)

```

## How to pass the caches into an analysis

The precompiled caches can be passed in via keyword arguments to the analysis functions. If they are not, they are generated as the first step in the analysis.

