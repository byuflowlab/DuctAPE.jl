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

You may run all these simultaneously using the `initialize_all_caches` function.

```@docs; canonical=false
DuctAPE.initialize_all_caches
```

## How to pass the caches into an analysis

The precompiled caches can be passed in via keyword arguments to the analysis functions. If they are not, they are generated as the first step in the analysis.

```@docs; canonical=false
DuctAPE.analyze(
    propulsor::Propulsor,
    options::Options=set_options())
```