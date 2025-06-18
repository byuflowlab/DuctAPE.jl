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

As an example of how to run this function, we'll grab [solver options](@ref "Aerodynamics Solvers") and [paneling constants](@ref "Paneling Constants") from previous examples

```@julia
# - grab an object of SolverOptionsType defined in a previous example - #
aero_solver_options = DuctAPE.NLsolveOptions(;
    algorithm=:newton,
    atol=1e-10,
    iteration_limite=30,
    linesearch_method=LineSearches.BackTracking, #don't include parentheses on method handle
    linesearch_kwargs=(; order=3, maxstep=1e6),
    additional_kwargs=(; autoscale=false),
)

# - grab an object of PanelingConstants type from the Getting Started tutorial - #
# number of panels for the duct inlet
nduct_inlet = 30

# number of panels for the center body inlet
ncenter_body_inlet = 30

# number of panels from:
#  - rotor to duct trailing edge
#  - duct trailing edge to center body trailing edge
#  - center body trailing edge to end of wake
npanels = [30, 1, 30]

# the duct trailing edge is ahead of the center_body trailing edge.
dte_minus_cbte = -1.0

# number of wake sheets (one more than blade elements to use)
nwake_sheets = 11

# non-dimensional wake length aft of rear-most trailing edge
wake_length = 0.8

# assemble paneling constants
paneling_constants = DuctAPE.PanelingConstants(
    nduct_inlet, ncenter_body_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
)

# DFDC-type airfoil object
afparams = DuctAPE.c4b.DFDCairfoil(;
    alpha0=0.0,
    clmax=1.5,
    clmin=-1.0,
    dclda=6.28,
    dclda_stall=0.5,
    dcl_stall=0.2,
    cdmin=0.012,
    clcdmin=0.1,
    dcddcl2=0.005,
    cmcon=0.0,
    Re_ref=2e5,
    Re_exp=0.35,
    mcrit=0.7,
)

# specify the airfoil array
airfoils = [fill(afparams, length(r))]

# - Initialize Caches - #
prepost_container_caching, solve_parameter_caching, solve_container_caching = DuctAPE.initialize_all_caches(aero_solver_options, paneling_constants, airfoils)
```

## How to pass the caches into an analysis

The precompiled caches can be passed in via keyword arguments to the analysis functions. If they are not, they are generated as the first step in the analysis.

```@docs; canonical=false
DuctAPE.analyze(
    ducted_rotor::DuctedRotor,
    operating_point::OperatingPoint,
    reference_parameters::ReferenceParameters,
    options::Options=set_options())
```
