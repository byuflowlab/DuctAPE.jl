# Public API

```@contents
Pages = ["public_api.md"]
Depth = 5
```

## Input Types
```@docs
DuctAPE.Propulsor
DuctAPE.RotorStatorParameters
DuctAPE.OperatingPoint
DuctAPE.PanelingConstants
DuctAPE.ReferenceParameters
```

## Preallocations
```@docs
DuctAPE.allocate_prepost_container_cache
DuctAPE.allocate_solve_parameter_cache
DuctAPE.allocate_solve_container_cache
```

## Options

### General Options
```@docs
DuctAPE.Options
DuctAPE.set_options
```

### Integration Options
```@docs
DuctAPE.IntegrationOptions
DuctAPE.GaussLegendre
DuctAPE.GaussKronrod
DuctAPE.Romberg
```

### Solver Options

#### Elliptic Grid Solve
```@docs
DuctAPE.SLORGridSolverOptions
DuctAPE.GridSolverOptions
```

#### Aerodynamics Solve
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
```
## Preprocess

```@docs
DuctAPE.setup_analysis
```

## Analysis
```@docs
DuctAPE.analyze
```

## Miscellaneous

### Airfoil/Geometry Manipulation

### NACA 6-Series Cascade Geometry Generation

