# API Reference

```@contents
Pages = ["api.md"]
Depth = 5
```

## Public API

### Input Types
```@docs
DuctAPE.Propulsor
DuctAPE.RotorStatorParameters
DuctAPE.OperatingPoint
DuctAPE.PanelingConstants
DuctAPE.ReferenceParameters
```

### Options

#### General Options
```@docs
DuctAPE.Options
DuctAPE.set_options
```

#### Integration Options
```@docs
DuctAPE.IntegrationOptions
DuctAPE.GaussLegendre
DuctAPE.GaussKronrod
DuctAPE.Romberg
```

#### Solver Options

##### Elliptic Grid Solve
```@docs
DuctAPE.SLORGridSolverOptions
DuctAPE.GridSolverOptions
```

##### Aero Solve
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

### Analysis
```@docs
DuctAPE.analyze
```

### Preprocess

```@docs
DuctAPE.setup_analysis
```
### Process

### Postprocess





----------------------------------------------------


## Private API

### Option Types
```@docs
DuctAPE.DFDC_options
DuctAPE.ConvergenceType
DuctAPE.Relative
DuctAPE.Absolute
DuctAPE.SolverOptionsType
DuctAPE.ExternalSolverOptions
DuctAPE.PolyAlgorithmOptions
DuctAPE.GridSolverOptionsType
DuctAPE.IntegrationMethod
```

### Bookkeeping
```@docs
DuctAPE.get_problem_dimensions
```

### Caching

#### Allocation

The following are various helper functions used in preallocating the various caches.

```@docs
DuctAPE.allocate_wake_panel_container!
DuctAPE.allocate_panel_containers!
DuctAPE.allocate_panel_container!
DuctAPE.allocate_body_panel_container!
DuctAPE.allocate_rotor_panel_container!
DuctAPE.allocate_solve_parameter_extras!
```

#### Reshaping

The following are used internally to reshape the cache vectors into more usable formats.

```@docs
DuctAPE.withdraw_prepost_container_cache
DuctAPE.withdraw_solve_parameter_cache
DuctAPE.withdraw_solve_container_cache
```

### Preprocess

#### General

```@docs
DuctAPE.set_index_maps
DuctAPE.precompute_parameters
DuctAPE.precompute_parameters!
```

#### Geometry
```@docs
DuctAPE.reinterpolate_geometry
DuctAPE.reinterpolate_geometry!
DuctAPE.generate_all_panels
DuctAPE.generate_all_panels!
```

##### Wake
```@docs
DuctAPE.discretize_wake
DuctAPE.generate_wake_grid
DuctAPE.generate_wake_grid!
DuctAPE.initialize_wake_grid
DuctAPE.initialize_wake_grid!
DuctAPE.relax_grid!
DuctAPE.generate_wake_panels
DuctAPE.generate_wake_panels!
DuctAPE.get_wake_k
DuctAPE.get_wake_k!
```

##### Bodies
```@docs
DuctAPE.reinterpolate_bodies!
DuctAPE.place_duct!
```

##### Rotors
```@docs
DuctAPE.interpolate_blade_elements
DuctAPE.interpolate_blade_elements!
```


#### Induced Velocities
```@docs
DuctAPE.calculate_unit_induced_velocities
DuctAPE.calculate_unit_induced_velocities!
DuctAPE.initialize_linear_system
DuctAPE.initialize_linear_system!
```

#### State Initialization



TODO: still need to do everything in initializestates as well as everything in the geometry and velocities directories

### Analysis
```@docs
DuctAPE.analyze_multipoint
```

### Post-process
```@docs
DuctAPE.post_process
```

### Convenience Functions
```@docs
DuctAPE.promote_propulosor_type
```

## Index

```@index
Modules=[DuctAPE, DuctAPE.C4Blade]
```
