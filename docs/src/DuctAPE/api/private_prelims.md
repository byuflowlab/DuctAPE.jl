## Option Types
```@docs
DuctAPE.DFDC_options
DuctAPE.ConvergenceType
DuctAPE.Relative
DuctAPE.Absolute
DuctAPE.SolverOptionsType
DuctAPE.InternalSolverOptions
DuctAPE.InternalPolyAlgorithmOptions
DuctAPE.ExternalSolverOptions
DuctAPE.ExternalPolyAlgorithmOptions
DuctAPE.GridSolverOptionsType
DuctAPE.IntegrationMethod
```

## Bookkeeping
```@docs
DuctAPE.get_problem_dimensions
```

## Caching

### Allocation

The following are various helper functions used in preallocating the various caches.

```@docs
DuctAPE.initialize_all_caches
DuctAPE.allocate_wake_panel_container!
DuctAPE.allocate_panel_containers!
DuctAPE.allocate_panel_container!
DuctAPE.allocate_body_panel_container!
DuctAPE.allocate_rotor_panel_container!
DuctAPE.allocate_solve_parameter_extras!
DuctAPE.allocate_grid_parameter_cache
DuctAPE.allocate_integration_containers
```

### Reshaping

The following are used internally to reshape the cache vectors into more usable formats.

```@docs
DuctAPE.withdraw_prepost_container_cache
DuctAPE.withdraw_solve_parameter_cache
DuctAPE.withdraw_solve_container_cache
DuctAPE.withdraw_grid_parameter_cache
```
