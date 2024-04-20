# API Reference

```@contents
Pages = ["api.md"]
Depth = 5
```

## Public API

```@docs
DuctAPE.setup_analysis
```

```@docs
DuctAPE.DFDC_options
```

The remainder of the public API elements are found throughout the remainder of this documentation in the context of their usage.  The index at the bottom of this page may be helpful in locating them.

## Private API

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


### Analysis
```@docs
DuctAPE.analyze
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
