# API Reference

```@contents
Pages = ["api.md"]
Depth = 5
```

## Public API

```@docs
DuctAPE.setup_analysis
```

The remainder of the public API elements are found throughout the remainder of this documentation in the context of their usage.  The index at the bottom of this page may be helpful in locating them.

## Private API

### Bookkeeping
```@docs
DuctAPE.get_problem_dimensions
```

### Caching
```@docs
DuctAPE.allocate_wake_panel_container!
DuctAPE.allocate_panel_containers!
DuctAPE.allocate_panel_container!
DuctAPE.allocate_body_panel_container!
DuctAPE.allocate_rotor_panel_container!
DuctAPE.allocate_solve_parameter_extras!
```


### Analysis
```@docs
DuctAPE.analyze
```

### Post-process
```@docs
DuctAPE.post_process
```

## Index

```@index
Modules=[DuctAPE, DuctAPE.C4Blade]
```
