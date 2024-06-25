## General

```@docs
DuctAPE.set_index_maps
DuctAPE.precompute_parameters
DuctAPE.precompute_parameters!
```

## Geometry
```@docs
DuctAPE.reinterpolate_geometry
DuctAPE.reinterpolate_geometry!
DuctAPE.generate_all_panels
DuctAPE.generate_all_panels!
```

### Wake
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

### Bodies
```@docs
DuctAPE.reinterpolate_bodies!
DuctAPE.place_duct!
```

### Rotors
```@docs
DuctAPE.interpolate_blade_elements
DuctAPE.interpolate_blade_elements!
DuctAPE.get_blade_ends_from_body_geometry
DuctAPE.get_blade_ends_from_body_geometry!
DuctAPE.get_local_solidity
DuctAPE.get_stagger
DuctAPE.generate_rotor_panels
DuctAPE.generate_rotor_panels!
```


## Induced Velocities
```@docs
DuctAPE.calculate_unit_induced_velocities
DuctAPE.calculate_unit_induced_velocities!
DuctAPE.initialize_linear_system
DuctAPE.initialize_linear_system!
```

### Unit Induced Velocities
```@docs
DuctAPE.calculate_xrm
DuctAPE.calculate_xrm!
DuctAPE.get_elliptics
DuctAPE.vortex_ring_vz
DuctAPE.vortex_ring_vz!
DuctAPE.smoke_ring_vz
DuctAPE.vortex_ring_vr
DuctAPE.vortex_ring_vr!
DuctAPE.source_ring_vz
DuctAPE.source_ring_vz!
DuctAPE.source_ring_vr
DuctAPE.source_ring_vr!
```

### Unit Induced Velocity Matrices
```@docs
DuctAPE.induced_velocities_from_vortex_panels_on_points
DuctAPE.induced_velocities_from_vortex_panels_on_points!
DuctAPE.induced_velocities_from_source_panels_on_points
DuctAPE.induced_velocities_from_source_panels_on_points!
DuctAPE.induced_velocities_from_trailing_edge_gap_panel!
```

### Panel Method Velocity Functions
```@docs
DuctAPE.vortex_aic_boundary_on_boundary
DuctAPE.vortex_aic_boundary_on_boundary!
DuctAPE.vortex_aic_boundary_on_field
DuctAPE.vortex_aic_boundary_on_field!
DuctAPE.add_kutta!
DuctAPE.add_te_gap_aic!
DuctAPE.source_aic
DuctAPE.source_aic!
DuctAPE.freestream_influence_vector
DuctAPE.freestream_influence_vector!
DuctAPE.assemble_lhs_matrix
DuctAPE.assemble_lhs_matrix!
DuctAPE.factorize_LHS
DuctAPE.factorize_LHS!
DuctAPE.assemble_rhs_matrix
DuctAPE.assemble_rhs_matrix!
DuctAPE.calculate_normal_velocity
DuctAPE.calculate_normal_velocity!
```

### Quadrature

#### Integrands
```@docs
DuctAPE.nominal_vortex_induced_velocity_sample
DuctAPE.nominal_vortex_induced_velocity_sample!
DuctAPE.subtracted_singular_vortex_influence
DuctAPE.subtracted_singular_vortex_influence!
DuctAPE.analytically_integrated_vortex_influence
DuctAPE.analytically_integrated_vortex_influence!
DuctAPE.self_vortex_induced_velocity_sample
DuctAPE.self_vortex_induced_velocity_sample!
DuctAPE.nominal_source_induced_velocity_sample
DuctAPE.nominal_source_induced_velocity_sample!
DuctAPE.subtracted_singular_source_influence
DuctAPE.subtracted_singular_source_influence!
DuctAPE.analytically_integrated_source_influence
DuctAPE.analytically_integrated_source_influence!
DuctAPE.self_source_induced_velocity_sample
DuctAPE.self_source_induced_velocity_sample!
```

#### Integrals
```@docs
DuctAPE.nominal_vortex_panel_integration
DuctAPE.nominal_vortex_panel_integration!
DuctAPE.self_vortex_panel_integration
DuctAPE.self_vortex_panel_integration!
DuctAPE.nominal_source_panel_integration
DuctAPE.nominal_source_panel_integration!
DuctAPE.self_source_panel_integration
DuctAPE.self_source_panel_integration!
DuctAPE.extrapolate!
```



## State Initialization
```@docs
DuctAPE.initialize_velocities
DuctAPE.initialize_velocities!
DuctAPE.initialize_strengths!
```
