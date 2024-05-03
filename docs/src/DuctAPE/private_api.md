# Private API

```@contents
Pages = ["private_api.md"]
Depth = 5
```

## Option Types
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

## Preprocess

### General

```@docs
DuctAPE.set_index_maps
DuctAPE.precompute_parameters
DuctAPE.precompute_parameters!
```

### Geometry
```@docs
DuctAPE.reinterpolate_geometry
DuctAPE.reinterpolate_geometry!
DuctAPE.generate_all_panels
DuctAPE.generate_all_panels!
```

#### Wake
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

#### Bodies
```@docs
DuctAPE.reinterpolate_bodies!
DuctAPE.place_duct!
```

#### Rotors
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


### Induced Velocities
```@docs
DuctAPE.calculate_unit_induced_velocities
DuctAPE.calculate_unit_induced_velocities!
DuctAPE.initialize_linear_system
DuctAPE.initialize_linear_system!
```

#### Unit Induced Velocities
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

#### Unit Induced Velocity Matrices
```@docs
DuctAPE.induced_velocities_from_vortex_panels_on_points
DuctAPE.induced_velocities_from_vortex_panels_on_points!
DuctAPE.induced_velocities_from_source_panels_on_points
DuctAPE.induced_velocities_from_source_panels_on_points!
DuctAPE.induced_velocities_from_trailing_edge_gap_panel!
```

#### Panel Method Velocity Functions
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

#### Quadrature

##### Integrands
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

##### Integrals
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



### State Initialization
```@docs
DuctAPE.initialize_velocities
DuctAPE.initialize_velocities!
DuctAPE.initialize_strengths!
```

## Analysis
```@docs
DuctAPE.analyze_multipoint
```

### Process
```@docs
DuctAPE.process
DuctAPE.solve
```

#### Residuals

##### CSOR
```@docs
DuctAPE.CSOR_residual!
DuctAPE.compute_CSOR_residual!
DuctAPE.relax_Gamr!
DuctAPE.relax_gamw!
DuctAPE.apply_relaxation_schedule
DuctAPE.update_CSOR_residual_values!
DuctAPE.check_CSOR_convergence!
```

##### External Solvers
```@docs
DuctAPE.system_residual
DuctAPE.system_residual!
DuctAPE.update_system_residual!
DuctAPE.estimate_states!
```

#### Solve Utilities
```@docs
DuctAPE.extract_initial_guess
DuctAPE.extract_state_variables
```

## Post-process

```@docs
DuctAPE.post_process
```

### Velocities
```@docs
DuctAPE.get_body_tangential_velocities
DuctAPE.get_body_tangential_velocities!
DuctAPE.calculate_vtheta
DuctAPE.calculate_induced_velocities_on_bodywake
```

### Pressures
```@docs
DuctAPE.steady_cp
DuctAPE.steady_cp!
DuctAPE.calculate_entropy_jumps
DuctAPE.calculate_rotor_jumps
DuctAPE.delta_cp
DuctAPE.calculate_body_delta_cp!
DuctAPE.calculate_bodywake_delta_cp
DuctAPE.get_body_cps
DuctAPE.get_body_cps!
DuctAPE.get_bodywake_cps
DuctAPE.forces_from_pressure
DuctAPE.forces_from_pressure!
DuctAPE.forces_from_TEpanels!
```

### Rotor Performance
```@docs
DuctAPE.inviscid_rotor_thrust
DuctAPE.inviscid_rotor_thrust!
DuctAPE.viscous_rotor_thrust
DuctAPE.viscous_rotor_thrust!
DuctAPE.inviscid_rotor_torque
DuctAPE.inviscid_rotor_torque!
DuctAPE.viscous_rotor_torque
DuctAPE.viscous_rotor_torque!
DuctAPE.rotor_power
DuctAPE.rotor_power!
DuctAPE.get_total_efficiency
DuctAPE.get_total_efficiency!
DuctAPE.get_induced_efficiency
DuctAPE.get_induced_efficiency!
DuctAPE.get_ideal_efficiency
DuctAPE.tqpcoeff
DuctAPE.tqpcoeff!
DuctAPE.get_blade_loads
DuctAPE.get_blade_loads!
```

## Utility Functions
```@docs
DuctAPE.promote_propulosor_type
DuctAPE.update_operating_point!
DuctAPE.isscalar
DuctAPE.dot
DuctAPE.norm
DuctAPE.cross2mag
DuctAPE.linear_transform
DuctAPE.extract_primals!
DuctAPE.lfs
DuctAPE.reset_containers!
DuctAPE.cache_dims!
DuctAPE.write_data
```
