```@contents
Pages = ["private_postprocess.md"]
Depth = 5
```

```@docs
DuctAPE.post_process
```

## Velocities
```@docs
DuctAPE.get_body_tangential_velocities
DuctAPE.get_body_tangential_velocities!
DuctAPE.calculate_vtheta
DuctAPE.calculate_induced_velocities_on_bodywake
```

## Pressures
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

## Rotor Performance
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

## Boundary Layer

### Thermodynamics
```@docs
DuctAPE.sa1
DuctAPE.sa2
DuctAPE.standard_atmosphere
DuctAPE.ideal_gas_rho
DuctAPE.sutherlands_law
DuctAPE.speed_of_sound
DuctAPE.calculate_mach
DuctAPE.total_temperature
DuctAPE.total_pressure
DuctAPE.static_temperature
DuctAPE.static_pressure
DuctAPE.static_density
DuctAPE.convert_temperature_to_kelvin
DuctAPE.convert_viscosity
```

### General Boundary Layer Functions

```@docs
DuctAPE.arc_lengths_from_panel_lengths
DuctAPE.split_at_stagnation_point
DuctAPE.bl_step_fun
DuctAPE.set_boundary_layer_steps
DuctAPE.RK2
DuctAPE.RK4
```

### Head's Method Specific Functions

```@docs
DuctAPE.setup_boundary_layer_functions_head
DuctAPE.calculate_H
DuctAPE.calculate_cf
DuctAPE.boundary_layer_residual_head
DuctAPE.boundary_layer_residual_head!
DuctAPE.solve_head_boundary_layer!
```

### Green's Method Specific Functions

#### Initialization
```@docs
DuctAPE.setup_boundary_layer_functions_green
DuctAPE.calculate_radius_of_curvature
DuctAPE.d2_init
DuctAPE.H12bar_init
DuctAPE.CE_init
DuctAPE.initialize_turbulent_boundary_layer_states_green
```

#### Intermediate Calculations

```@docs
DuctAPE.Fc
DuctAPE.FR
DuctAPE.calculate_Re
DuctAPE.calculate_Cf0
DuctAPE.calculate_H12bar0
DuctAPE.calculate_Cf
DuctAPE.calculate_H12
DuctAPE.calculate_Ctau
DuctAPE.calculate_F
DuctAPE.calculate_H1
DuctAPE.calculate_dH12bardH1
DuctAPE.calculate_richardson_number
DuctAPE.calculate_beta
DuctAPE.longitudinal_curvature_influence
DuctAPE.lateral_strain_influence
DuctAPE.dilation_influence
DuctAPE.calculate_lambda
DuctAPE.calculate_d2dUedsUeeq0
DuctAPE.calculate_CEeq0
DuctAPE.calculate_Ctaueq0
DuctAPE.calculate_CEeq
DuctAPE.calculate_d2dUedsUeeq
```

#### Residual Functions

```@docs
DuctAPE.boundary_layer_residual_green
DuctAPE.boundary_layer_residual_green!
DuctAPE.solve_green_boundary_layer!
```

### Viscous Drag

```@docs
DuctAPE.squire_young
DuctAPE.total_viscous_drag_duct
DuctAPE.compute_viscous_drag_duct
DuctAPE.compute_single_side_drag_coefficient_head
DuctAPE.compute_single_side_drag_coefficient_green
DuctAPE.compute_viscous_drag_duct_schlichting
```
