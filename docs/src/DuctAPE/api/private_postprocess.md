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
DuctAPE.total_temperature
DuctAPE.total_pressure
DuctAPE.static_temperure
DuctAPE.static_pressure
DuctAPE.static_density
```

### Boundary Layer Functions

#### Set Up
```@docs
DuctAPE.arc_lengths_from_panel_lengths
DuctAPE.split_at_stagnation_point
DuctAPE.bl_step_fun
DuctAPE.set_boundary_layer_steps
DuctAPE.calculate_radius_of_curvature
DuctAPE.setup_boundary_layer_functions
```

#### Initialization
```@docs
DuctAPE.d2_init
DuctAPE.H12bar_init
DuctAPE.CE_init
DuctAPE.initialize_turbulent_boundary_layer_states
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
DuctAPE.boundary_layer_residual
DuctAPE.boundary_layer_residual!
```

#### Solver Functions

```@docs
DuctAPE.RK2
DuctAPE.RK4
DuctAPE.solve_turbulent_boundary_layer_rk!
```

#### Viscous Drag

```@docs
DuctAPE.squire_young
DuctAPE.total_viscous_drag_duct
DuctAPE.compute_single_side_drag_coefficient
DuctAPE.compute_viscous_drag_duct
```
