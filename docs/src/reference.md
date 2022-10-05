# API Reference

## Public

### Additional Airfoil Functions
These functions are nearly identical to CCBlade implementations, but include a solidity factor as one of the inputs.

```@autodocs
Modules = [DuctTAPE]
Pages = ["airfoils.jl"]
```

---

---

## Private

### Wake Grid Geometry Functions

The wake grid initalization function actually calls two functions. The first is a method for initializing the grid geometry based on conservation of mass:

```@docs
DuctTAPE.generate_grid_points
```

The second is an elliptic grid solver using successive line over relaxation (SLOR) to relax the grid, meaning more accurately place the radial grid points along streamlines:

```@docs
DuctTAPE.relax_grid
```

### Rotor Geometric Functions
```@docs
DuctTAPE.reinterpolate_rotor!
```

### Rotor & Wake Grid Aerodynamic Functions
```@docs
DuctTAPE.set_grid_aero!
DuctTAPE.set_rotor_velocities
DuctTAPE.calc_gamma_i
```

### Utility Functions

The following functions are used as various utility functions throughout the code.

```@autodocs
Modules = [DuctTAPE]
Pages   = ["utils.jl"]
```
