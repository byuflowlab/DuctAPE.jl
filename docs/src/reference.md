# API Reference

## Public


## Private

### Grid Functions

The wake grid initalization function actually calls two functions. The first is a method for initializing the grid geometry based on conservation of mass:

```@docs
DuctTAPE.generate_grid_points
```

The second is an elliptic grid solver using successive line over relaxation (SLOR) to relax the grid, meaning more accurately place the radial grid points along streamlines:

```@docs
DuctTAPE.relax_grid
```
