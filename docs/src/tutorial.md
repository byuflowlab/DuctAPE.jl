# Quick Start Guide

## Setting up Duct Geometry

You can define whatever geometry you'd like, and you could do so manually if you wanted.
For this example, however, we're going to use the airfoil geometry generation and manimulation functions available in FLOWFoil.jl.
In this example, we'll simply use the NACA 4-series airfoil function.

```@setup geom
using FLOWFoil
using DuctTAPE
using Plots
include("plots_default.jl")
plot(xlabel="x", ylabel="r")
```

For the wall geometry, we'll use a NACA 4420 airfoil that we'll flip over and then adjust the angle of attack and shift up from the center axis.
We'll also split the inner and outer portions of the duct wall for later use.
```@example geom
# --- WALL GEOMETRY DEFINITION
xwall, ywall = FLOWFoil.naca4(4, 4, 20)
wallcoords = [xwall -ywall]
wallangle = 8.0
walllocation = [0.0; 0.75]
FLOWFoil.position_coordinates(wallcoords, 1.0, wallangle, walllocation)
wallx = wallcoords[:, 1]
wallr = wallcoords[:, 2]
outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(wallx, wallr)
plot!(innerwallx,innerwallr,aspectratio=:equal)
plot!(outerwallx,outerwallr,linestyle=:dash, color=1)
```

Similarly for the hub geometry, we'll use a NACA 2410, but just use the top half of the airfoil and keep it at the center axis by using the `split_wall()` function:

```@docs
DuctTAPE.split_wall
```

We'll also scale it down and shift it back using FLOWFoil functionality.

```@example geom
# --- HUB GEOMETRY DEFINITION
_, hubxcoordinates, _, hubrcoordinates = FLOWFoil.naca4(2, 4, 10; split=true)
hubcoordinates = [hubxcoordinates hubrcoordinates]
FLOWFoil.position_coordinates(hubcoordinates, 0.8, 0.0, [0.25; 0.0])
hubx = hubcoordinates[:, 1]
hubr = hubcoordinates[:, 2]

plot!(hubx,hubr,color=2)
```

Up to this point, we have not used any functionality of the DuctTAPE package.  It should be remembered, therefore, that the user can define any wall and hub geometry they want, using whatever methods they want, as long as the geometry can be input into the DuctGeometry struct:

```@docs
DuctTAPE.DuctGeometry
DuctTAPE.defineDuctGeometry
```

Using the `defineDuctGeometry` contructor function, we can input our wall and hub geometries and let the leading and trailing edges and chord length be calculated automatically.

```@example geom
# --- DEFINE DUCT OBJECT
duct = DuctTAPE.defineDuctGeometry(
    innerwallx,
    innerwallr,
    outerwallx,
    outerwallr,
    hubcoordinates[:, 1],
    hubcoordinates[:, 2],
)
```

The `duct` object now contains all the geometry information we'll need.

## Wake Grid Initialization

For the sake of this example, we're going to define a rotor object that has a location only (we won't need the rest of the rotor defnition yet).

```@docs
DuctTAPE.Rotor
```

Note that we want to create an array, even if we only have one rotor.  When we initialize the grid, it will expect an array.
```@example geom
# --- GENERATE ROTOR OBJECT
rotor_x_location = 0.5
rotors = [DuctTAPE.Rotor(rotor_x_location, nothing, nothing, nothing, nothing, nothing)]
```

We can next define some grid options

```@docs
DuctTAPE.GridOptions
DuctTAPE.defineGridOptions
```

```@example geom
# --- SET GRID OPTIONS
num_radial_stations = 10
grid_options = DuctTAPE.defineGridOptions(num_radial_stations)
```

With options set, rotor location chosen, and the wall and hub geometry available, we're finally ready to initialize the wake grid object.

```@docs
DuctTAPE.initialize_grid
DuctTAPE.Grid
```

```@example geom
# --- INITIALIZE GRID
grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

xg = grid.x_grid_points
rg = grid.r_grid_points
nx = grid.nx
nr = grid.nr

plot!(xg, rg, color=3, linewidth=0.5)
plot!(xg', rg', color=3, linewidth=0.5)
```

!!! note
    This process should be relatively robust in that most combinations of wall and hub and rotor position should generate a grid.
    However, in the case that there is no overlap in the x-positions of the wall and hub, problems may occur later in the solver.
    In addition if the rotor is positioned somewhere not between the wall and hub, for example, if the rotor is out in front of the duct, then if a solution is found, it will likely be inaccurate.
    And finally, things will break if the hub and wall overlap in the radial direction, in other words, if the duct is blocked.
