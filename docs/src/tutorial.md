# Quick Start Guide

## Setting up Duct Geometry

For this example, we are going to set up an example available in the DFDC source files.
We have taken one of the case files and transcribed it into a julia file located in the [data directory](../../data/dfdc/dstestr2_case.jl).

```@setup geom
using DuctTAPE
using Plots

default()
default(;
    fontfamily="Palatino Roman",
    size=(800, 600), #it appears that 100 ≈ 1inch in LaTeX
    fillalpha=0.125,
    fillcolor=RGB(128 / 255, 128 / 255, 128 / 255),
    linewidth=1.0,
    markerstrokealpha=0,
    annotationfontfamily="Palatino Roman",
    annotationfontsize=10,
    background_color_inside=nothing,
    background_color_legend=nothing,
    background_color_subplot=nothing,
    color_palette=[
        RGB(0.0, 46.0 / 255.0, 93.0 / 255.0), #BYU Blue
        RGB(155.0 / 255.0, 0.0, 0.0), #"BYU" Red
        RGB(128.0 / 255.0, 128.0 / 255.0, 128.0 / 255.0), #Middle Gray
        RGB(162.0 / 255.0, 227.0 / 255.0, 162.0 / 255.0), #Light Green
        RGB(243.0 / 255.0, 209.0 / 255.0, 243.0 / 255.0), #Pink
        RGB(205.0 / 255.0, 179.0 / 255.0, 0.0), #Yellow
        RGB(161.0 / 255.0, 161.0 / 255.0, 226.0 / 255.0), #Purple
    ],
    foreground_color_legend=nothing,
    legend=false, # include legend true/false
    grid=false, # background grid true/false
    gridlinewidth=0.5,
)

plot(xlabel="x", ylabel="r")
```

We'll go ahead and load that file and grab pieces of it as we go.

```@example geom
# --- WALL GEOMETRY DEFINITION
include("../../data/dfdc/dstestr2_case.jl");

```

The geometry is defined for a complete airfoil as the duct wall.
We actually want to split the duct wall coordinates for easier navigation as we set up the flow field.

```@docs
DuctTAPE.split_wall
```

```@example geom
# - Split Wall Coordinates
outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(ductx, ductr)

# - Plot Geometry
plot!(innerwallx, innerwallr, aspectratio=:equal)
plot!(outerwallx, outerwallr, linestyle=:dash, color=1)
plot!(hubx, hubr, color=2)
```


Now we want to put all the geometry together in a `DuctGeometry` object.

```@docs
DuctTAPE.DuctGeometry
DuctTAPE.defineDuctGeometry
```

Using the `defineDuctGeometry` contructor function, we can input our wall and hub geometries and let the leading and trailing edges and chord length be calculated automatically.
Note that this function also outputs a spline object for the inner duct wall and the hub wall.
These splines are used throughout the initialization process to help with rotor placement.

```@example geom
# --- DEFINE DUCT OBJECT
ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
    innerwallx,
    innerwallr,
    outerwallx,
    outerwallr,
    hubx,
    hubr
)
```

The `ductgeometry` and `ductsplines` objects now contains all the geometry information we'll need.

## Wake Grid Initialization

For the sake of this example, we're going to define a rotor object that has a location only (we won't need the rest of the rotor defnition yet).

```@docs
DuctTAPE.RotorGeometry
```

Note that we want to create an array, even if we only have one rotor.  When we initialize the grid, it will expect an array.
Also, our rotor object has more fields than are used in the original dfdc, for now, we'll set the section skew, rake, reynolds, airfoil, solidity, and mach to nothing.

```@example geom
# --- GENERATE ROTOR OBJECT ARRAY

#generate rotor object
rotor1 = DuctTAPE.RotorGeometry(xdisk1, nblade1, rnondim1, 0.0, chord1, beta1, nothing, nothing, nothing, nothing, rpm)

#generate stator object (rpm is zero for stator)
rotor2 = DuctTAPE.RotorGeometry(xdisk2, nblade2, rnondim2, 0.0, chord2, beta2, nothing, nothing, nothing, nothing, 0.0)

#assemble array
rotors = [rotor1; rotor2]
```

We can next define some grid options.

!!! note
    The number of radial stations set for grid options must match the number of radial stations defined for the rotor objects.
    Furthermore, the rotors must have the same number of radial stations defined.


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
grid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

xg = grid.x_grid_points
rg = grid.r_grid_points
nx = grid.nx
nr = grid.nr

xg, rg, nx, nr = DuctTAPE.generate_grid_points(ductgeometry, ductsplines, rotors, grid_options)

plot!(xg, rg, color=3, linewidth=0.5)
plot!(xg', rg', color=3, linewidth=0.5)
```


!!! note
    This process should be relatively robust in that most combinations of wall and hub and rotor position should generate a grid.
    However, in the case that there is no overlap in the x-positions of the wall and hub, problems may occur later in the solver.
    In addition if the rotor is positioned somewhere not between the wall and hub, for example, if the rotor is out in front of the duct, then if a solution is found, it will likely be inaccurate.
    And finally, things will break if the hub and wall overlap in the radial direction, in other words, if the duct is blocked.


## System Paneling


