# Quick Start Guide

## Setting up Duct Geometry

For this example, we are going to set up an example available in the DFDC source files.
We have taken one of the case files and transcribed it into a julia file located in the [data directory](../../data/dfdc/dstestr2_case.jl).

```@setup geom
using DuctTAPE
using Plots
using Measures

default()
default(;
    fontfamily="Palatino Roman",
    size=(800, 600), #it appears that 100 ≈ 1inch in LaTeX
    fillalpha=0.125,
    fillcolor=RGB(128 / 255, 128 / 255, 128 / 255),
    linewidth=1.0,
    annotationfontfamily="Palatino Roman",
    markerstrokewidth=0.1,
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
    margin = 10mm,
)

plot(xlabel="x", ylabel="r", aspectratio=:equal)
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
plot!(innerwallx, innerwallr, linewidth=2)
plot!(outerwallx, outerwallr, linestyle=:dash, color=1, linewidth=2)
plot!(hubx, hubr, color=2, linewidth=2)
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

## Defining Rotors

Along with the duct geometry, we will need to define an array of `Rotor` objects.

```@docs
DuctTAPE.RotorGeometry
```

Note that we want to create an array, even if we only have one rotor.  When we initialize the grid, it will expect an array.
Also, our rotor object has more fields than are used in the original dfdc, for now, we'll set the section skew, rake, and solidity to nothing.
For the airfoils, we will use the [CCBlade.jl](https://flow.byu.edu/CCBlade.jl/stable/reference/#CCBlade.AlphaReAF) functionality and for our case here, define our airfoils as `ccb.AlphaReAF` objects (where CCBlade has been included in the DuctTAPE package and renamed 'ccb') using data files from digitized dfdc plots.

```@example geom
# -- GENERATE ROTOR OBJECT ARRAY

# set datapath for airfoil data files. Remember to grab these files and change the data path for your case.
datapath = "../../data/dfdc/airfoils/"

# generate airfoils for rotor
af1 = DuctTAPE.ccb.AlphaReAF([
    datapath * "disk1_re5e5.dat",
    datapath * "disk1_re1e6.dat",
    datapath * "disk1_re1.5e6.dat",
    datapath * "disk1_re2e6.dat",
])

# generate airfoils for stator
af2 = DuctTAPE.ccb.AlphaReAF([
    datapath * "disk2_re5e5.dat",
    datapath * "disk2_re1e6.dat",
    datapath * "disk2_re1.5e6.dat",
    datapath * "disk2_re2e6.dat",
])

#generate rotor object
rotor1 = DuctTAPE.RotorGeometry(
    xdisk1, #x position of rotor
    nblade1, #number of blades
    rnondim1, #radial stations
    0.0, #tip gap
    chord1, #chords
    beta1, #twists
    nothing, #skews
    nothing, #rakes
    fill(af1,length(rnondim1)), #airfoils
    nothing, #solidities
    rpm, #RPM
)

#generate stator object (rpm is zero for stator)
rotor2 = DuctTAPE.RotorGeometry(
    xdisk2, #x position of rotor
    nblade2, #number of blades
    rnondim2, #radial stations
    0.0, #tip gap
    chord2, #chords
    beta2, #twists
    nothing, #skews
    nothing, #rakes
    fill(af2,length(rnondim2)), #airfoils
    nothing, #solidities
    0.0, #RPM
)

#assemble array
rotors = [rotor1; rotor2]
```

Take note that the rotor objects are defined using non-dimensional parameters.
If we want to visualize things at this point, we'll have to dimensionalize things, which is done internally in the code using `BladeDimensions` objects.

```@docs
DuctTAPE.BladeDimensions
DuctTAPE.initialize_blade_dimensions
```

```@example geom
blade1 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor1)
blade2 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor2)

nr = length(blade1.rdim)

plot!(rotor1.xlocation.*ductgeometry.chord*ones(nr), blade1.rdim, color=4, linewidth=2)
plot!(rotor2.xlocation.*ductgeometry.chord*ones(nr), blade2.rdim, color=4, linewidth=2)

```

## Wake Grid Initialization

With the duct geometry and rotor data defined, we can now initialize the rotor wake grid.
We'll begin by setting some grid options.

!!! note
    The number of radial stations set for grid options must match the number of radial stations defined for the rotor objects.
    Furthermore, the rotors must have the same number of radial stations defined.


```@docs
DuctTAPE.GridOptions
DuctTAPE.defineGridOptions
```

```@example geom
# --- SET GRID OPTIONS
num_radial_stations = length(rnondim1)
grid_options = DuctTAPE.defineGridOptions(num_radial_stations)
```

With options set, rotor location chosen, and the wall and hub geometry available, we're finally ready to initialize the wake grid object.

```@docs
DuctTAPE.initialize_grid
DuctTAPE.WakeGridGeometry
```

```@example geom
# --- INITIALIZE GRID
wakegrid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

xg = wakegrid.x_grid_points
rg = wakegrid.r_grid_points
nx = wakegrid.nx
nr = wakegrid.nr

plot(; xlabel="x", ylabel="r", aspectratio=:equal)
plot!(
    ductgeometry.wallinnerxcoordinates,
    ductgeometry.wallinnerrcoordinates;
    color=1,
    linewidth=2,
)
plot!(
    ductgeometry.wallouterxcoordinates,
    ductgeometry.wallouterrcoordinates;
    color=1,
    linestyle=:dash,
    linewidth=2,
)
plot!(ductgeometry.hubxcoordinates, ductgeometry.hubrcoordinates; color=2, linewidth=2)

plot!(rotor1.xlocation .* ductgeometry.chord * ones(nr), blade1.rdim; color=4, linewidth=2)
plot!(rotor2.xlocation .* ductgeometry.chord * ones(nr), blade2.rdim; color=4, linewidth=2)

plot!(xg, rg; color=3, linewidth=0.5)
plot!(xg', rg'; color=3, linewidth=0.5)
```


!!! note
    This process should be relatively robust in that most combinations of wall and hub and rotor position should generate a grid.
    However, in the case that there is no overlap in the x-positions of the wall and hub, problems may occur later in the solver.
    In addition if the rotor is positioned somewhere not between the wall and hub, for example, if the rotor is out in front of the duct, then if a solution is found, it will likely be inaccurate.
    And finally, things will break if the hub and wall overlap in the radial direction, in other words, if the duct is blocked.

As can be seen in the plot above, there are a few features of the wake grid that deserve mention for the user's information.
First, the grid is designed such that the grid spacing in the axial direction will line up with the rotor positions as well as the trailing edge positions of the duct wall and hub (as well as the leading edge positions if they happen to be behind a rotor, which is not advised).
In addition, the grid axial spacing is taken to be as close as possible to the radial spacing, which is defined directly from the rotor radial station positions.
The wake spacing is started at the average of the axial spacing inside the duct area and then expanded by an expansion factor that can be defined by the user and is set to 1.1 by default.  This means that the end of the wake will actually not lie directly at the length input by the user (default 2x duct chord), but should be close enough.


## System Paneling

With the geometry defined and the grid initialized, we can define the system panels and their control points (centers).
To do so, we use the `generate_paneling` function:

```@docs
DuctTAPE.generate_paneling
```

For our case here, we are going to use a rough sampling of the geometry data in order to more clearly see the paneling.

```@setup geom
# --- DEFINE DUCT OBJECT
ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
    innerwallx[1:4:end],
    innerwallr[1:4:end],
    outerwallx[1:4:end],
    outerwallr[1:4:end],
    hubx[1:4:end],
    hubr[1:4:end]
)
# -- GENERATE ROTOR OBJECT ARRAY

#generate rotor object
rotor1 = DuctTAPE.RotorGeometry(
    xdisk1,
    nblade1,
    rnondim1[1:2:end],
    0.0,
    chord1[1:2:end],
    beta1[1:2:end],
    nothing,
    nothing,
    nothing,
    nothing,
    rpm,
)

#generate stator object (rpm is zero for stator)
rotor2 = DuctTAPE.RotorGeometry(
    xdisk2,
    nblade2,
    rnondim2[1:2:end],
    0.0,
    chord2[1:2:end],
    beta2[1:2:end],
    nothing,
    nothing,
    nothing,
    nothing,
    0.0,
)

#assemble array
rotors = [rotor1; rotor2]

# --- SET GRID OPTIONS
num_radial_stations = length(rnondim1[1:2:end])
grid_options = DuctTAPE.defineGridOptions(num_radial_stations)


# --- INITIALIZE GRID
wakegrid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

xg = wakegrid.x_grid_points
rg = wakegrid.r_grid_points
nx = wakegrid.nx
nr = wakegrid.nr

```



```@example geom
# Get paneling of various objects

wall_panels, hub_panels, wake_panels, rotor_source_panels = DuctTAPE.generate_paneling(
    ductgeometry, ductsplines, rotors, wakegrid
)
```

```@example geom

# PLOT PANELS

plot(; xlabel="x", ylabel="r", aspectratio=:equal, legend=true, label="")

# wall panels:
for i in 1:length(wall_panels.panel_edges_x)
    plot!(
        [wall_panels.panel_edges_x[i][1]; wall_panels.panel_edges_x[i][2]],
        [wall_panels.panel_edges_r[i][1]; wall_panels.panel_edges_r[i][2]];
        color=1,
        linewidth=0.5,
        markershape=:diamond,
        markersize=2,
        label="",
    )
end

scatter!(
    getindex.(wall_panels.panel_centers, 1),
    getindex.(wall_panels.panel_centers, 2);
    color=1,
    markersize=3,
    markershape=:circle,
    label="wall panel centers",
)

#hub panels:
for i in 1:length(hub_panels.panel_edges_x)
    plot!(
        [hub_panels.panel_edges_x[i][1]; hub_panels.panel_edges_x[i][2]],
        [hub_panels.panel_edges_r[i][1]; hub_panels.panel_edges_r[i][2]];
        markersize=2,
        markershape=:diamond,
        color=2,
        linewidth=0.5,
        label="",
    )
end

scatter!(
    getindex.(hub_panels.panel_centers, 1),
    getindex.(hub_panels.panel_centers, 2);
    markersize=3,
    markershape=:circle,
    color=2,
    label="hub panel centers",
)

#vortex sheet panels
for i in 1:length(wake_panels.panel_centers)
    plot!(
        [wake_panels.panel_edges_x[i][1]; wake_panels.panel_edges_x[i][2]],
        [wake_panels.panel_edges_r[i][1]; wake_panels.panel_edges_r[i][2]];
        markersize=2,
        markershape=:diamond,
        color=3,
        linewidth=0.5,
        label="",
    )
end

scatter!(
    getindex.(wake_panels.panel_centers, 1),
    getindex.(wake_panels.panel_centers, 2);
    markersize=3,
    markershape=:circle,
    color=3,
    label="vortex sheet panel centers",
)

#rotor source panels:
for i in 1:length(rotor_source_panels[1].panel_centers)
    plot!(
        [
            rotor_source_panels[1].panel_edges_x[i][1]
            rotor_source_panels[1].panel_edges_x[i][2]
        ],
        [
            rotor_source_panels[1].panel_edges_r[i][1]
            rotor_source_panels[1].panel_edges_r[i][2]
        ];
        markersize=2,
        markershape=:diamond,
        color=4,
        linewidth=0.5,
        label="",
    )
end

scatter!(
    getindex.(rotor_source_panels[1].panel_centers, 1),
    getindex.(rotor_source_panels[1].panel_centers, 2);
    markersize=3,
    markershape=:circle,
    color=4,
    label="rotor source panel centers",
)

# stator source panels
for i in 1:length(rotor_source_panels[2].panel_centers)
    plot!(
        [
            rotor_source_panels[2].panel_edges_x[i][1]
            rotor_source_panels[2].panel_edges_x[i][2]
        ],
        [
            rotor_source_panels[2].panel_edges_r[i][1]
            rotor_source_panels[2].panel_edges_r[i][2]
        ];
        markersize=2,
        markershape=:diamond,
        color=4,
        linewidth=0.5,
        label="",
    )
end

scatter!(
    getindex.(rotor_source_panels[2].panel_centers, 1),
    getindex.(rotor_source_panels[2].panel_centers, 2);
    markersize=3,
    markershape=:circle,
    color=4,
    label="",
)
```


The paneling for the rotor sources and the vortex wake sheets is based directly on the wake grid initialized in the previous section (remember that we used a rougher dataset here for visual clarity).
The paneling of the duct wall and hub are set such that the panels aft of the foremost rotor also align perfectly with the wake grid.
The points in front of the foremost rotor are set using cosine spacing such that the last panel before the foremost rotor is roughly similar in length to the average of the panel lenghts in the remainder of the duct.
(Note that the estimation process for this is not particularly robust at this point.)
