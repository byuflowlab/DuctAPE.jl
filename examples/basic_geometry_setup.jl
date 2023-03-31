using DuctTAPE
using FLOWMath
include("../plots_default.jl")
# pyplot()
gr()

#---------------------------------#
#      Duct and Hub Geometry      #
#---------------------------------#
# duct geometry
duct_x = [1.25; 0.375; 0.0; 0.375; 1.25]
duct_r = 0.75 .+ [1.0; 0.875; 1.0; 1.375; 1.0]
duct_coordinates = [duct_x duct_r]

# hub geometry
hub_x = [0.0; 0.35; 1.0]
hub_r = [0.0; 0.25; 0.0]
hub_coordinates = [hub_x hub_r]

# Required rotor information for duct and wake generation
xrotors = [0.5; 0.75]
Rtip = 1.0 # leading rotor tip radius

# non-dimensional wake length
wake_length = 1.0

# number of wakes for each rotor (synonymous with number of blade elements)
nwake_sheets = 10

# number of panels between discrete points
# in this case, 5 panels between rotors, 5 panels between last rotor and hub te, 3 panels between hub and duct te's, and 20 panels from duct TE to end of wake
npanels = [5; 5; 4; 20]

# discretize the wake x-coordinates
xwake, rotor_indices = DuctTAPE.discretize_wake(
    duct_coordinates, hub_coordinates, xrotors, wake_length, nwake_sheets, npanels
)

# number of panels between hub leading edge and first rotor
nhub_inlet = 5

# number of panels between duct leading edge and first rotor
nduct_inlet = 5

# number of panels on duct outer surface
nduct_outer = 10

# update the body paneling to match the wake discretization
new_duct_xr, new_hub_xr = DuctTAPE.update_body_geometry(
    duct_coordinates,
    hub_coordinates,
    xwake,
    nhub_inlet,
    nduct_inlet,
    nduct_outer;
    finterp=FLOWMath.linear,
)

# shift the duct geometry according to the leading rotor tip radius, and return the rotor hub and tip dimensions for all rotors
trans_duct_xr, Rtips, Rhubs = DuctTAPE.place_duct(new_duct_xr, new_hub_xr, Rtip, xrotors)

# generate the body panels
body_panels = DuctTAPE.generate_body_panels(trans_duct_xr, new_hub_xr)

## -- PLOT -- ##
plot(duct_x, duct_r; label="input geometry", color=mycolors[1], marker=true, linestyle=:dot)
plot!(hub_x, hub_r; label="", color=mycolors[1], marker=true, linestyle=:dot)
plot!(
    xwake,
    -0.0625 * ones(length(xwake));
    markershape=:vline,
    label="wake x-coordinates",
)
plot!(
    trans_duct_xr[:, 1],
    trans_duct_xr[:, 2];
    color=mycolors[2],
    markershape=:vline,
    label="re-paneled, shifted bodies",
)
plot!(new_hub_xr[:, 1], new_hub_xr[:, 2]; color=mycolors[2], markershape=:vline, label="")

for i in 1:length(body_panels)
    mylabel = i == 1 ? "body panel centers" : ""
    plot!(
        body_panels[i].panel_center[:, 1],
        body_panels[i].panel_center[:, 2];
        seriestype=:scatter,
        markersize=3,
        markershape=:square,
        color=mycolors[1],
        label=mylabel,
    )
end

#---------------------------------#
#              WAKE               #
#---------------------------------#
# get discretization of wakes at leading rotor position
rwake = range(Rhubs[1], Rtips[1], nwake_sheets)

# initialize wake grid
xgrid, rgrid = DuctTAPE.initialize_wake_grid(trans_duct_xr, new_hub_xr, xwake, rwake)
for i in 1:length(rwake)
    mylabel = i == 1 ? "initial wake panel edges" : ""
    plot!(xgrid[:, i], rgrid[:, i]; markershape=:vline, color=mycolors[5], label=mylabel)
end

# Relax Grid
DuctTAPE.relax_grid!(xgrid, rgrid; max_iterations=100, tol=1e-9, verbose=false)
for i in 1:length(rwake)
    mylabel = i == 1 ? "relaxed wake panel edges" : ""
    plot!(xgrid[:, i], rgrid[:, i]; markershape=:vline, color=mycolors[4], label=mylabel)
end

# create wake panels
wake_panels = DuctTAPE.generate_wake_panels(xgrid, rgrid)

## -- PLOT -- ##
for i in 1:length(wake_panels)
    mylabel = i == 1 ? "wake panel centers" : ""
    plot!(
        wake_panels[i].panel_center[:, 1],
        wake_panels[i].panel_center[:, 2];
        seriestype=:scatter,
        markersize=2,
        markershape=:square,
        color=mycolors[7],
        label=mylabel,
    )
end

#---------------------------------#
#             ROTORS              #
#---------------------------------#

# rotor parameters
rotor1_parameters = (;
    B=2,
    omega=50,
    xrotor=xrotors[1],
    rblade=[0.0, 1.0],
    chords=[0.5, 0.25],
    twists=[50.0, 10.0],
    airfoils=[nothing, nothing],
)

# stator parameters
rotor2_parameters = (; rotor1_parameters..., xrotor=xrotors[2])

# array with rotor and stator parameters
rotor_parameters = [rotor1_parameters, rotor2_parameters]

# generate rotor source panel objects
rotor_source_panels = [
    DuctTAPE.generate_rotor_panels(xrotors[i], rgrid[rotor_indices[i], :]) for
    i in 1:length(xrotors)
]

## -- PLOT -- ##
plot!(xrotors[1] * ones(2), [Rhubs[1]; Rtips[1]]; color=:black, label="rotor locations")
plot!(xrotors[2] * ones(2), [Rhubs[2]; Rtips[2]]; color=:black, label="")
plot!(
    rotor_source_panels[1].panel_center[:, 1],
    rotor_source_panels[1].panel_center[:, 2];
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=mycolors[6],
    label="rotor blade element locations",
)
plot!(
    rotor_source_panels[2].panel_center[:, 1],
    rotor_source_panels[2].panel_center[:, 2];
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=mycolors[6],
    label="",
)

# blade elements for rotors
blade_elements = [
    DuctTAPE.generate_blade_elements(
        rotor_parameters[i].B,
        rotor_parameters[i].omega,
        rotor_parameters[i].xrotor,
        rotor_parameters[i].rblade,
        rotor_parameters[i].chords,
        rotor_parameters[i].twists,
        rotor_parameters[i].airfoils,
        Rtips[i],
        Rhubs[i],
        rotor_source_panels[i].panel_center[:, 2],
    ) for i in 1:length(xrotors)
]

savefig("examples/basic_geometry_setup.pdf")
