using DuctTAPE
const dt = DuctTAPE
using CCBlade
const ccb = CCBlade
using FLOWMath
const fm = FLOWMath
include("../plots_default.jl")

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

#---------------------------------#
#             ROTORS              #
#---------------------------------#
# Required rotor information for duct and wake generation
xrotors = [0.5; 0.75]

airfoil = fill(ccb.AlphaAF("test/data/naca4412.dat"), 2)
# rotor parameters
rotor1_parameters = (;
    Rtip=1.5,
    B=2,
    Omega=50,
    xrotor=xrotors[1],
    r=[0.0, 1.0],
    chords=[0.5, 0.25],
    twists=[50.0, 10.0],
    airfoils=airfoil,
)

# stator parameters
rotor2_parameters = (; rotor1_parameters..., xrotor=xrotors[2])

# array with rotor and stator parameters
# rotor_parameters = [rotor1_parameters]#, rotor2_parameters]
rotor_parameters = [rotor1_parameters, rotor2_parameters]

nwake_sheets = 10

paneling_constants = (
    wake_length=1.0,
    nwake_sheets=nwake_sheets,
    nhub_inlet=4,
    nduct_inlet=5,
    # npanels=[5; 4; 20],
    npanels=[5; 5; 4; 20], #this is a vector of number of panels between discrete points after the first rotor, e.g. between the rotor and the duct trailing edge and between the duct trailing edge and the end of the wake
)

inputs = dt.precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters, #vector of named tuples
    (Vinf=5.0,);#freestream;
    finterp=fm.linear,
)

# - Plot Geometry - #
plot(; aspectratio=:equal, xlabel="x", ylabel="r")
plot!(
    inputs.t_duct_coordinates[:, 1],
    inputs.t_duct_coordinates[:, 2];
    marker=true,
    markersize=3,
    label="Bodies",
    color=1,
)
plot!(
    inputs.rp_hub_coordinates[:, 1],
    inputs.rp_hub_coordinates[:, 2];
    marker=true,
    markersize=3,
    label="",
    color=1,
)
for i in 1:(inputs.num_rotors)
    plot!(
        rotor_parameters[i].xrotor .* ones(length(inputs.rotor_panel_centers[:, i])),
        inputs.rotor_panel_centers[:, i];
        label="RPC $i",
        color=2,
        seriestype=:scatter,
    )
    plot!(
        rotor_parameters[i].xrotor .* ones(length(inputs.rotor_panel_edges[:, i])),
        inputs.rotor_panel_edges[:, i];
        label="RPE $i",
        color=3,
        seriestype=:scatter,
    )
end

for i in 1:nwake_sheets
    mylabel = i == 1 ? "WPC" : ""
    plot!(
        inputs.wake_vortex_panels[i].panel_center[:, 1],
        inputs.wake_vortex_panels[i].panel_center[:, 2];
        seriestype=:scatter,
        color=4,
        markersize=2,
        label=mylabel,
    )
end

savefig("examples/initialized-geometry.pdf")

