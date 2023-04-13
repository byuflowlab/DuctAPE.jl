#=

Sanity check using rotor only rotor and duct only duct geometry to see if the solutions together make sense.

=#

#---------------------------------#
#             Includes            #
#---------------------------------#
using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include("run_ccblade.jl")

# using Plots
# pyplot()
# using LaTeXStrings
include("../plots_default.jl")

#---------------------------------#
#         ROTOR Geometry          #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

# Blade Hub radius, in meters
Rhub = 0.10 * Rtip

# number of blades
B = 2

# x position of rotor
xrotor = 0.25

# Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
propgeom = [
    0.15 0.130 32.76
    0.20 0.149 37.19
    0.25 0.173 33.54
    0.30 0.189 29.25
    0.35 0.197 25.64
    0.40 0.201 22.54
    0.45 0.200 20.27
    0.50 0.194 18.46
    0.55 0.186 17.05
    0.60 0.174 15.97
    0.65 0.160 14.87
    0.70 0.145 14.09
    0.75 0.128 13.39
    0.80 0.112 12.84
    0.85 0.096 12.25
    0.90 0.081 11.37
    0.95 0.061 10.19
    # 1.00 0.041 8.99
]

# extract non-dimensional radial positions
rnondim = propgeom[:, 1]
# Dimensionalize chords
chords = propgeom[:, 2] * Rtip
# convert twists to radians
twists = propgeom[:, 3] * pi / 180

# use a NACA 4412 airfoils
#=
Note here we are using the CCBlade functionality to define the airfoils data function.
In addition, we are using the airfoils data file available from the CCBlade repository that has been extrapolated using the Viterna method as well as corrected for rotational effects as described in the CCBlade documentation.
=#
airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(r))

#---------------------------------#
#         Paneling Options        #
#---------------------------------#
#=
Note: the solver with interpolate the rotor data using the given number of blade element inputs
=#
nwake_sheets = 15

# non-dimensional wake length
wake_length = 1.0

# number of panels between discrete points
# in this case, 5 panels between rotors, 5 panels between last rotor and hub te, 3 panels between hub and duct te's, and 20 panels from duct TE to end of wake
npanels = [10; 3; 25]

nhub_inlet = 20
nduct_inlet = 20
#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

# freestream conditions
rho = 1.225 #kg/m^3
mu = 1.81e-5 # kg/(m⋅s)
asound = 341.0 #m/s

#---------------------------------#
#      Define BODY Coordinates    #
#---------------------------------#

# - Duct Coordinates - #
# use duct coordinates from FLOWFoil validation cases
#=
The file containing the duct coordinates contains the following geometry items:
- x_duct : x-coordinates of duct geometry defined from trailing edge to trailing edge clockwise
- r_duct : r-coordinates of duct geometry defined from trailing edge to trailing edge clockwise
note in this case, we include the radial offset of the duct since we have no rotor tip radius to define the duct radial location.
=#
include("../test/data/naca_662-015.jl")
duct_coordinates = [x_duct r_duct]./2.0

# - Hub Coordinates - #
# use hub coordinates from FLOWFoil validation cases
#=
Similarly, the hub coordinates here contain x_hub and r_hub, which contain the x and r coordinates from the leading edge to the trailing edge of the center body (hub), thus the coordinates are also effectively clockwise.
=#
include("../test/data/bodyofrevolutioncoords.jl")
hub_coordinates = [x_hub[1:(end - 1)]./2.0 r_hub[1:(end - 1)] * Rhub / maximum(r_hub)]

#---------------------------------#
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    xrotor, nwake_sheets, r=rnondim, chords, twists, airfoils, Rtip, Rhub, B, Omega
)]

# Paneling Parameters
paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

# Freestream Parameters
freestream = (; rho, mu, asound, Vinf=5.0)

######################################################################
#                                                                    #
#                       CHECK SEPARATE PIECES                        #
#                                                                    #
######################################################################

# initialize various inputs used in analysis
inputs = dt.precomputed_inputs(
    duct_coordinates, hub_coordinates, paneling_constants, rotor_parameters, freestream
)

### --- Sanity Plots --- ###
plot(; aspectratio=1)

for ib in 1:2
    plot!(
        inputs.body_panels[ib].panel_center[:, 1],
        inputs.body_panels[ib].panel_center[:, 2];
        color=mycolors[1],
        label="",
    )
end

plot!(
    inputs.rotor_source_panels[1].panel_center[:, 1],
    inputs.rotor_source_panels[1].panel_center[:, 2];
    color=mycolors[2],
    label="",
)

for iw in 1:nwake_sheets
    plot!(
        inputs.wake_vortex_panels[iw].panel_center[:, 1],
        inputs.wake_vortex_panels[iw].panel_center[:, 2];
        color=mycolors[3],
        label="",
    )
end

savefig("examples/all-together-geometry.pdf")

### --- Initialize States --- ###
states = dt.initialize_states(inputs)
gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)

dp = inputs.body_panels[1].panel_center[:, 1]
hp = inputs.body_panels[2].panel_center[:, 1]
gamd = gamb[1:length(dp)]
gamh = gamb[(length(dp) + 1):end]
plot(dp, gamd / 5.0; xlabel="x", ylabel="vs/vinf", label="initial duct suface velocity")
plot!(hp, gamh / 5.0; label="initial hub suface velocity")
savefig("examples/body-init-sanity-check.pdf")

#---------------------------------#
#           Run Analysis          #
#---------------------------------#
strengths = dt.analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream;
    tol=1e-8,
    maxiter=100,
)
