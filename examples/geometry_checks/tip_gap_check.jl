#=

checking geometry using tip gap

=#

#---------------------------------#
#             Includes            #
#---------------------------------#
using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include("run_ccblade.jl")

# using Plots
# pyplot()
# using LaTeXStrings
include("../plots_default.jl")

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
# npanels = [5; 5; 4; 20]
npanels = [5; 4; 20] # no hub case

Rtip = 1.0

#! Tip Gap
#=
we have this dimensional such that it's easy for the user to know how the limits of the input (want to keep larger than 1e-4 or so to avoid near singularities between the wake and body panels)
=#
tip_gap = 1.0

# Blade Hub radius, in meters
Rhub = 0.10 * Rtip

# number of blades
B = 2

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
airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(rnondim))

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

#Vinf
Vinf = 5.0

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

# freestream conditions
rho = 1.225 #kg/m^3
mu = 1.81e-5 # kg/(m⋅s)
asound = 341.0 #m/s

#---------------------------------#
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters1 = (;
    xrotor=xrotors[1],
    nwake_sheets,
    r=rnondim,
    chords,
    twists,
    airfoils,
    Rtip,
    Rhub,
    tip_gap,
    B,
    Omega,
)
rotor_parameters2 = (; rotor_parameters1..., xrotor=xrotors[2])

rotor_parameters = [rotor_parameters1; rotor_parameters2]

# Paneling Parameters
paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

# Freestream Parameters
freestream = (; rho, mu, asound, Vinf)

######################################################################
#                                                                    #
#                       CHECK SEPARATE PIECES                        #
#                                                                    #
######################################################################

# initialize various inputs used in analysis
inputs = dt.precomputed_inputs(
    duct_coordinates,
    # hub_coordinates,
    nothing,
    paneling_constants,
    rotor_parameters,
    freestream;
    finterp=FLOWMath.linear,
)

### --- Sanity Plots --- ###
plot(; ylims=(0.0, 2.5), aspectratio=1)

# for ib in 1:2
for ib in 1:1
    plot!(
        inputs.body_panels[ib].panel_center[:, 1],
        inputs.body_panels[ib].panel_center[:, 2];
        color=mycolors[1],
        label="",
    )
end

for ir in 1:length(xrotors)
    plot!(
        inputs.rotor_source_panels[ir].panel_center[:, 1],
        inputs.rotor_source_panels[ir].panel_center[:, 2];
        color=mycolors[2],
        markershape=:square,
        markersize=1,
        label="",
    )
end

for iw in 1:nwake_sheets
    plot!(
        inputs.wake_vortex_panels[iw].panel_center[:, 1],
        inputs.wake_vortex_panels[iw].panel_center[:, 2];
        linewidth=0.5,
        color=mycolors[3],
        label="",
    )
end

savefig("examples/tip-gap-check-geometry.pdf")
