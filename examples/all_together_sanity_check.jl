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
include("rotor_only/run_ccblade.jl")

using FLOWMath

# using Plots
# pyplot()
# using LaTeXStrings
include("../plots_default.jl")

#---------------------------------#
#         ROTOR Geometry          #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

#! Tip Gap
#=
we have this dimensional such that it's easy for the user to know how the limits of the input (want to keep larger than 1e-4 or so to avoid near singularities between the wake and body panels)
=#
# tip_gap = 0.0
# tip_gap = 0.1*Rtip #10 percent tip gap
tip_gap = 100.0

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
airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(rnondim))

#---------------------------------#
#         Paneling Options        #
#---------------------------------#
#=
Note: the solver with interpolate the rotor data using the given number of blade element inputs
=#
nwake_sheets = 18

# non-dimensional wake length
wake_length = 1.0

# number of panels between discrete points
# in this case, 5 panels between rotors, 5 panels between last rotor and hub te, 3 panels between hub and duct te's, and 20 panels from duct TE to end of wake
# npanels = [10; 5; 25]
npanels = [40; 80]

nhub_inlet = 40
nduct_inlet = 40
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
duct_coordinates = [x_duct r_duct] ./ 2.0

# - Hub Coordinates - #
# use hub coordinates from FLOWFoil validation cases
#=
Similarly, the hub coordinates here contain x_hub and r_hub, which contain the x and r coordinates from the leading edge to the trailing edge of the center body (hub), thus the coordinates are also effectively clockwise.
=#
# include("../test/data/bodyofrevolutioncoords.jl")
# hub_coordinates = [x_hub[1:(end - 1)] ./ 2.0 r_hub[1:(end - 1)] * Rhub / maximum(r_hub)]
hub_coordinates = nothing

#---------------------------------#
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    xrotor, nwake_sheets, r=rnondim, chords, twists, airfoils, Rtip, Rhub, tip_gap, B, Omega
)]

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
    # finterp=FLOWMath.linear,
)

### --- Sanity Plots --- ###
plot(; aspectratio=1)

# for ib in 1:2
for ib in 1:1
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
        linewidth=0.5,
        color=mycolors[3],
        label="",
    )
end

savefig("examples/all-together-geometry.pdf")

### --- Initialize States --- ###
states = dt.initialize_states(inputs)
gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)

## -- check body surface velocity initialiation -- ##
dp = inputs.body_panels[1].panel_center[:, 1]
# hp = inputs.body_panels[2].panel_center[:, 1]
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

pb = plot(dp, gamd; xlabel="x", ylabel="cp", label="initial duct surface pressure")
# plot!(pb, hp, gamh ; label="initial hub surface pressure")
# savefig(pb, "examples/body-init-sanity-check.pdf")

## -- check rotor circulation and source initial strengths -- ##
pG = plot(Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="initial")
# savefig(pG, "examples/rotorcirculation-init-sanity-check.pdf")

ps = plot(sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="initial")
# savefig(ps, "examples/rotorsources-init-sanity-check.pdf")

pw = plot(
    gamw, inputs.rotor_panel_edges; xlabel=L"\gamma_\theta", ylabel="r", label="initial"
)
# savefig(pw, "examples/wakegamma-init-sanity-check.pdf")

#---------------------------------#
#           Run Analysis          #
#---------------------------------#
ccbouts = run_ccblade(Vinf)

for i in 1:3
    # - Extract commonly used items from precomputed inputs - #
    blade_elements = inputs.blade_elements
    rpc = inputs.rotor_panel_centers
    Vinf = inputs.Vinf

    # - Fill out wake strengths - #
    wake_vortex_strengths = dt.fill_out_wake_strengths(
        gamw, inputs.rotor_indices, inputs.num_wake_x_panels
    )

    # - Calculate body vortex strengths - #
    dt.calculate_body_vortex_strengths!(
        gamb,
        inputs.A_bb,
        inputs.b_bf,
        inputs.kutta_idxs,
        inputs.A_bw,
        wake_vortex_strengths,
        inputs.A_br,
        sigr,
    )

    # - Get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = dt.calculate_induced_velocities_on_rotors(
        blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        wake_vortex_strengths,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ inputs.Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* rpc

    # meridional component
    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # Get the inflow magnitude at the rotor as the combination of all the components
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    dt.calculate_gamma_sigma!(
        Gamr, sigr, inputs.blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor
    )

    # - Calculate net circulation and enthalpy jumps - #
    # TODO: check that your get property override works here for inputting an array of number of blades and rotation rates
    Gamma_tilde = dt.calculate_net_circulation(Gamr, blade_elements.B)
    H_tilde = dt.calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

    # - update wake strengths - #
    dt.calculate_wake_vortex_strengths!(
        gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
    )

    ## -- check body surface velocity initialiation -- ##
    # gamd = gamb[1:length(dp)]
    # gamh = gamb[(length(dp) + 1):end]
    gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
    # gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

    plot!(pb, dp, gamd; xlabel="x", ylabel="cp", label="iter$i duct surface pressure")
    # plot!(pb, hp, gamh ; label="iter hub surface pressure")

    ## -- check rotor circulation and source initial strengths -- ##
    plot!(
        pG, Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="iter$i"
    )

    plot!(
        ps, sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="iter$i"
    )

    plot!(
        pw,
        gamw,
        inputs.rotor_panel_edges;
        xlabel=L"\gamma_\theta",
        ylabel="r",
        label="iter$i",
    )
end

strengths = dt.analyze_propulsor(
    duct_coordinates,
    # hub_coordinates,
    nothing,
    paneling_constants,
    rotor_parameters,
    freestream;
)
gamb, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

plot!(pb, dp, gamd; xlabel="x", ylabel="cp", label="converged duct surface pressure")
# plot!(pb, hp, gamh ; label="iter hub surface pressure")

## -- check rotor circulation and source initial strengths -- ##
plot!(pG, Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="converged")
plot!(pG, ccbouts.circ, inputs.rotor_panel_centers; label="CCBlade, rotor only")

plot!(ps, sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="converged")

plot!(
    pw,
    gamw,
    inputs.rotor_panel_edges;
    xlabel=L"\gamma_\theta",
    ylabel="r",
    label="converged",
)

savefig(pb, "examples/body-init-sanity-check.pdf")
savefig(pG, "examples/rotorcirculation-init-sanity-check.pdf")
savefig(ps, "examples/rotorsources-init-sanity-check.pdf")
savefig(pw, "examples/wakegamma-init-sanity-check.pdf")

