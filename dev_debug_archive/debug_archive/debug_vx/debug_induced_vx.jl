#---------------------------------#
#             Includes            #
#---------------------------------#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

include(project_dir * "/plots_default.jl")
pyplot()

#---------------------------------#
#            Constants            #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

# Blade Hub radius, in meters
Rhub = 0.10 * Rtip

# number of blades
B = 2

# rotor x position
xrotor = 0.25

# duct distance from tip
tip_gap = 0.0

# Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
propgeom = [
    # 0.15 0.130 32.76
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
#
twists = propgeom[:, 3] * pi / 180

# use a basic airfoil model data
airfoil_file = project_dir * "/test/data/xrotor_af_test.dat"
airfoils = fill(ccb.AlphaAF(airfoil_file), length(rnondim))

#---------------------------------#
#      Define BODY Coordinates    #
#---------------------------------#
scale = 0.5

# - Duct Coordinates - #
duct_file = project_dir * "/test/data/naca_662-015.jl"
include(duct_file)
duct_coordinates = [x_duct r_duct] .* scale

# - Hub Coordinates - #
# include("../test/data/bodyofrevolutioncoords.jl")
# hub_coordinates = [x_hub[1:(end - 1)] ./ 2.0 r_hub[1:(end - 1)] * Rhub / maximum(r_hub)]
hub_coordinates = nothing

#---------------------------------#
#         Paneling Options        #
#---------------------------------#

nwake_sheets = 10
npanels_inlet = 10
discscale = 1.0

# non-dimensional wake length
wake_length = 1.0

ductle = minimum(duct_coordinates[:, 1])
ductte = maximum(duct_coordinates[:, 1])
ductchord = maximum(duct_coordinates[:, 1]) - minimum(duct_coordinates[:, 1])
outletinletratio = (ductte - xrotor) / (xrotor - ductle)

nhub_inlet = round(Int, npanels_inlet * discscale)

nduct_inlet = round(Int, npanels_inlet * discscale)

nduct_outlet = round(Int, nduct_inlet * outletinletratio)

nwake = round(Int, (nduct_inlet + nduct_outlet) * wake_length)

npanels = [nduct_outlet, nwake]

nducttot = (nduct_inlet + nduct_outlet) * 2
npersheet = nduct_outlet + nwake

#Vinf
Vinf = 5.0

# rotor rotation rate in rad/s
Omega = 5400 * pi / 30  # convert from RPM to rad/s

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
freestream = (; Vinf)

#---------------------------------#
#        Precompute Inputs        #
#---------------------------------#

# initialize various inputs used in analysis
inputs = dt.precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream;
    debug=true,
)

######################################################################
#                                                                    #
#           check isolated rotor/wake-induced velocities             #
#                                                                    #
######################################################################

# solve isolated rotor/wake system

roinputs, roparams = dt.initialize_rotor_states(
    rotor_parameters, freestream; wake_length=wake_length
)

rostates = dt.solve_rotor_only(roinputs, roparams)

roGamr, rogamw, rosigr = dt.extract_rotor_states(roinputs, roparams)

# fill out wake strengths
rowakegamma = dt.fill_out_wake_strengths(rogamw, [1], inputs.num_wake_x_panels)

#----------------------------------------------#
# look at wake-induced velocity on rotor plane #
#----------------------------------------------#

vx = similar(roGamr) .= 0 # axial induced velocity
vr = similar(roGamr) .= 0 # radial induced velocity
# add wake induced velocities
for jwake in 1:nwake_sheets
    vx .+= inputs.vx_rw[1, jwake] * rowakegamma[jwake, :]
    vr .+= inputs.vr_rw[1, jwake] * rowakegamma[jwake, :]
end

plot(; xlabel="wake-induced velocity", ylabel="r")
plot!(vx, inputs.rotor_source_panels[1].panel_center[:, 2]; label=L"v_x")
plot!(vr, inputs.rotor_source_panels[1].panel_center[:, 2]; label=L"v_r")
savefig(
    project_dir * "/dev_debug_archive/debug_vx/isolated-wake-induced-velocities-on-rotor-plane.pdf"
)

plot(; xlabel="wake strengths", ylabel="r")
plot!(rogamw, inputs.rotor_panel_edges; label="")
savefig(project_dir * "/dev_debug_archive/debug_vx/isolated-wake-strengths.pdf")

#---------------------------------------#
# look at wake-induced velocity on body #
#---------------------------------------#

pii = plot(; xlabel="x", ylabel="wake-induced inner surface velocity")
pio = plot(; xlabel="x", ylabel="wake-induced outer surface velocity")
# pci = plot(; xlabel="x", ylabel="wake-induced unit velocities on inner surface")
# pco = plot(; xlabel="x", ylabel="wake-induced unit velocities on outer surface")
#split upper and lower to make sure you can see what's going on
dp = inputs.body_panels[1].panel_center[:, 1]
_, leidx = findmin(dp)
#split into inner and outer surfaces
dpinner = dp[1:leidx]
dpouter = dp[(leidx + 1):end]

#set up colors
myr = range(165.0, 165.0, nwake_sheets) / 255.0
myg = range(209.0, 0.0, nwake_sheets) / 255.0
myb = range(255.0, 0.0, nwake_sheets) / 255.0

b_bw = similar(body_gammas) .= 0.0
for jwake in 1:nwake_sheets
    # get induced velocity in the x-direction
    b = inputs.A_bw[jwake] * rowakegamma[jwake, :]
    b_bw .+= b

    bwakeinner = b[1:leidx]
    bwakeouter = b[(leidx + 1):end]

    # awakeinner = inputs.A_bw[jwake][1:leidx]
    awakeouter = inputs.A_bw[jwake][(leidx + 1):end]

    col = RGB(myr[jwake], myg[jwake], myb[jwake])

    plot!(
        pii,
        dpinner,
        bwakeinner;
        color=col,
        label="wake sheet $(nwake_sheets -jwake) away from duct wall",
    )
    plot!(
        pio,
        dpouter,
        bwakeouter;
        color=col,
        label="wake sheet $(nwake_sheets -jwake) away from duct wall",
    )

    # plot!(
    #     pci,
    #     dpinner,
    #     awakeinner;
    #     color=col,
    #     label="wake sheet $(nwake_sheets -jwake) away from duct wall",
    # )
    # plot!(
    #     pco,
    #     dpouter,
    #     awakeouter;
    #     color=col,
    #     label="wake sheet $(nwake_sheets -jwake) away from duct wall",
    # )

    #attempt a countour plot of A_bw's
    pcont = contourf(
        [1:(inputs.num_wake_x_panels)],
        [1:length(body_gammas)],
        inputs.A_bw[jwake];
        color=:coolwarm,
        # clims=(minimum(inputs.A_bw[jwake]), maximum(inputs.A_bw[jwake])),
        title="wake sheet $(nwake_sheets-jwake) away from duct",
        xlabel="influencing wake panel from rotor to end",
        ylabel="affected body panels from inner TE clockwise",
    )
    savefig(
        pcont, project_dir * "/dev_debug_archive/debug_vx/wake$(jwake)-on-body-influence-matrix.pdf"
    )
end

savefig(
    pii,
    project_dir *
    "/dev_debug_archive/debug_vx/individual-wake-induced-velocities-on-inner-body-surface.pdf",
)
savefig(
    pio,
    project_dir *
    "/dev_debug_archive/debug_vx/individual-wake-induced-velocities-on-outer-body-surface.pdf",
)

# savefig(
#     pci,
#     project_dir *
#     "/dev_debug_archive/debug_vx/individual-UNIT-wake-induced-velocities-on-inner-body-surface.pdf",
# )
# savefig(
#     pco,
#     project_dir *
#     "/dev_debug_archive/debug_vx/individual-UNIT-wake-induced-velocities-on-outer-body-surface.pdf",
# )

bwakeinner = b_bw[1:leidx]
bwakeouter = b_bw[(leidx + 1):end]
plot(; xlabel="x", ylabel="wake-induced surface velocity")
plot!(dpinner, bwakeinner; label="inner surface, wake")
plot!(dpouter, bwakeouter; label="outer surface, wake")
savefig(
    project_dir * "/dev_debug_archive/debug_vx/isolated-wake-induced-velocities-on-body-surface.pdf"
)

######################################################################
#                                                                    #
#                 check isolated body-induced velocities             #
#                                                                    #
######################################################################

# solve isolated body system
body_gammas = dt.solve_body_system(inputs.A_bb, inputs.b_bf, inputs.kutta_idxs)

# sanity check on the body surface velocities
plot(inputs.body_panels[1].panel_center[:, 1], body_gammas; xlabel="x", ylabel=L"\gamma")
savefig(project_dir * "/dev_debug_archive/debug_vx/isolated-body-panel-strength-sanity-check.pdf")

#attempt a countour plot of A_bb
pcont = contourf(
    [1:length(body_gammas)],
    [1:length(body_gammas)],
    inputs.A_bb;
    color=:coolwarm,
    # clims=(minimum(inputs.A_bw[jwake]), maximum(inputs.A_bw[jwake])),
    xlabel="influencing body panels from inner TE clockwise",
    ylabel="affected body panels from inner TE clockwise",
)
savefig(pcont, project_dir * "/dev_debug_archive/debug_vx/body-on-body-influence-matrix.pdf")

#----------------------------------------------#
# look at body-induced velocity on rotor plane #
#----------------------------------------------#

vx = inputs.vx_rb[1] * body_gammas
vr = inputs.vr_rb[1] * body_gammas

plot(; xlabel="body-induced velocity", ylabel="r")
plot!(vx, inputs.rotor_source_panels[1].panel_center[:, 2]; label=L"v_x")
plot!(vr, inputs.rotor_source_panels[1].panel_center[:, 2]; label=L"v_r")
savefig(
    project_dir * "/dev_debug_archive/debug_vx/isolated-body-induced-velocities-on-rotor-plane.pdf"
)


#attempt a countour plot of A_rb
nrp, nbp = size(inputs.vx_rb[1])
pcont = contourf(
    [1:nbp],
    [1:nrp],
    inputs.vx_rb[1];
    color=:coolwarm,
    # clims=(minimum(inputs.A_bw[jwake]), maximum(inputs.A_bw[jwake])),
    xlabel="influencing body panels from inner TE clockwise",
    ylabel="affected rotor panels from hub to tip",
)
savefig(pcont, project_dir * "/dev_debug_archive/debug_vx/body-on-rotor-influence-matrix.pdf")

