#=
Verification and Validation of rotor/wake without duct and center body
=#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# TODO: change datapath to test directory or validation data directory
datapath = project_dir * "/validation/rotor_only/"
savepath = project_dir * "/validation/rotor_only/figures/"
dispath =
    project_dir *
    "/../../../../Writing/dissertation/src/ductsolvercontents/ductsolverfigures/"

#---------------------------------#
#             Includes            #
#---------------------------------#
using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade
include(datapath * "run_ccblade.jl")

include(project_dir * "/visualize/plots_default.jl")

#---------------------------------#
#             Geometry            #
#---------------------------------#

# Blade Tip Radius, in meters
Rtip = 10 / 2.0 * 0.0254  # inches to meters

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
r = propgeom[:, 1]
# Dimensionalize chords
chords = propgeom[:, 2] * Rtip
# convert twists to radians
twists = propgeom[:, 3] * pi / 180

# use a NACA 4412 airfoils
#=
Note here we are using the CCBlade fuCCBladectionality to define the airfoils data function.
In addition, we are using the airfoils data file available from the CCBlade repositCCBladery that has been extrapolated using the Viterna method as well as corrected for rotational effects as described in the CCBlade documentatioCCBlade.
=#
airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(r))

##### ----- User Options ----- #####
# number of blade elements to use in analysis
#=
Note: the solver will interpolate the rotor data using the given number of blade element inputs
=#
nwake_sheets = 18

# x position of rotor
xrotor = 0.0

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
#          Define Inputs          #
#---------------------------------#

# Rotor Parameters
rotor_parameters = [(;
    xrotor, nwake_sheets, r, chords, twists, airfoils, Rtip, Rhub, B, Omega
)]

# Wake Parameters
wake_length = 1.0
wake_x_refine = 0.05

# Freestream Parameters
Vinf = 5.0
freestream = (; rhoinf=rho, muinf=mu, asound, Vinf)

#---------------------------------#
#        Initialize Rotor         #
#---------------------------------#

nrotor = length(rotor_parameters)

##### ----- Panels ----- #####
# - Rotor Panel Edges - #
rotor_panel_edges = [
    range(rotor_parameters[i].Rhub, rotor_parameters[i].Rtip, nwake_sheets) for i in 1:nrotor
]
rotor_panel_edges = reduce(hcat, rotor_panel_edges)

# - Generate Panels - #
# note there will be nbe panels, but we require nbe+1 panel edges.
rotor_source_panels = [
    dt.generate_rotor_panels(rotor_parameters[i].xrotor, rotor_panel_edges[:, i]) for
    i in 1:nrotor
]

#get first rotor panel radial center points for convenience
rotor_panel_centers = rotor_source_panels[1].controlpoint[:, 2]
nbe = length(rotor_panel_centers)

##### ----- Blade Elements ----- #####
blade_elements = [
    dt.generate_blade_elements(
        rotor_parameters[i].B,
        rotor_parameters[i].Omega,
        rotor_parameters[i].xrotor,
        rotor_parameters[i].r,
        rotor_parameters[i].chords,
        rotor_parameters[i].twists,
        rotor_parameters[i].airfoils,
        rotor_parameters[i].Rtip,
        rotor_parameters[i].Rhub,
        rotor_source_panels[i].controlpoint[:, 2],
    ) for i in 1:nrotor
]

#---------------------------------#
#        Initialize Wake          #
#---------------------------------#

#=
Trailing wake sheets extend from between the blade elements, in other words, from the rotor panel edge locations.
=#

##### ----- General Geometry ----- #####

# - Choose blade element spacing - #
dr = wake_x_refine * rotor_parameters[1].Rtip

# - Rotor Diameter - #
D = 2.0 * rotor_parameters[1].Rtip

# - Define x grid spacing - #
#note that by using step rather than length, we aren't guarenteed to have the wake as long as we want, but it will be close.
xrange = range(
    rotor_parameters[1].xrotor, rotor_parameters[1].xrotor + wake_length * D; step=dr
)

num_wake_x_nodes = length(xrange)

# - Put together the initial grid point matrices - #
#=
For the grid relaxation (which adjusts the grid points such that the grids lie along streamlines) the function expects the grid points to be formatted such that there are x rows, and r columns.
In addition, we pass in 2 matrices, one for x and r such that all the x's are in one matrix, and all the r's in another.

note, however, that we are not yet using the wake relaxation because we do not have any duct bodies to deform the streamlines, so we will assume that the wake extends straight back for now.
NOTE THAT THIS MEANS WE AREN'T MODELING WAKE CONTRACTION, WHICH WILL LIKELY CAUSE DISCREPANCIES BETWEEN EXPERIMENT AND OUR MODEL, BUT THEY SHOULD BE REASONABLE.
=#
x_grid_points = repeat(xrange; outer=(1, nbe + 1))
r_grid_points = repeat(rotor_panel_edges'; outer=(length(xrange), 1))

##### ----- Panels ----- #####

#=
In order to make sure there aren't erronious panels from the end of one wake sheet to the beginning of the next, we have a panel object for each wake sheet, rather than for the entire wake.
Thus the wake_vortex_panels object is a vector of Panel structs, where each struct is comprised of a single row (sheet) of wake panels.
We will need to remember this in indexing such that we call wake_vortex_panels[row].property[column]
Note that there are nbe+1 trailing wake surfaces
=#
wake_vortex_panels = dt.generate_wake_panels(x_grid_points, r_grid_points)

wakeK = dt.get_wake_k(wake_vortex_panels)

# Go through the wake panels and determine the index of the aftmost rotor infront and the blade node from which the wake strength is defined.
rotorwakeid = ones(Int, wake_vortex_panels.totnode, 2)
for i in 1:(rotor_parameters[1].nwake_sheets)
    rotorwakeid[(1 + (i - 1) * num_wake_x_nodes):(i * num_wake_x_nodes), 1] .= i
end
for (i, wn) in enumerate(eachrow(wake_vortex_panels.node))
    rotorwakeid[i, 2] = findlast(x -> x <= wn[1], rotor_parameters.xrotor)
end

#---------------------------------#
# Calculate Influence Coefficients#
#---------------------------------#

# # TODO: update to linear panel influence functions
# function vortex_aic_boundary_on_field(
#     controlpoint, normal, tangent, node, nodemap, influence_length
# )

# unit induced velocity on rotor from wake sheets
v_rw = [
    -dt.induced_velocities_from_vortex_panels_on_points(
        rotor_source_panels[i].controlpoint,
        wake_vortex_panels.node,
        wake_vortex_panels.nodemap,
        wake_vortex_panels.influence_length,
        ones(wake_vortex_panels.totpanel, 2),
    ) for i in 1:length(rotor_source_panels)
]

# axial components
vx_rw = [v_rw[i][:, :, 1] for i in 1:length(rotor_source_panels)]

# radial components
vr_rw = [v_rw[i][:, :, 2] for i in 1:length(rotor_source_panels)]

##### ----- Rotor to Rotor ----- #####
# - rotor to rotor - #
v_rr = [
    dt.induced_velocities_from_source_panels_on_points(
        rotor_source_panels[i].controlpoint,
        rotor_source_panels[j].node,
        rotor_source_panels[j].nodemap,
        rotor_source_panels[j].influence_length,
        ones(rotor_source_panels[j].totpanel, 2),
    ) for i in 1:length(rotor_source_panels), j in 1:length(rotor_source_panels)
]

# axial components
vx_rr = [
    v_rr[i, j][:, :, 1] for i in 1:length(rotor_source_panels),
    j in 1:length(rotor_source_panels)
]

# radial components
vr_rr = [
    v_rr[i, j][:, :, 2] for i in 1:length(rotor_source_panels),
    j in 1:length(rotor_source_panels)
]

# - rotor to wake - #
v_wr = [
    dt.induced_velocities_from_source_panels_on_points(
        wake_vortex_panels.controlpoint,
        rotor_source_panels[j].node,
        rotor_source_panels[j].nodemap,
        rotor_source_panels[j].influence_length,
        ones(rotor_source_panels[j].totpanel, 2),
    ) for j in 1:length(rotor_source_panels)
]

# axial components
vx_wr = [v_wr[j][:, :, 1] for j in 1:length(rotor_source_panels)]

# radial components
vr_wr = [v_wr[j][:, :, 2] for j in 1:length(rotor_source_panels)]

# - wake to wake - #
v_ww = -dt.induced_velocities_from_vortex_panels_on_points(
    wake_vortex_panels.controlpoint,
    wake_vortex_panels.node,
    wake_vortex_panels.nodemap,
    wake_vortex_panels.influence_length,
    ones(wake_vortex_panels.totpanel, 2),
)

# axial components
vx_ww = v_ww[:, :, 1]

# radial components
vr_ww = v_ww[:, :, 2]

#---------------------------------#
#         Assemble inputs         #
#---------------------------------#
inputs = (;
    converged=[false],
    Vinf=freestream.Vinf,
    freestream,
    # nxwake=length(xrange) - 1,
    num_wake_x_nodes,
    wakeK,
    rotorwakeid,
    vx_rw=vx_rw,
    vr_rw=vr_rw,
    vx_rr=vx_rr,
    vr_rr=vr_rr,
    vx_wr=vx_wr,
    vr_wr=vr_wr,
    vx_ww=vx_ww,
    vr_ww=vr_ww,
    blade_elements,
    num_rotors=1,
    rotor_panel_edges=rotor_panel_edges,
    rotor_panel_centers=rotor_panel_centers,
    wake_vortex_panels,
    rotor_source_panels,
    xrange,
    x_grid_points,
    r_grid_points,
    nrotor_nodes=sum(rotor_source_panels.totnode),
    nwake_nodes=wake_vortex_panels.totnode,
    ductwakeinterfaceid=nothing,
    hubwakeinterfaceid=nothing,
)

#---------------------------------#
#         Set Up for Solves       #
#---------------------------------#

# TODO: need to update the various rotor/wake aero functions to account for linear panels.
#=
will need to update how things are averaged, etc.
for example, wake strengths at nodes are based on mean velocity of control points on either meridional side of the node, with no averaging at first and last wake nodes.
similarly, rotor source strengths are based on cd values at blade elements (panel control points)
states will have an additional source and wake value for each rotor and wake sheet, respectively.
=#

# - Initialize with freestream only - #
Wtheta = -inputs.rotor_panel_centers .* inputs.blade_elements.Omega'
# use freestream magnitude as meridional velocity at each blade section
Wm = similar(Wtheta) .= inputs.Vinf
# magnitude is simply freestream and rotation
W = sqrt.(Wtheta .^ 2 .+ Wm .^ 2)

# initialize circulation and source panel strengths
Gamr, sigr = dt.calculate_gamma_sigma(inputs.blade_elements, Wm, Wtheta, W, inputs.freestream)

nwakenode = inputs.wake_vortex_panels.totnode
gamw = zeros(nwakenode)
dt.calculate_wake_vortex_strengths!(
    gamw, Gamr, inputs.Vinf * ones(length(gamw)), inputs; debug=false
)

# get values needed for backing out freestream velocity from advance ratio
n = Omega / (2 * pi) #get revolutions per second
D = 2 * Rtip #rotor tip diameter

J = collect(range(0.1, 0.6; step=0.025))  # advance ratio
nJ = length(J)

rbe = inputs.blade_elements[1].rbe

# initialize outputs
eff = zeros(nJ)
CT = zeros(nJ)
CQ = zeros(nJ)
effccb = zeros(nJ)
CTccb = zeros(nJ)
CQccb = zeros(nJ)

# Loop through advance ratios
# for i in 1:nJ
for i in 1:1
    println()
    println("Running for J = $(J[i])")
    println()

    # calculate freestream velocity for given advance ratio
    Vinf_sweep = J[i] * D * n

    @time begin
        # run solver
        # remember to update freestream velocity in parameters
        states = dt.solve_rotor_only([Gamr; gamw; sigr], (; inputs..., Vinf=Vinf_sweep))
        Gamrconv, gamwconv, sigrconv = dt.extract_rotor_states(states, inputs)

        # - Post Process - #
        dtout = dt.states_to_outputs_rotor_only(states, (; inputs..., Vinf=Vinf_sweep))

        aero = dt.get_rotor_loads(
            dtout.W,
            dtout.phi,
            dtout.cl,
            dtout.cd,
            inputs.blade_elements[1],
            (; freestream..., Vinf=Vinf_sweep),
        )
    end

    # pG = plot(
    #     Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="initial"
    # )
    # plot!(pG, Gamrconv, inputs.rotor_panel_centers; label="converged")
    # # savefig(pG, savepath*"Circulation_J$(J[i]).pdf")

    CT[i] = aero.CT
    CQ[i] = aero.CQ
    eff[i] = aero.eff

    # Run CCBlade:
    ccbouts = run_ccblade(Vinf_sweep; airfoil="test/data/naca4412.dat")
    effccb[i] = ccbouts.eff
    CTccb[i] = ccbouts.CT
    CQccb[i] = ccbouts.CQ
    out = ccbouts.out

    ##### --- PLOTS --- ###
    ###Uncomment to see all the details
    #println("Plotting...")
    #println()

    # plot!(pG, ccbouts.circ, ccbouts.r; label="BEMT")
    #savefig(savepath*"circulation_J$(J[i]).pdf")
end

    # # double check geometry
    # plot(
    #     dtout.r,
    #     dtout.chord;
    #     xlabel="r",
    #     ylabel="chords",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(r * Rtip, chords; label="BEMT")
    # # savefig(savepath*"chord_J$(J[i]).pdf")

    # plot(
    #     dtout.r,
    #     dtout.twist * 180 / pi;
    #     ylabel="twists (deg)",
    #     xlabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(r * Rtip, twists * 180 / pi; label="BEMT")
    # savefig(savepath*"twist_J$(J[i]).pdf")

    # # axial induced velocity
    # plot(
    #     dtout.vx_rotor,
    #     rbe / Rtip;
    #     xlabel=L"v_x",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.u, r; label="BEMT")
    # savefig(savepath*"vx_J$(J[i]).pdf")

    # # tangential induced velocity
    # plot(
    #     dtout.vtheta_rotor,
    #     rbe / Rtip;
    #     xlabel=L"v_\theta",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.v, r; label="BEMT")
    # savefig(savepath*"vtheta_J$(J[i]).pdf")

    # # tangential total velocity
    # plot(
    #     dtout.Wθ,
    #     rbe / Rtip;
    #     xlabel=L"W_\theta",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.v .- Omega * r * Rtip, r; label="BEMT")
    # savefig(savepath*"Wtheta_J$(J[i]).pdf")

    # # meridional total velocity
    # plot(
    #     dtout.Wm,
    #     rbe / Rtip;
    #     xlabel=L"W_m",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.u .+ Vinf_sweep, r; label="BEMT")
    # savefig(savepath*"Wm_J$(J[i]).pdf")

    # # inflow angle
    # plot(
    #     dtout.phi * 180 / pi,
    #     rbe / Rtip;
    #     xlabel=L"\phi~(deg)",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.phi * 180 / pi, r; label="BEMT")
    # savefig(savepath*"phi_J$(J[i]).pdf")

    # # angle of attack
    # plot(
    #     dtout.alpha * 180 / pi,
    #     rbe / Rtip;
    #     xlabel=L"\alpha~(deg)",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.alpha * 180 / pi, r; label="BEMT")
    # savefig(savepath*"alpha_J$(J[i]).pdf")

    # # inflow magnitude
    # plot(
    #     dtout.W, rbe / Rtip; xlabel=L"W", ylabel="r", label="DuctAPE", title="J = $(J[i])"
    # )
    # plot!(out.W, r; label="BEMT")
    # savefig(savepath*"W_J$(J[i]).pdf")

    # # Lift
    # plot(
    #     dtout.cl,
    #     rbe / Rtip;
    #     xlabel=L"c_\ell",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(dtout.clin, rbe / Rtip; label="inner", linestyle=:dash, color=1)
    # plot!(dtout.clout, rbe / Rtip; label="outer", linestyle=:dot, color=1)
    # plot!(out.cl, r; label="BEMT")
    # savefig(savepath*"cl_J$(J[i]).pdf")

    # Drag
    # plot(
    #     dtout.cd,
    #     rbe / Rtip;
    #     xlabel=L"c_d",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(dtout.cdin, rbe / Rtip; label="inner", linestyle=:dash, color=1)
    # plot!(dtout.cdout, rbe / Rtip; label="outer", linestyle=:dot, color=1)
    # plot!(out.cd, r; label="BEMT")
    # savefig(savepath*"cd_J$(J[i]).pdf")

    # # normal coeff
    # plot(
    #     aero.cn,
    #     rbe / Rtip;
    #     xlabel=L"c_n",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.cn, r; label="BEMT")
    # savefig(savepath*"cn_J$(J[i]).pdf")

    # # tangential coeff
    # plot(
    #     aero.ct,
    #     rbe / Rtip;
    #     xlabel=L"c_t",
    #     ylabel="r",
    #     label="DuctAPE",
    #     title="J = $(J[i])",
    # )
    # plot!(out.ct, r; label="BEMT")
    # savefig(savepath*"ct_J$(J[i]).pdf")

    # # Distributed Loads
    # plot(
    #     rbe / Rtip,
    #     aero.Np;
    #     xlabel="r/Rtip",
    #     label="Normal load per unit length",
    #     title="J = $(J[i])",
    # )
    # plot!(r, out.Np; label="BEMT Np")
    # plot!(rbe / Rtip, aero.Tp; label="Tangential load per unit length")
    # plot!(r, out.Tp; label="BEMT Tp")
    # savefig(savepath*"distributed_loads_J$(J[i]).pdf")
# end

#---------------------------------#
#        Experimental Data        #
#---------------------------------#
exp = [
    0.113 0.0912 0.0381 0.271
    0.145 0.0890 0.0386 0.335
    0.174 0.0864 0.0389 0.387
    0.200 0.0834 0.0389 0.429
    0.233 0.0786 0.0387 0.474
    0.260 0.0734 0.0378 0.505
    0.291 0.0662 0.0360 0.536
    0.316 0.0612 0.0347 0.557
    0.346 0.0543 0.0323 0.580
    0.375 0.0489 0.0305 0.603
    0.401 0.0451 0.0291 0.620
    0.432 0.0401 0.0272 0.635
    0.466 0.0345 0.0250 0.644
    0.493 0.0297 0.0229 0.640
    0.519 0.0254 0.0210 0.630
    0.548 0.0204 0.0188 0.595
    0.581 0.0145 0.0162 0.520
]
# extract advance ratios
Jexp = exp[:, 1]
# extract thrust coefficients
CTexp = exp[:, 2]
# extract power coefficients
CPexp = exp[:, 3]
# extract efficiencies
etaexp = exp[:, 4]

#---------------------------------#
#              Plots              #
#---------------------------------#

plot(J, CT; xlabel=L"\mathrm{Advance~Ratio~}(J)", label=L"DuctAPE", color=1)
plot!(J, CQ * 2 * pi; label="", color=2)
plot!(J, CTccb; label=L"BEMT", color=1, linestyle=:dash)
plot!(J, CQccb * 2 * pi; label="", color=2, linestyle=:dash)
plot!(
    Jexp, CTexp; seriestype=:scatter, markershape=:utriangle, label="Experimental", color=1
)
plot!(Jexp, CPexp; seriestype=:scatter, markershape=:utriangle, label="", color=2)
annotate!(0.4, 0.07, text(L"C_T", 11, :right, myblue))
annotate!(0.4, 0.0225, text(L"C_P", 11, :right, myred))
savefig(savepath * "rotor-only-thrust-and-power-validation.pdf")

plot(
    J,
    eff;
    xlabel=L"\mathrm{Advance~Ratio~}(J)",
    ylabel=L"\mathrm{Efficiency~}(\eta)",
    color=1,
    label="DuctAPE",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(J, effccb; label="BEMT", linestyle=:dash, color=1)
plot!(
    Jexp, etaexp; seriestype=:scatter, markershape=:utriangle, color=1, label="experimental"
)
savefig(savepath * "rotor-only-efficiency-validation.pdf")

# TODO:
#=
Want to show the rotor/wake induced velocities, both meridional and swirl at points swept from "far" upstream, across the rotor plane, to "far" downstream.
Probably pick the representative radial position that Farokhi uses in simplified models as a good point at which to sample.  Will want to go 1 diamter upstream to 1 diameter downstream, making sure the wake is sufficiently long downstream.
prop is probably different, so will likely want to do this in a different file.
actually is just comparision with near and far field from BEMT, still probably want to do in a separate file though, just because this one is super long.
=#
