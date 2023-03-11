#=

Iterated wake using the addition of rotor sources.

Also moving to expresssion for wake panel vortex strength using Vm_avg, which requires velocity sampling throughout the wake, and therefore the precomputation of the induced velocity terms.

=#

#---------------------------------#
#             Includes            #
#---------------------------------#

# - Include Plots and my default settings - #
include("../../../plots_default.jl")

# - use CCBlade for comparison - #
using CCBlade

# - Rename DuctTAPE for convenience - #
using DuctTAPE
const dt = DuctTAPE

# - Rename FLOWFoil for convenience - #
using FLOWFoil
const ff = FLOWFoil

######################################################################
#                                                                    #
#                      SET UP ALL THE GEOMETRY                       #
#                                                                    #
######################################################################

#---------------------------------#
#              Rotor              #
#---------------------------------#
#=
Rotor Blade Sections DO NOT extend to hub and tip
Rotor panels start and stop at the midpoints between blade elements, with the blade tip and hub being the outermost panel edges.
TODO: actually need to think about this more.  if panels start and stop and hub and tip, then outer most panels are weird and don't line up with input data in a way that guarentees that interpolation can work.
=#

##### ----- General Geometry ----- #####
#=
For this test case, we use the geometry presented in the documentation of CCBlade, with some slight modifications since we can't handle blade element definitions at the hub and tip as of yet.
=#

## -- Tip and Hub Radii -- ##
Rtip = 10 / 2.0 * 0.0254  # inches to meters
Rhub = 0.10 * Rtip

# - r, chord, twist of CCBlade example, minus tip section data - #
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
]

# - Radial Stations - #
r = propgeom[:, 1] * Rtip

# - Rename the number of blade elements for convenience - #
nbe = length(r)

# - Chord Distribution - #
chord = propgeom[:, 2] * Rtip

# - Twist Distribution - #
theta = propgeom[:, 3] * pi / 180

# - Number of Blades - #
B = 2

# - Airfoils - #
af = AlphaAF("test/data/naca4412.dat")

##### ----- Panels ----- #####
#=
For now, we just set the x-location of the rotor to zero since we don't have anywhere better to place it.
Note that the panels on the rotor are the source panels modeling the blade element profile drag, centered at the blade element locations.
=#
# - Rotor Panel Edges - #
# start at hub, and edges are midpoints between blade element locations
rotor_panel_edges = [Rhub; (r[2:end] + r[1:(end - 1)]) / 2.0; Rtip]

# For the rotor panels, we want to use Constant Source Panels.
# The boundary condition and body of revolution boolean don't matter, should probably TODO: add a constructor that doesn't require those when not needed.
rotor_method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])
# - Generate Panels - #
# note there will be nbe panels, but we require nbe+1 panel edges.
rotor_panels = ff.generate_panels(rotor_method, [zeros(nbe + 1) rotor_panel_edges])

##### ----- SANITY PLOTS ----- #####

# - plot the original radial stations - #
plot(
    zeros(nbe),
    r;
    seriestype=:scatter,
    markershape=:square,
    markersize=3,
    label="Input Blade Element Positions",
)

# - plot the panel edges - #
plot!(
    zeros(nbe + 1),
    rotor_panel_edges;
    markershape=:rect,
    markersize=1,
    linewidth=0.5,
    label="Rotor Panels",
)

# - Plot panel centers - #
plot!(
    rotor_panels.panel_center[:, 1],
    rotor_panels.panel_center[:, 2];
    seriestype=:scatter,
    markershape=:circle,
    markersize=1.5,
    label="Panel Centers",
)

# - Plot Rtip and Rhub - #
plot!([-0.1; 0.1], [Rtip; Rtip]; linestyle=:dash, color=:black, label="Rtip and Rhub")
plot!([-0.1; 0.1], [Rhub; Rhub]; linestyle=:dash, color=:black, label="")

# - save figure - #
savefig("test/manual_tests/rotor_wake_tests/sanity_check_rotor_geometry.pdf")

#---------------------------------#
#              Wake               #
#---------------------------------#
#=
Trailing wake sheets extend from between the blade elements, in other words, from the rotor panel edge locations.
In general, attempt to keep wake "grids" somewhat square by setting the initial x-spacing to be similar to that of the rotor panel lengths (assuming the rotor panels are spaced linearly across the rotor).
We will also choose the wake to be 2 times the rotor diameter for now. (when adding the duct, we will take the larger of that and extending 1 duct chord behind the duct trailing edge, whichever is longer).
=#

##### ----- General Geometry ----- #####

# - Get blade element spacing - #
#note that the spacing is constant above, if not, we may have selected the average spacing.
dr = r[end] - r[end - 1]

# - Rotor Diameter - #
D = 2.0 * Rtip

# - Define x grid spacing - #
#note that by using step rather than length, we aren't guarenteed to have the wake as long as we want, but it will be close.
xrange = range(0.0, 2.0 * D; step=dr)

# - Put together the initial grid point matrices - #
#=
For the grid relaxation (which adjusts the grid points such that the grids lie along streamlines) the function expects the grid points to be formatted such that there are x rows, and r columns.
In addition, we pass in 2 matrices, one for x and r such that all the x's are in one matrix, and all the r's in another.

note, however, that we are not yet using the wake relaxation because we do not have any duct bodies to deform the streamlines, so we will assume that the wake extends straight back for now.
NOTE THAT THIS MEANS WE AREN'T MODELING WAKE CONTRACTION, WHICH WILL LIKELY CAUSE DISCREPANCIES BETWEEN CCBLADE AND OUR MODEL, BUT THEY SHOULD BE REASONABLE.
=#
x_grid_points = repeat(xrange; outer=(1, nbe + 1))
r_grid_points = repeat(rotor_panel_edges'; outer=(length(xrange), 1))

##### ----- Panels ----- #####

# For the wake panels, we want to use Constant Vortex Panels.
# The boundary condition and body of revolution boolean don't matter, should probably TODO: add a constructor that doesn't require those when not needed.
wake_method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])

#=
In order to make sure there aren't erronious panels from the end of one wake sheet to the beginning of the next, we have a panel object for each wake sheet, rather than for the entire wake.
Thus the wake_panels object is a vector of Panel structs, where each struct is comprised of a single row (sheet) of wake panels.
We will need to remember this in indexing such that we call wake_panels[row].property[column]
Note that there are nbe+1 trailing wake surfaces
=#
wake_panels = dt.generate_wake_panels(
    x_grid_points, r_grid_points, nbe + 1; method=wake_method
)

##### ----- SANITY PLOTS ----- #####

# - plot the original radial stations - #
for i in 1:length(x_grid_points[1, :])
    if i == 1
        plot(
            x_grid_points[:, i],
            r_grid_points[:, i];
            markershape=:rect,
            markersize=0.5,
            linewidth=0.25,
            color=mycolors[2],
            aspectratio=:equal,
            label="Input Wake Grid (panel edges)",
        )
    else
        plot!(
            x_grid_points[:, i],
            r_grid_points[:, i];
            markershape=:rect,
            markersize=0.5,
            linewidth=0.25,
            color=mycolors[2],
            aspectratio=:equal,
            label="",
        )
    end
end

# - Plot panel centers - #
for i in 1:length(wake_panels)
    if i == 1
        plot!(
            wake_panels[i].panel_center[:, 1],
            wake_panels[i].panel_center[:, 2];
            seriestype=:scatter,
            markershape=:circle,
            markersize=0.5,
            color=mycolors[3],
            label="Panel Centers",
        )
    else
        plot!(
            wake_panels[i].panel_center[:, 1],
            wake_panels[i].panel_center[:, 2];
            seriestype=:scatter,
            markershape=:circle,
            markersize=0.5,
            color=mycolors[3],
            label="",
        )
    end
end

# - plot where 2D point is - #
plot!(
    [2.0 * D; 2.0 * D], [Rhub, Rtip]; linestyle=:dash, color=:black, label="2D Downstream"
)

# - save figure - #
savefig("test/manual_tests/rotor_wake_tests/sanity_check_wake_geometry.pdf")

######################################################################
#                                                                    #
#                  PRECOMPUTE UNIT INDUCED VELOCITIES                #
#                                                                    #
######################################################################

#---------------------------------#
#             MESHES              #
#---------------------------------#

##### ----- Rotor to Wake ----- #####
rotor_to_wake_mesh = [
    dt.generate_one_way_mesh(rotor_panels, wake_panels[i]) for i in 1:length(wake_panels)
]

##### ----- Rotor to Rotor ----- #####
rotor_to_rotor_mesh = dt.generate_one_way_mesh(rotor_panels, rotor_panels)

##### ----- Wake to Rotor ----- #####
wake_to_rotor_mesh = [
    dt.generate_one_way_mesh(wake_panels[i], rotor_panels) for i in 1:length(wake_panels)
]

##### ----- Wake to Wake ----- #####
wake_to_wake_mesh = [
    dt.generate_one_way_mesh(wake_panels[i], wake_panels[j]) for i in 1:length(wake_panels),
    j in 1:length(wake_panels)
]

#---------------------------------#
#           Velocities            #
#---------------------------------#

##### ----- Rotor to Wake ----- #####
A_rotor_to_wake = [
    dt.assemble_induced_velocity_matrices(
        rotor_to_wake_mesh[i], rotor_panels, wake_panels[i]
    ) for i in 1:length(wake_panels)
]

# - Axial - #
vxd_rotor_to_wake = [A_rotor_to_wake[i][1] for i in 1:length(wake_panels)]

# - Radial - #
vrd_rotor_to_wake = [A_rotor_to_wake[i][2] for i in 1:length(wake_panels)]

##### ----- Rotor to Rotor ----- #####
A_rotor_to_rotor = dt.assemble_induced_velocity_matrices(
    rotor_to_rotor_mesh, rotor_panels, rotor_panels
)

# - Axial - #
vxd_rotor_to_rotor = A_rotor_to_rotor[1]

# - Radial - #
vrd_rotor_to_rotor = A_rotor_to_rotor[2]

##### ----- Wake to Rotor ----- #####
A_wake_to_rotor = [
    dt.assemble_induced_velocity_matrices(
        wake_to_rotor_mesh[i], wake_panels[i], rotor_panels
    ) for i in 1:length(wake_panels)
]

# - Axial - #
vxd_wake_to_rotor = [A_wake_to_rotor[i][1] for i in 1:length(wake_panels)]

# - Radial - #
vrd_wake_to_rotor = [A_wake_to_rotor[i][2] for i in 1:length(wake_panels)]

##### ----- Wake to Wake ----- #####
A_wake_to_wake = [
    dt.assemble_induced_velocity_matrices(
        wake_to_wake_mesh[i], wake_panels[i], wake_panels[j]
    ) for i in 1:length(wake_panels), j in 1:length(wake_panels)
]

# - Axial - #
vxd_wake_to_wake = [
    A_wake_to_wake[i, j][1] for i in 1:length(wake_panels), j in 1:length(wake_panels)
]

# - Radial - #
vrd_wake_to_wake = [
    A_wake_to_wake[i, j][2] for i in 1:length(wake_panels), j in 1:length(wake_panels)
]

######################################################################
#                                                                    #
#                        OPERATING CONDITIONS                        #
#                                                                    #
######################################################################

# - Freestream Velocity - #
Vinf = 5.0

# - Rotor Rotation Rate in radians/second - #
Omega = 5400 * pi / 30  # convert to rad/s

# - Freestream Properties - #
rho = 1.225
mu = 1.81e-5
asound = 343.0

######################################################################
#                                                                    #
#                        RUN CCBLADE EXAMPLE                         #
#                                                                    #
######################################################################
rotor = Rotor(Rhub, Rtip, B; tip=nothing)
sections = Section.(r, chord, theta, Ref(af))
op = simple_op.(Vinf, Omega, r, rho)
out = CCBlade.solve.(Ref(rotor), sections, op)

## -- Get Circulation from outputs -- ##
function get_gamma_sigma(ccbouts, chord, r)
    cl = ccbouts.cl
    cd = ccbouts.cd
    W = ccbouts.W

    Gamma = 0.5 .* W .* cl .* chord
    Sigma = B .* W .* cd .* chord ./ (4.0 * pi * r)

    return Gamma, Sigma
end

ccbGamma, ccbSigma = get_gamma_sigma(out, chord, r)

######################################################################
#                                                                    #
#                        DuctTAPE FUNCTIONS                          #
#                                                                    #
######################################################################

#---------------------------------#
#       Vmavg_wake Functions      #
#---------------------------------#
#=
Functions for obtaining the average meridional velocity on the wake panel centers

We need these velocities to calculate the wake panel strengths (along with contributions from the circulation and enthalpy jumps from the rotor)

total meridional velocity will include the freestream as well as the induced axial and radial velocities from the other wake panels and the rotor source panels.
=#

"""
Use the induced velocities of the wake and rotor along with the freestream to calculate the meridional velocity on each of the wake panel centers.

use this version when you need to initialize the wake velocities

TODO: may need to just initialize to Vinf, rather than getting everything. need to think about the order of operations here.
"""
function get_wake_vmavg(
    Vinf,
    wake_vortex_strengths,
    vxd_wake_to_wake,
    vrd_wake_to_wake,
    rotor_source_strengths,
    vxd_rotor_to_wake,
    vrd_rotor_to_wake,
)

    # - Initialize output with r rows and x columns starting with freestream influence - #
    #(TODO: check that this will work with what you need)
    Vmavg_wake =
        Vinf .*
        ones(length(wake_vortex_strengths[:, 1]), length(wake_vortex_strengths[1, :]))

    get_wake_vmavg!(
        Vmavg_wake,
        wake_vortex_strengths,
        vxd_wake_to_wake,
        vrd_wake_to_wake,
        rotor_source_strengths,
        vxd_rotor_to_wake,
        vrd_rotor_to_wake,
    )

    return Vmavg_wake
end

"""
use this version after Vmavg_wake has been initialized
"""
function get_wake_vmavg!(
    Vmavg_wake,
    wake_vortex_strengths,
    vxd_wake_to_wake,
    vrd_wake_to_wake,
    rotor_source_strengths,
    vxd_rotor_to_wake,
    vrd_rotor_to_wake,
)

    # - add wake on wake contributions to Vmavg on wake panels
    add_wake_on_wake_inducement!(
        Vmavg_wake, wake_vortex_strengths, vxd_wake_to_wake, vrd_wake_to_wake
    )

    # - add rotor on wake contributions to Vmavg on wake panels
    add_rotor_on_wake_inducement!(
        Vmavg_wake, rotor_source_strengths, vxd_rotor_to_wake, vrd_rotor_to_wake
    )

    return nothing
end

"""
"""
function add_wake_on_wake_inducement!(
    Vmavg_wake, wake_vortex_strengths, vxd_wake_to_wake, vrd_wake_to_wake
)

    ##### ----- Wake on Wake Influence ----- #####
    # - Loop through the affected wakes - #
    for i in 1:length(vxd_wake_to_wake[:, 1])
        # - Loop through the wakes doing the influencing - #
        for j in 1:length(vxd_wake_to_wake[1, :])
            # - compute the affect of the entire wake row on the entire wake row - #

            # axial direction
            vx_row_j_on_i = vxd_wake_to_wake[i, j] * wake_vortex_strengths[j, :]
            # radial direction
            vr_row_j_on_i = vrd_wake_to_wake[i, j] * wake_vortex_strengths[j, :]

            # - Put the axial and radial components together to get the meridional velocity contribution from the jth wake on the ith wake- #
            Vmavg_wake[i, :] .+= sqrt.(vx_row_j_on_i .^ 2 .+ vr_row_j_on_i .^ 2)
        end
    end

    return nothing
end

"""
"""
function add_rotor_on_wake_inducement!(
    Vmavg_wake, rotor_source_strengths, vxd_rotor_to_wake, vrd_rotor_to_wake
)

    ##### ----- Rotor on Wake Influence ----- #####
    # - Loop through the affected wakes - #
    for i in 1:length(vxd_rotor_to_wake)
        # - compute the affect of all the rotor panels on the entire wake row - #

        # axial direction
        vx_rotor_on_wake_i = vxd_rotor_to_wake[i] * rotor_source_strengths
        # radial direction
        vr_rotor_on_wake_i = vrd_rotor_to_wake[i] * rotor_source_strengths

        # - Put the axial and radial components together to get the meridional velocity contribution from the jth wake on the ith wake- #
        Vmavg_wake[i, :] .+= sqrt.(vx_rotor_on_wake_i .^ 2 .+ vr_rotor_on_wake_i .^ 2)
    end

    return nothing
end

#---------------------------------#
#     Rotor Inflow Functions      #
#---------------------------------#
#=
Functions for obtaining the meridional, tangential, and total inflow velocities to the rotor plane.

We use these velocities to calculate the blade element information such as angle of attack and local circulation on the blade.

total meridional velocity will include the freestream as well as the induced axial and radial velocities from the wake panels and the rotor source panels.
Tangential velocity will come from the rotor self-inducement directly without any sort of panel computations.
=#

"""
Note that Vm = Wm, but Wtheta includes the rotational component where Vtheta does not.
"""
function get_rotor_inflow(
    Vinf,
    Omega,
    radial_positions,
    BGamma,
    wake_vortex_strengths,
    vxd_wake_to_rotor,
    vrd_wake_to_rotor,
    rotor_source_strengths,
    vxd_rotor_to_rotor,
    vrd_rotor_to_rotor,
)

    # - Initialize Outputs - #
    Vm = Vinf .* ones(length(radial_positions))
    Wtheta = -Omega .* radial_positions

    # - add the wake induced velocities at the rotor plane to Vm - #
    add_wake_on_rotor_vm_inducement!(
        Vm, wake_vortex_strengths, vxd_wake_to_rotor, vrd_wake_to_rotor
    )

    # - and the rotor source panel induced velocities at the rotor plane to Vm - #
    add_rotor_on_rotor_vm_inducement!(
        Vm, rotor_source_strengths, vxd_rotor_to_rotor, vrd_rotor_to_rotor
    )

    # - add the rotor rotational self-induction directly - #
    add_rotor_self_induced_vtheta!(Wtheta, radial_positions, BGamma)

    # - Get the inflow magnitude - #
    Wmag = sqrt.(Vm .^ 2 .+ Wtheta .^ 2)

    # - return everything as a convenient named tuple (not necessary, may change later) - #
    return (Wm=Vm, Wtheta=Wtheta, Wmag=Wmag)
end

"""
TODO: will need to think about how to change things for multiple rotors
"""
function add_rotor_on_rotor_vm_inducement!(
    Vm, rotor_source_strengths, vxd_rotor_to_rotor, vrd_rotor_to_rotor
)

    ##### ----- Rotor on Rotor Influence ----- #####

    # axial direction
    vx_rotor_on_rotor = vxd_rotor_to_rotor * rotor_source_strengths
    # radial direction
    vr_rotor_on_rotor = vrd_rotor_to_rotor * rotor_source_strengths

    # - Put the axial and radial components together to get the meridional velocity contribution from the jth wake on the ith wake- #
    Vm .+= sqrt.(vx_rotor_on_rotor .^ 2 .+ vr_rotor_on_rotor .^ 2)

    return nothing
end

function add_wake_on_rotor_vm_inducement!(
    Vm, wake_vortex_strengths, vxd_wake_to_rotor, vrd_wake_to_rotor
)

    ##### ----- Wake on Rotor Influence ----- #####
    # - Loop through the wakes - #
    for i in 1:length(vxd_wake_to_rotor)

        # axial direction
        vx_wake_i_on_rotor = vxd_wake_to_rotor[i] * wake_vortex_strengths[i, :]
        # radial direction
        vr_wake_i_on_rotor = vrd_wake_to_rotor[i] * wake_vortex_strengths[i, :]


        # - Put the axial and radial components together to get the meridional velocity contribution from the jth wake on the ith wake- #
        Vm .+= sqrt.(vx_wake_i_on_rotor .^ 2 .+ vr_wake_i_on_rotor .^ 2)
    end

    return nothing
end

"""
TODO: eventually want to add upstream rotor contributions as well.
"""
function add_rotor_self_induced_vtheta!(Wtheta, radial_positions, BGamma)

    # the rotor adds half of its own circulation to the self-induced tangential velocity
    Wtheta .+= 1.0 ./ (2.0 .* pi .* r) .* (0.5 .* BGamma)

    return nothing
end

#---------------------------------#
# Calculate Wake Vortex Strenghts #
#---------------------------------#
#=
Here we use the expression for vortex strength using the average meridional velocity computed directly on the panel rather than marching the velocity jumps from the outer edge of the wake. This way it's easier to just add all the contributions together from the wake and source panels.
=#

"""
"""
function gamma_theta_from_vmavg(Vmavg, Gamma_tilde, H_tilde, r)

    # - Initialize Output - #
    gamma_theta = similar(Vmavg)

    # - Loop through the radial positions - #
    nw = length(Vmavg[:, 1])
    for i in 1:nw

        # - Loop through the axial positions - #
        for j in 1:length(Vmavg[1, :])

            # - if at the hub - #
            # for now, just set circulation etc. to zero outside of wake
            if i == 1
                gamma_theta[i, j] =
                    1.0 / (2.0 * Vmavg[i, j]) *
                    (-(1.0 / (2.0 * pi * r[i]))^2 * (Gamma_tilde[i]^2) + 2.0 * (H_tilde[i]))
                # - if at the tip - #
            elseif i == nw
                gamma_theta[i, j] =
                    1.0 / (2.0 * Vmavg[i, j]) * (
                        -(1.0 / (2.0 * pi * r[i]))^2 * (-Gamma_tilde[i - 1]^2) +
                        2.0 * (-H_tilde[i - 1])
                    )
                # - otherwise - #
            else
                gamma_theta[i, j] =
                    1.0 / (2.0 * Vmavg[i, j]) * (
                        -(1.0 / (2.0 * pi * r[i]))^2 *
                        (Gamma_tilde[i]^2 - Gamma_tilde[i - 1]^2) +
                        2.0 * (H_tilde[i] - H_tilde[i - 1])
                    )
            end
        end
    end

    return gamma_theta
end

######################################################################
#                                                                    #
#                           INITIALIZE                               #
#                                                                    #
######################################################################

# # - Define some basic blade element named tuples as place holders
# blade_elements = [(
#     omega=Omega, radial_positions=rpc, chords=chord
# )]

# # initialize induced velocities at rotors to zeros
# vm_init = zeros(nbe)
# vtheta_init = zeros(nbe)

# # calculate inflow from freestream and rotation
# inflow_init = dt.calculate_inflow_velocities(
#     blade_elements, Ref(Vinf), vm_init, vtheta_init
# )

# calculate angle of attack
# alpha_init = dt.calculate_angle_of_attack(theta, inflow_init.Wm, inflow_init.Wtheta)

# # - Look up lift and drag data - #
# cl_init = zeros(nbe)
# cd_init = zeros(nbe)
# for a in 1:nbe
#     cl_init[a], cd_init[a] = dt.search_polars(af, alpha_init[a])
# end

# - Calculate Circulation from lift and inflow_init - #
# Gamma = 0.5 .* inflow_init.Wmag .* chord .* cl_init
Gamma = deepcopy(ccbGamma)
rotor_source_strengths = deepcopy(ccbSigma)

# get net circulation and B*Gamma
BGamma_init, Gamma_tilde_init = dt.calculate_net_circulation(Gamma, B)

# - Calculate Enthalpy Jumps - #
H_tilde_init = dt.calculate_enthalpy_jumps(Gamma, Omega, B)

# - Calculate Vm in wake - #
Vmavg_wake = Vinf .* ones(length(wake_panels), length(wake_panels[1].panel_center[:, 1]))

# - Calculate gamma_thetas - #
wake_vortex_strengths = gamma_theta_from_vmavg(
    Vmavg_wake, Gamma_tilde_init, H_tilde_init, rotor_panel_edges
)

## -- Set up some plots for iteration -- ##
rpc = rotor_panels.panel_center[:, 2]

pcirc = plot(ccbGamma, r; xlabel=L"\Gamma = 0.5 W c c_\ell", ylabel="r", label="CCBlade")

psigma = plot(
    ccbSigma, r; xlabel=L"\sigma = \frac{B}{4\pi r}W c c_d", ylabel="r", label="CCBlade"
)

palpha = plot(
    out.alpha * 180.0 / pi, r; xlabel=L"\alpha (degrees)", ylabel="r", label="CCBlade"
)

pvm = plot(
    out.u,
    r;
    xlabel="total axial (meridional) velocity induced at rotor",
    ylabel="r",
    label="CCBlade",
)

pvt = plot(
    out.v, r; xlabel="tangential velocity induced at rotor", ylabel="r", label="CCBlade"
)

pw = plot(out.W, r; xlabel="inflow magnitude (Wmag)", ylabel="r", label="CCBlade")
pcl = plot(out.cl, r; xlabel="lift coefficient", ylabel="r", label="CCBlade")
pcd = plot(out.cd, r; xlabel="drag coefficient", ylabel="r", label="CCBlade")

pgammatheta = plot(
    wake_vortex_strengths[:, 41],
    rotor_panel_edges;
    xlabel="wake vortex strengths at midwake",
    ylabel="r",
    label="Init",
)

# pvmwake = plot(; xlabel="wake-induced meridional velocity at rotor plane", ylabel="r")

######################################################################
#                                                                    #
#                             ITERATE                                #
#                                                                    #
######################################################################
Gamma_temp = 99 * ones(length(Gamma))
iter = [0]
while abs(Gamma_temp[5] - Gamma[5]) > 1e-5
# for i in 1:5
    iter[1] += 1
    println("iter $(iter[1])")

    # println("Before Iteration\n")
    # println(" Vinf")
    # println(Vinf)
    # println(" Omega")
    # println(Omega)
    # println(" B*Gamma")
    # println(B * Gamma)
    # println(" rotor_panels.panel_center[:,2]")
    # println(rpc)
    # println(" wake_vortex_strengths")
    # display(wake_vortex_strengths)
    # # println(" vxd_wake_to_rotor")
    # # println(" vrd_wake_to_rotor")
    # println(" rotor_source_strengths")
    # println(rotor_source_strengths)
    # # println(" vxd_rotor_to_rotor")
    # # println(" vrd_rotor_to_rotor")

    inflow = get_rotor_inflow(
        Vinf,
        Omega,
        rpc,
        B * Gamma,
        wake_vortex_strengths,
        vxd_wake_to_rotor,
        vrd_wake_to_rotor,
        rotor_source_strengths,
        vxd_rotor_to_rotor,
        vrd_rotor_to_rotor,
    )

    # calculate angle of attack
    alpha = dt.calculate_angle_of_attack(theta, inflow.Wm, inflow.Wtheta)

    # - Look up lfit and drag data - #
    cl = zeros(length(alpha))
    cd = zeros(length(alpha))
    for a in 1:length(alpha)
        cl[a], cd[a] = dt.search_polars(af, alpha[a])
    end

    # - swap old to new - #
    Gamma_temp .= Gamma

    # - Calculate Circulation from lift and inflow - #
    @. Gamma = 0.5 * inflow.Wmag * chord * cl
    @. rotor_source_strengths = B * inflow.Wmag * chord * cd / (4.0 * pi * rpc)

    # get net circulation and B*Gamma
    BGamma = Gamma * B
    Gamma_tilde = B * Gamma #only because we only have one rotor

    # - Calculate Enthalpy Jumps - #
    H_tilde = dt.calculate_enthalpy_jumps(Gamma, Omega, B)

    # - Calculate Vm in wake - #
    Vmavg_wake = get_wake_vmavg(
        Vinf,
        wake_vortex_strengths,
        vxd_wake_to_wake,
        vrd_wake_to_wake,
        rotor_source_strengths,
        vxd_rotor_to_wake,
        vrd_rotor_to_wake,
    )

    # - Calculate gamma_thetas - #
    wake_vortex_strengths .= gamma_theta_from_vmavg(
        Vmavg_wake, Gamma_tilde, H_tilde, rotor_panel_edges
    )

    # - PLOT - #
    if iter[1]%10 == 0
    plot!(pcirc, Gamma, rpc; label="iter #$(iter[1])")
    plot!(psigma, rotor_source_strengths, rpc; label="iter #$(iter[1])")
    plot!(palpha, alpha * 180.0 / pi, rpc; label="iter #$(iter[1])")
    plot!(pvm, inflow.Wm .- Vinf, rpc; label="iter #$(iter[1])")
    plot!(pvt, inflow.Wtheta.+Omega.*rpc, rpc; label="iter #$(iter[1])")
    plot!(pw, inflow.Wmag, rpc; label="iter #$(iter[1])")
    plot!(pcl,cl, rpc; label="iter #$(iter[1])")
    plot!(pcd, cd, rpc; label="iter #$(iter[1])")
    plot!(pgammatheta, wake_vortex_strengths[:,41], rotor_panel_edges; label="iter #$(iter[1])")
    end
end

## -- Save Figures -- ##
savefig(pcirc, "test/manual_tests/rotor_wake_tests/circulation_with_sources.pdf")
savefig(psigma, "test/manual_tests/rotor_wake_tests/rotor_source_strengths_with_sources.pdf")
savefig(palpha, "test/manual_tests/rotor_wake_tests/alpha_with_sources.pdf")
savefig(pvm, "test/manual_tests/rotor_wake_tests/vm_with_sources.pdf")
savefig(pvt, "test/manual_tests/rotor_wake_tests/vtheta_with_sources.pdf")
savefig(pw, "test/manual_tests/rotor_wake_tests/Wmag_with_sources.pdf")
savefig(pcl, "test/manual_tests/rotor_wake_tests/cl_with_sources.pdf")
savefig(pcd, "test/manual_tests/rotor_wake_tests/cd_with_sources.pdf")
savefig(
    pgammatheta, "test/manual_tests/rotor_wake_tests/wake_vortex_strength_iteration.pdf"
)

plot(ccbGamma, r; xlabel=L"\Gamma", ylabel="r", label="CCBlade")
plot!(Gamma, rpc; label="DuctTAPE")
savefig("test/manual_tests/rotor_wake_tests/final_circulation.pdf")
