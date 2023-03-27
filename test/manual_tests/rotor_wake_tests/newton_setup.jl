
#---------------------------------#
#             Includes            #
#---------------------------------#

# - Include Plots and my default settings - #
include("../../../plots_default.jl")
# include("plots_default.jl")

# - use CCBlade for comparison - #
using CCBlade

# - Rename DuctTAPE for convenience - #
using DuctTAPE
const dt = DuctTAPE

# - Rename FLOWFoil for convenience - #
using FLOWFoil
const ff = FLOWFoil

#SOLVER STUFF
using NLsolve
using ImplicitAD

using LineSearches: BackTracking

using SpecialFunctions

function setup_stuff()
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
    Note that the panels on the rotor are the source panels modeling the blade element profile drag, centered at the blade element locations for the most part.
    Again, the sections close to the hub and tip don't perfectly align with the panel centers.
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

    #get rotor panel radial center points for convenience
    rpc = rotor_panels.panel_center[:, 2]

    # ##### ----- SANITY PLOTS ----- #####

    # # - plot the original radial stations - #
    # plot(
    #     zeros(nbe),
    #     r;
    #     seriestype=:scatter,
    #     markershape=:square,
    #     markersize=3,
    #     label="Input Blade Element Positions",
    # )

    # # - plot the panel edges - #
    # plot!(
    #     zeros(nbe + 1),
    #     rotor_panel_edges;
    #     markershape=:rect,
    #     markersize=1,
    #     linewidth=0.5,
    #     label="Rotor Panels",
    # )

    # # - Plot panel centers - #
    # plot!(
    #     rotor_panels.panel_center[:, 1],
    #     rotor_panels.panel_center[:, 2];
    #     seriestype=:scatter,
    #     markershape=:circle,
    #     markersize=1.5,
    #     label="Panel Centers",
    # )

    # # - Plot Rtip and Rhub - #
    # plot!([-0.1; 0.1], [Rtip; Rtip]; linestyle=:dash, color=:black, label="Rtip and Rhub")
    # plot!([-0.1; 0.1], [Rhub; Rhub]; linestyle=:dash, color=:black, label="")

    # # - save figure - #
    # savefig("test/manual_tests/rotor_wake_tests/sanity_check_rotor_geometry.pdf")

    #---------------------------------#
    #              Wake               #
    #---------------------------------#
    #=
    Trailing wake sheets extend from between the blade elements, in other words, from the rotor panel edge locations.
    =#

    ##### ----- General Geometry ----- #####

    # - Get blade element spacing - #
    #arbitrarily pick spacing close to radial spacing?
    # dr = r[end - 1] - r[end - 2]
    #or just pick something small? (in practice, just match with the body panels
    dr = 0.001

    # - Rotor Diameter - #
    D = 2.0 * Rtip

    # - Define x grid spacing - #
    #note that by using step rather than length, we aren't guarenteed to have the wake as long as we want, but it will be close.
    xrange = range(0.0, 5.0 * D; step=dr)

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
    Thus the wake_vortex_panels object is a vector of Panel structs, where each struct is comprised of a single row (sheet) of wake panels.
    We will need to remember this in indexing such that we call wake_vortex_panels[row].property[column]
    Note that there are nbe+1 trailing wake surfaces
    =#
    wake_vortex_panels = dt.generate_wake_panels(
        x_grid_points, r_grid_points; method=wake_method
    )

    # ##### ----- SANITY PLOTS ----- #####

    # # - plot the original radial stations - #
    # for i in 1:length(x_grid_points[1, :])
    #     if i == 1
    #         plot(
    #             x_grid_points[:, i],
    #             r_grid_points[:, i];
    #             markershape=:rect,
    #             markersize=0.5,
    #             linewidth=0.25,
    #             color=mycolors[2],
    #             aspectratio=:equal,
    #             label="Input Wake Grid (panel edges)",
    #         )
    #     else
    #         plot!(
    #             x_grid_points[:, i],
    #             r_grid_points[:, i];
    #             markershape=:rect,
    #             markersize=0.5,
    #             linewidth=0.25,
    #             color=mycolors[2],
    #             aspectratio=:equal,
    #             label="",
    #         )
    #     end
    # end

    # # - Plot panel centers - #
    # for i in 1:length(wake_vortex_panels)
    #     if i == 1
    #         plot!(
    #             wake_vortex_panels[i].panel_center[:, 1],
    #             wake_vortex_panels[i].panel_center[:, 2];
    #             seriestype=:scatter,
    #             markershape=:circle,
    #             markersize=0.5,
    #             color=mycolors[3],
    #             label="Panel Centers",
    #         )
    #     else
    #         plot!(
    #             wake_vortex_panels[i].panel_center[:, 1],
    #             wake_vortex_panels[i].panel_center[:, 2];
    #             seriestype=:scatter,
    #             markershape=:circle,
    #             markersize=0.5,
    #             color=mycolors[3],
    #             label="",
    #         )
    #     end
    # end

    # # - plot where 2D point is - #
    # plot!(
    #     [2.0 * D; 2.0 * D], [Rhub, Rtip]; linestyle=:dash, color=:black, label="2D Downstream"
    # )

    # # - save figure - #
    # savefig("test/manual_tests/rotor_wake_tests/sanity_check_wake_geometry.pdf")

    ######################################################################
    #                                                                    #
    #                  PRECOMPUTE UNIT INDUCED VELOCITIES                #
    #                                                                    #
    ######################################################################

    #---------------------------------#
    #             MESHES              #
    #---------------------------------#
    #=
    Meshes contain the relative geometry used to calculate the unit induced velocities of the panels doing the influencing to the panels being affected.
    the first input is the influencing panel arrays
    the second input is the affected panel arrays
    the output is a matrix of [influence index, affect index] for the various panel arrays
    =#

    ##### ----- Wake to Rotor ----- #####
    wake_to_rotor_mesh = [
        dt.generate_one_way_mesh(wake_vortex_panels[i], rotor_panels) for
        i in 1:length(wake_vortex_panels)
    ]

    ##### ----- Wake to Wake ----- #####

    #---------------------------------#
    #           Velocities            #
    #---------------------------------#
    #=
    To get the induced velocities, we take the meshes just created, then input the influencing panel arrays again, then the affected panel arrays.
    The output is a matrix again with index [influencing, affected]

    We then separate out the inputs into the x and r components
    =#

    ##### ----- Wake to Rotor ----- #####
    A_wake_to_rotor = [
        dt.assemble_induced_velocity_matrices(
            wake_to_rotor_mesh[i], wake_vortex_panels[i], rotor_panels
        ) for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    # - Axial - #
    vxd_wake_to_rotor = [
        A_wake_to_rotor[i,j][1] for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    # - Radial - #
    vrd_wake_to_rotor = [
        A_wake_to_rotor[i,j][2] for i in 1:1, j in 1:length(wake_vortex_panels)
    ]

    ##### ----- Wake to Wake ----- #####

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

    #################################

    # - Initialize with freestream only - #
    W = sqrt.(Vinf^2 .+ (Omega * rpc) .^ 2)
    cl, cd = dt.search_polars(af, atan.(Vinf, Omega * rpc))
    Gamma_init = 0.5 .* W .* cl .* chord
    Vmr = Vinf * ones(length(W))

    Gamma_tilde_init = B .* Gamma_init
    H_tilde_init = Omega * B * Gamma_init / (2.0 * pi)

    # # - Initialize from freestream - #
    # gamma_wake_init = zeros(length(Gamma_init) + 1, 1)
    # dt.calculate_wake_vortex_strengths!(
    #     gamma_wake_init, rotor_panel_edges, Vmr, Gamma_tilde_init, H_tilde_init
    # )
    # println("gamma wake init: ")
    # display(gamma_wake_init)

    function vm_from_vinf(Vinf, Gamma_tilde, H_tilde, radial_stations)

        #rename for convenience
        nr = length(Gamma_tilde)

        #initialize
        Vm = zeros(nr)

        #march from just outside the tip to the hub
        for i in nr:-1:1
            if i == nr

                #this is the change at the tip, so we need to set Vm2 to Vinf, and set Gamma dn H 2's to zeros
                radical =
                    Vinf^2 +
                    (1.0 / (2.0 * pi * radial_stations[i]))^2 * (-Gamma_tilde[i]^2) -
                    2.0 * (-H_tilde[i])
                # if radical > 0.0
                Vm[i] = sqrt(radical)
                # else
                #     Vm[i] = 0.0
                # end
            else

                #otherwise, we just take the differences inside as-is
                radical =
                    Vm[i + 1]^2 +
                    (1.0 / (2.0 * pi * radial_stations[i]))^2 *
                    (Gamma_tilde[i + 1]^2 - Gamma_tilde[i]^2) -
                    2.0 * (H_tilde[i + 1] - H_tilde[i])
                # if radical > 0.0
                Vm[i] = sqrt(radical)
                # else
                #     Vm[i] = 0.0
                # end
            end
        end

        #Vm should now be in the order of hub to tip
        return Vm
    end
    # - initialize from open rotor - #
    Vm_init = vm_from_vinf(Vinf, Gamma_tilde_init, H_tilde_init, rpc)
    gamma_wake_init = [
        Vm_init[1] - Vinf
        Vm_init[2:end] .- Vm_init[1:(end - 1)]
        Vinf - Vm_init[end]
    ]

    #### PLOTS ####

    #Plot vx comparison
    pu = plot(
        out.u,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="u, vx",
        ylabel="r",
    )

    #plot vtheta comparisons
    pv = plot(
        out.v,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="v, vtheta",
        ylabel="r",
    )

    #plot inflow velocity comparison
    pW = plot(
        out.W,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="W",
        ylabel="r",
    )

    #plot angle of attack comparison
    pa = plot(
        out.alpha,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="alpha",
        ylabel="r",
    )

    #plot cl comparison
    pcl = plot(
        out.cl,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="cl",
        ylabel="r",
    )

    #plot cd comparison
    pcd = plot(
        out.cd,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="cd",
        ylabel="r",
    )

    #plot circulation comparison
    pg = plot(
        ccbGamma,
        r;
        label="CCBlade",
        linewidth=2,
        linestyle=:dash,
        color=:black,
        xlabel="Gamma",
        ylabel="r",
    )

    pgw = plot(
        gamma_wake_init,
        rotor_panel_edges;
        xlabel="wake gammas",
        ylabel="r",
        label="dVm Init",
        linewidth=2,
        linestyle=:dash,
        color=:black,
    )

    pvm = plot(
        Vmr,
        rpc;
        xlabel="Vm",
        ylabel="r",
        linestyle=:dash,
        linewidth=2,
        label="Vm march init",
    )

    # - Set Up States and Parameters - #
    states = [Gamma_init; gamma_wake_init]
    # states = [ccbGamma; gamma_wake_init]
    params = (
        pu=pu,
        pv=pv,
        pW=pW,
        pa=pa,
        pcl=pcl,
        pcd=pcd,
        pg=pg,
        pgw=pgw,
        pvm=pvm,
        plotting=true,
        converged=[false],
        Gammaidx=1:length(ccbGamma),
        gamma_theta_idx=(length(ccbGamma) + 1):(length(gamma_wake_init) + length(ccbGamma)),
        nxwake=length(xrange) - 1,
        vx_rw=vxd_wake_to_rotor,
        vr_rw=vrd_wake_to_rotor,
        blade_elements=[(
            omega=Omega,
            num_blades=B,
            chords=chord,
            twists=theta,
            airfoils=fill(af, length(rpc)),
            radial_positions=rpc,
        )],
        num_rotors=1,
        num_blades=B,
        Vinf=Vinf,
        Omega=Omega,
        iter=[1],
        af=af,
        chord=chord,
        twist=theta,
        rotor_panel_edges=rotor_panel_edges,
        rotor_panel_centers=rpc,
    )

    return states, params
end
