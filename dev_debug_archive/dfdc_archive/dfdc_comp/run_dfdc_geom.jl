
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

savepath = project_dir * "/dev_debug_archive/dfdc_comp/"

include(savepath * "gather_dfdc_data.jl")
include(savepath * "run_dfdc_example.jl")

using DuctAPE
const dt = DuctAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

using NLsolve

function init(dfdc)
    #---------------------------------#
    #         Rotor Geometry          #
    #---------------------------------#
    #=
    ROTOR
    !       Xdisk        Nblds       NRPdef
      0.12000                5           11
    !  #stations
        10
    !           r        Chord         Beta
      0.50491E-01  0.89142E-01   69.012
      0.61567E-01  0.79785E-01   59.142
      0.72644E-01  0.71300E-01   51.825
      0.83721E-01  0.63979E-01   46.272
      0.94798E-01  0.57777E-01   41.952
      0.10587      0.52541E-01   38.509
      0.11695      0.48103E-01   35.699
      0.12803      0.44316E-01   33.354
      0.13911      0.41061E-01   31.349
      0.15018      0.38243E-01   29.596
    ENDROTOR
    =#

    xrotor = 0.12
    B = 5

    # rotor panel nodes
    rotor_panel_edges = dfdc.rotor_node_r
    nr = length(rotor_panel_edges) - 1 #number of blade elements

    # rotor chord lengths
    chords = dfdc.chord

    #twist in degrees converted to radians
    twists = dfdc.twist

    #local solidity
    solidity = dfdc.solidity

    # local stagger
    stagger = dfdc.stagger

    #Airfoil Data:
    alpha0 = 0.0
    clmax = 1.5
    clmin = -1.0
    dclda = 2.0 * pi
    dclda_stall = 0.5
    dcl_stall = 0.2
    cdmin = 0.012
    cldmin = 0.1
    dcdcl2 = 0.005
    cmcon = 0.0
    Re_ref = 2e5
    Re_exp = 0.35
    mcrit = 0.7

    afparams = dt.DFDCairfoil(
        alpha0,
        clmax,
        clmin,
        dclda,
        dclda_stall,
        dcl_stall,
        cdmin,
        cldmin,
        dcdcl2,
        cmcon,
        Re_ref,
        Re_exp,
        mcrit,
    )

    outer_airfoil = fill(afparams, nr)
    inner_airfoil = fill(afparams, nr)
    inner_fraction = ones(nr)

    #---------------------------------#
    #       Duct and Hub Geometry     #
    #---------------------------------#

    duct_coordinates = [reverse(dfdc.duct_node_x) reverse(dfdc.duct_node_r)]
    hub_coordinates = [reverse(dfdc.hub_node_x) reverse(dfdc.hub_node_r)]

    _, duct_leidx = findmin(duct_coordinates[:, 1])
    ductxin = reverse(duct_coordinates[1:duct_leidx, 1])
    ductrin = reverse(duct_coordinates[1:duct_leidx, 2])

    # load in duct and hub geometry, spline, and find out what the duct and hub radii are at the rotor positions to figure out what Rtip and Rhub are.
    Rhub = FLOWMath.akima(hub_coordinates[:, 1], hub_coordinates[:, 2], xrotor)
    Rtip = FLOWMath.akima(ductxin, ductrin, xrotor)

    #---------------------------------#
    #        Wake Coordinates         #
    #---------------------------------#

    inner_wake_x = dfdc.wake_node_x
    inner_wake_r = dfdc.wake_node_r
    hub_wake_x = dfdc.hubwake_node_x
    hub_wake_r = dfdc.hubwake_node_r
    duct_wake_x = dfdc.ductwake_node_x
    duct_wake_r = dfdc.ductwake_node_r
    duct_wall_x = dfdc.ductinterface_node_x
    hub_wall_x = dfdc.hubinterface_node_x
    duct_wall_r = dfdc.ductinterface_node_r
    hub_wall_r = dfdc.hubinterface_node_r
    wake1_x = [reverse(hub_wall_x); hub_wake_x[2:end]]
    wakeend_x = [duct_wall_x; duct_wake_x[2:end]]
    wake1_r = [reverse(hub_wall_r); hub_wake_r[2:end]]
    wakeend_r = [duct_wall_r; duct_wake_r[2:end]]

    xgrid = [wake1_x'; inner_wake_x; wakeend_x']
    rgrid = [wake1_r'; inner_wake_r; wakeend_r']

    wake_coordinates = (; xgrid, rgrid)

    #---------------------------------#
    #      Operating Conditions       #
    #---------------------------------#
    #=
    OPER
    !        Vinf         Vref          RPM
       0.0000       50.000       8000.0
    !         Rho          Vso          Rmu           Alt
       1.2260       340.00      0.17800E-04   0.0000
    !       XDwake        Nwake
      0.80000               20
    !       Lwkrlx
                F
    ENDOPER
    =#

    Vinf = 0.0
    Vref = 50.0
    RPM = 8000
    Omega = RPM * pi / 30  # convert from RPM to rad/s
    rho = 1.226
    asound = 340.0
    mu = 1.78e-5

    #--------------------------------#
    #      Assemble Named Tuples      #
    #---------------------------------#

    # Rotor Parameters
    rotor_parameters = [(;
        xrotor,
        rotor_panel_edges,
        chords,
        twists,
        stagger,
        solidity,
        outer_airfoil,
        inner_airfoil,
        inner_fraction,
        Rtip,
        Rhub,
        tip_gap=0.0,
        B,
        Omega,
    )]

    # Freestream Parameters
    freestream = (; rho, mu, asound, Vinf)

    reference_parameters = (; Vref, Rref=Rtip)

    return duct_coordinates,
    hub_coordinates, wake_coordinates, rotor_parameters, freestream,
    reference_parameters
end

function plot_geom(inputs)
    pg = plot(; xlabel="x", ylabel="r", aspectratio=1)

    # - Plot Duct - #
    plot!(
        pg,
        inputs.duct_coordinates[:, 1],
        inputs.duct_coordinates[:, 2];
        markershape=:rect,
        markersize=1.0,
        linewidth=0.5,
        color=mycolors[1],
        label="body coordinates",
    )

    plot!(
        pg,
        inputs.body_panels[1].panel_center[:, 1],
        inputs.body_panels[1].panel_center[:, 2];
        seriestype=:scatter,
        markersize=1.5,
        color=mycolors[4],
        label="body panel centers",
    )

    # - Plot Hub - #
    plot!(
        pg,
        inputs.hub_coordinates[:, 1],
        inputs.hub_coordinates[:, 2];
        markershape=:rect,
        markersize=1.0,
        linewidth=0.5,
        color=mycolors[1],
        label="",
    )

    plot!(
        pg,
        inputs.body_panels[2].panel_center[:, 1],
        inputs.body_panels[2].panel_center[:, 2];
        seriestype=:scatter,
        markersize=1.5,
        color=mycolors[4],
        label="",
    )

    # - Plot Rotor - #
    plot!(
        pg,
        inputs.blade_elements[1].xrotor * ones(length(inputs.rotor_panel_edges)),
        inputs.rotor_panel_edges;
        linewidth=0.5,
        markershape=:rect,
        markersize=1.0,
        color=mycolors[2],
        label="rotor coordinates",
    )

    plot!(
        pg,
        inputs.blade_elements[1].xrotor * ones(length(inputs.rotor_panel_centers)),
        inputs.rotor_panel_centers;
        seriestype=:scatter,
        markersize=1.5,
        color=mycolors[7],
        label="rotor panel centers",
    )

    # - Plot Wakes - #
    for i in 1:length(inputs.wake_vortex_panels)
        lab1 = i == 1 ? "wake coordinates" : ""
        lab2 = i == 1 ? "wake panel centers" : ""

        plot!(
            pg,
            inputs.wakexgrid[i, :],
            inputs.wakergrid[i, :];
            markershape=:rect,
            markersize=1.0,
            linewidth=0.5,
            color=mycolors[3],
            label=lab1,
        )

        plot!(
            pg,
            inputs.wake_vortex_panels[i].panel_center[:, 1],
            inputs.wake_vortex_panels[i].panel_center[:, 2];
            seriestype=:scatter,
            markersize=1.5,
            color=mycolors[6],
            label=lab2,
        )
    end

    savefig(savepath * "dfdc-manual-init-geometry.pdf")

    return nothing
end

function run_geom()
    dfdc = get_dfdc()

    duct_coordinates, hub_coordinates, wake_coordinates, rotor_parameters, freestream, reference_parameters = init(
        dfdc
    )

    inputs = dt.manual_precomputed_inputs(
        duct_coordinates,      # panel node locations
        hub_coordinates,       # panel node locations
        wake_coordinates,      # tuple containing xgrid[x,r], rgrid[x,r]
        rotor_parameters,      # tuple with rotor paramters
        freestream,            # tuple containing Vinf, rho, mu, asound
        reference_parameters;  # tuple containing Vref, Rref
        debug=false,
    )

    plot_geom(inputs)

    states, initials, conv = dt.analyze_propulsor(inputs)
    out = dt.post_process(states, inputs)

    cdump = dt.dump(states, inputs)

    println("Plotting Rotor/Wake States")
    plotstates(cdump.Gamr, cdump.sigr, cdump.gamw, inputs, true)

    println("Plotting Blade Velocities")
    plotbladevelocities(cdump, inputs)

    println("Plotting Blade Angles")
    plotangles(cdump, inputs)

    println("Plotting Lift/Drag on Blade")
    plotclcd(cdump, inputs)

    println("Plotting Enthalpy and Entropy Disk Jumps")
    ploths(cdump, inputs)

    println("Plotting Body Aerodynamics")
    plotbodyaero(states, inputs, out, true)

    return nothing
end
