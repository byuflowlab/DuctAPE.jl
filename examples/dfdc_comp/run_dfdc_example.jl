
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

using DuctTAPE
const dt = DuctTAPE

# CCBlade used for it's airfoils function objects here.
using CCBlade
const ccb = CCBlade

using FLOWMath
const fm = FLOWMath

#---------------------------------#
# - Include DFDC DATA for Comparison
#---------------------------------#
#surface velocity and pressure
include(project_dir * "/examples/dfdc_comp/dfdccp.jl")
#test to make sure you have the right indices
include(project_dir * "/examples/dfdc_comp/CPR_duct.jl")
# a bunch of the things to compare from dfdc,
include(project_dir * "/examples/dfdc_comp/ALL_THE_STUFF.jl")
# final printout stuff
# r/R    c/R     beta deg alfa     CL     CD    REx10^3   Mach    B*Gam
include(project_dir * "/examples/dfdc_comp/DFDC_printout.jl")

function runandplot()
    println("Setting Up")
    rotor_parameters, paneling_constants, freestream, duct_coords, hub_coords, Vref = init()

    println("Running Analysis")
    converged_states, inputs, initial_states, convergeflag = rundt(
        duct_coords, hub_coords, rotor_parameters, paneling_constants, freestream
    )

    println("Plotting Geometry")
    plotgeom(inputs, paneling_constants)


    println("Dumping Everything")
    cdump = dt.dump(converged_states, inputs)

    println("Plotting Rotor/Wake States")
    plotstates(cdump.Gamr, cdump.sigr, cdump.gamw, inputs, convergeflag)

    println("Plotting Blade Velocities")
    plotbladevelocities(cdump, inputs)

    println("Plotting Blade Angles")
    plotangles(cdump, inputs)

    println("Plotting Lift/Drag on Blade")
    plotclcd(cdump, inputs)

    println("Plotting Enthalpy and Entropy Disk Jumps")
    ploths(cdump, inputs)

    println("Plotting Body Aerodynamics")
    plotbodyaero(converged_states, inputs, initial_states, convergeflag, Vref)

    return cdump
end

function init()
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

    #dimensional radius
    r = [
        0.50491E-01
        0.61567E-01
        0.72644E-01
        0.83721E-01
        0.94798E-01
        0.10587
        0.11695
        0.12803
        0.13911
        0.15018
    ]

    #dimensional chord
    chords = [
        0.89142E-01
        0.79785E-01
        0.71300E-01
        0.63979E-01
        0.57777E-01
        0.52541E-01
        0.48103E-01
        0.44316E-01
        0.41061E-01
        0.38243E-01
    ]

    #twist in degrees converted to radians
    twists =
        [
            69.012
            59.142
            51.825
            46.272
            41.952
            38.509
            35.699
            33.354
            31.349
            29.596
        ] * pi / 180.0

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

    # airfoils = fill(ccb.AlphaAF("examples/dfdc_comp/dfdc_polar.dat"), length(r))
    airfoils = fill(afparams, length(r))

    #---------------------------------#
    #       Duct and Hub Geometry     #
    #---------------------------------#

    include(project_dir * "/examples/dfdc_comp/dfdc_ductgeom.jl")
    _, duct_leidx = findmin(duct_coords[:, 1])
    ductxin = reverse(duct_coords[1:duct_leidx, 1])
    ductrin = reverse(duct_coords[1:duct_leidx, 2])

    # load in duct and hub geometry, spline, and find out what the duct and hub radii are at the rotor positions to figure out what Rtip and Rhub are.
    Rhub = FLOWMath.akima(hub_coords[:, 1], hub_coords[:, 2], xrotor)
    Rtip = FLOWMath.akima(ductxin, ductrin, xrotor)

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

    #---------------------------------#
    #        Paneling Options         #
    #---------------------------------#

    wake_length = 0.8 #times duct chord?

    nwake_sheets = 11 #note that nwake in the dfdc file is how many panels to have in the wake

    npanels_inlet = 10
    discscale = 1.0

    ductle = minimum(duct_coords[:, 1])
    ductte = maximum(duct_coords[:, 1])
    ductchord = maximum(duct_coords[:, 1]) - minimum(duct_coords[:, 1])
    outletinletratio = (ductte - xrotor) / (xrotor - ductle)

    nhub_inlet = round(Int, npanels_inlet * discscale)

    nduct_inlet = round(Int, npanels_inlet * discscale)

    nduct_outlet = round(Int, nduct_inlet * outletinletratio)

    nwake = round(Int, (nduct_inlet + nduct_outlet) * wake_length)

    npanels = [nduct_outlet, nwake]

    #--------------------------------#
    #      Assemble Named Tuples      #
    #---------------------------------#

    # Rotor Parameters
    rotor_parameters = [(;
        xrotor=xrotor,
        nwake_sheets,
        r=r ./ Rtip, #non-dimensionalize
        chords,
        twists,
        airfoils,
        Rtip,
        Rhub,
        tip_gap=0.0,
        B,
        Omega,
    )]

    # Paneling Parameters
    paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

    # Freestream Parameters
    freestream = (; rho, mu, asound, Vinf)

    return rotor_parameters, paneling_constants, freestream, duct_coords, hub_coords, Vref
end

function rundt(duct_coords, hub_coords, rotor_parameters, paneling_constants, freestream)
    #---------------------------------#
    #           Run Solver            #
    #---------------------------------#

    converged_states, inputs, initial_states, convergeflag = dt.analyze_propulsor(
        duct_coords,
        hub_coords,
        paneling_constants,
        rotor_parameters,
        freestream;
        debug=false,
        maximum_linesearch_step_size=1e6,
        iteration_limit=100,
    )

    return converged_states, inputs, initial_states, convergeflag
end

function plotgeom(inputs, paneling_constants)
    ##### ----- Plot GEOMETRY ----- #####
    #initialize plot
    pgeom = plot(; aspectratio=1, xlabel="x", ylabel="r")
    plot!(
        pgeom,
        inputs.blade_elements[1].xrotor * ones(length(inputs.rotor_panel_edges)),
        inputs.rotor_panel_edges;
        color=mycolors[2],
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        label="",
    )

    plot!(
        pgeom,
        inputs.rotor_source_panels[1].panel_center[:, 1],
        inputs.rotor_source_panels[1].panel_center[:, 2];
        color=mycolors[3],
        seriestype=:scatter,
        markersize=0.75,
        markershape=:circle,
        label="",
    )

    for iw in 1:(paneling_constants.nwake_sheets)
        plot!(
            pgeom,
            inputs.wakexgrid[:, iw],
            inputs.wakergrid[:, iw];
            linewidth=0.25,
            markersize=0.5,
            markershape=:rect,
            color=:black,
            label="",
        )

        plot!(
            pgeom,
            inputs.wake_vortex_panels[iw].panel_center[:, 1],
            inputs.wake_vortex_panels[iw].panel_center[:, 2];
            seriestype=:scatter,
            markersize=0.75,
            markershape=:circle,
            color=mycolors[2],
            label="",
        )
    end

    savefig(project_dir * "/examples/dfdc_comp/precomputed-rotor-wake-geometry.pdf")

    # plot body panels
    plot!(
        pgeom,
        inputs.duct_coordinates[:, 1],
        inputs.duct_coordinates[:, 2];
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        color=mycolors[3],
        label="",
    )

    plot!(
        pgeom,
        inputs.hub_coordinates[:, 1],
        inputs.hub_coordinates[:, 2];
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        color=mycolors[3],
        label="",
    )

    # Plot body panel centers
    for ib in 1:2
        # for ib in 1:1
        plot!(
            pgeom,
            inputs.body_panels[ib].panel_center[:, 1],
            inputs.body_panels[ib].panel_center[:, 2];
            color=mycolors[1],
            seriestype=:scatter,
            markersize=0.75,
            label="",
        )
    end

    savefig(pgeom, project_dir * "/examples/dfdc_comp/precomputed-full-geometry.pdf")

    return nothing
end

function plotstates(Gamr, sigr, gamw, inputs, convergeflag)

    #TODO: need sigr and gamw from dfdc

    if convergeflag
        convlabel = "Converged"
    else
        convlabel = "NOT converged"
    end

    ##### ----- Plot rotor circulation distribution ----- #####
    # initialize plot
    pG = plot(; xlabel=L"\Gamma", ylabel="r")

    # plot solution
    plot!(pG, Gamr, inputs.rotor_panel_centers; label=convlabel)
    plot!(pG, dfdcGamr, dfdcr; label="DFDC Converged")

    #save
    savefig(pG, project_dir * "/examples/dfdc_comp/rotorcirculation-check.pdf")

    ##### ----- Plot Wake Strengths ----- #####

    #pg = plot(; xlabel=L"\gamma_\theta^{wake}", ylabel="r")

    ## plot solution
    #plot!(pg, gamw, inputs.rotor_panel_edges; label=convlabel)

    ##save
    #savefig(pg, project_dir * "/examples/dfdc_comp/wake-strength.pdf")

    ##### ----- Plot Source Strengths ----- #####
    ps = plot(; xlabel=L"\sigma", ylabel="r")

    # plot solution
    plot!(ps, sigr, inputs.rotor_panel_centers; label=convlabel)
    plot!(ps, dfdcsigma, dfdcr; label="DFDC")

    #save
    savefig(ps, project_dir * "/examples/dfdc_comp/source-strength-check.pdf")

    return nothing
end

function plotbodyaero(cstates, inputs, istates, convergeflag, Vref)
    if convergeflag
        convlabel = "Converged"
    else
        convlabel = "NOT converged"
    end

    #---------------------------------#
    #    Extract Outputs and Plot     #
    #---------------------------------#
    gamb, gamw, Gamr, sigr = dt.extract_state_variables(cstates, inputs)

    ##### ----- Plot duct surface velocity ----- #####

    #prepare outputs
    dp = inputs.body_panels[1].panel_center[:, 1]
    hp = inputs.body_panels[2].panel_center[:, 1]
    _, leidx = findmin(dp)
    #split into inner and outer surfaces
    dpinner = dp[1:leidx]
    dpouter = dp[(leidx + 1):end]

    gamdinner = gamb[1:leidx] / Vref
    gamdouter = gamb[(leidx + 1):length(dp)] / Vref
    gamh = gamb[(length(dp) + 1):end]

    # initialize plot
    pb = plot(; xlabel="x", ylabel="Q/Qref")

    # plot solution
    # plot!(pb, dpinner, gamdinner; color=mycolors[1], label=convlabel * " inner surface, with rotor")
    # plot!(pb, dpouter, gamdouter; label=convlabel * " outer surface, with rotor")

    plot!(
        pb,
        dp,
        abs.(gamb[1:length(dp)]) ./ Vref;
        color=mycolors[1],
        label=convlabel * " DuctTAPE Duct",
    )
    plot!(
        pb,
        duct_cpv[:, 1],
        duct_cpv[:, 4] ./ Vref;
        color=mycolors[1],
        linestyle=:dash,
        label="DFDC Duct",
    )

    plot!(
        pb,
        hp,
        abs.(gamb[(length(dp) + 1):end]) ./ Vref;
        color=mycolors[2],
        label=convlabel * " DuctTAPE Hub",
    )
    plot!(
        pb,
        hub_cpv[:, 1],
        hub_cpv[:, 4] ./ Vref;
        color=mycolors[2],
        linestyle=:dash,
        label="DFDC Hub",
    )

    #plot rotor location
    plot!(
        pb,
        inputs.blade_elements[1].xrotor * ones(2),
        abs.([0.0; maximum(duct_cpv[:, 4] ./ Vref)]);
        # linewidth=0.25,
        linestyle=:dot,
        color=mycolors[3],
        label="rotor location",
    )

    #save
    savefig(pb, project_dir * "/examples/dfdc_comp/body-velocity.pdf")

    ##### ----- Plot Surface Pressure ----- #####

    cpductinner, cpductouter, cphub, xdi, xdo, xh = dt.get_cps(cstates, inputs, Vref)

    pcp = plot(; xlabel="x", ylabel=L"C_p", yflip=true)

    # plot solution
    # plot!(pcp, xdi, cpductinner; color=mycolors[1], label=convlabel * " inner duct surface")

    # plot!(pcp, xdo, cpductouter; label=convlabel * " outer duct surface")

    plot!(
        pcp,
        [xdi; xdo],
        [cpductinner; cpductouter];
        color=mycolors[1],
        label=convlabel * " duct surface",
    )

    plot!(
        pcp,
        duct_cpv[:, 1],
        duct_cpv[:, 3];
        color=mycolors[1],
        linestyle=:dash,
        label="DFDC Duct",
    )
    plot!(
        pcp,
        duct_cpv[:, 1],
        CPR_DUCT;
        color=mycolors[4],
        linestyle=:dot,
        label="DFDC Duct sanity check on indices",
    )

    plot!(pcp, xh, cphub; color=mycolors[2], label=convlabel * " hub surface")
    plot!(
        pcp,
        hub_cpv[:, 1],
        hub_cpv[:, 3];
        color=mycolors[2],
        linestyle=:dash,
        label="DFDC Hub",
    )

    #save
    savefig(pcp, project_dir * "/examples/dfdc_comp/body-pressure.pdf")

    return nothing
end

function plotbladevelocities(cdump, inputs)

    # - extract ducttape radial positions - #
    rpc = inputs.rotor_panel_centers
    omega = inputs.blade_elements[1].Omega

    # - initialize plots - #
    pvx = plot(; xlabel="axial velocities", ylabel="radial positions")
    pvr = plot(; xlabel="radial velocities", ylabel="radial positions")
    pvt = plot(; xlabel="tangential velocities", ylabel="radial positions")
    pw = plot(; xlabel="inflow magnitude", ylabel="radial positions")
    por = plot(; xlabel=L"\Omega r", ylabel="radial positions")

    plot!(pvx, cdump.Wx_rotor, rpc; label="Wx")
    plot!(pvx, cdump.vx_rotor, rpc; label="vx")
    plot!(pvx, cdump.vxfrombody, rpc; label="vxi from body")
    plot!(pvx, cdump.vxfromwake, rpc; label="vxi from wake")
    plot!(pvx, cdump.vxfromrotor, rpc; label="vxi from rotor")
    plot!(pvx, dfdcvx, dfdcr; label="DFDC Wx")
    savefig(pvx, project_dir * "/examples/dfdc_comp/VXcomp.pdf")

    plot!(pvr, cdump.vr_rotor, rpc; label="Wr")
    plot!(pvr, cdump.vrfrombody, rpc; label="vri from body")
    plot!(pvr, cdump.vrfromwake, rpc; label="vri from wake")
    plot!(pvr, cdump.vrfromrotor, rpc; label="vri from rotor")
    plot!(pvr, dfdcvr, dfdcr; label="DFDC Wr")
    savefig(pvr, project_dir * "/examples/dfdc_comp/VRcomp.pdf")

    plot!(pvt, cdump.Wtheta_rotor, rpc; label=L"W_\theta")
    plot!(pvt, cdump.vtheta_rotor, rpc; label=L"v_\theta")
    plot!(pvt, dfdcvtheta, dfdcr; label=L"DFDC~v_\theta")
    plot!(pvt, dfdcwtheta, dfdcr; label=L"DFDC~W_\theta")
    savefig(pvt, project_dir * "/examples/dfdc_comp/Vthetacomp.pdf")

    plot!(pw, cdump.Wmag_rotor, rpc; label="W")
    plot!(pw, dfdcWmag, dfdcr; label="DFDC W")
    savefig(pw, project_dir * "/examples/dfdc_comp/Wmagcomp.pdf")

    plot!(por, omega * rpc, rpc; label="DuctTAPE")
    plot!(por, omega * dfdcr, dfdcr; linewidth=2, linestyle=:dash, label="DFDC")
    savefig(por, project_dir * "/examples/dfdc_comp/Omegarcomp.pdf")

    return nothing
end

function plotangles(cdump, inputs)
    # - Extract DFDC Stuff - #
    dfdctwist = dfdc_printout[:, 3]
    dfdcalpha = dfdc_printout[:, 4]

    # - extract ducttape radial positions - #
    rpc = inputs.rotor_panel_centers
    twist = inputs.blade_elements[1].twists

    # - Plot Angles - #
    pa = plot(; xlabel="Angles (deg)", ylabel="radial positions")

    plot!(pa, twist * 180.0 / pi, rpc; color=mycolors[1], label="DuctTAPE blade twist")
    plot!(
        pa, dfdctwist, dfdcr; color=mycolors[1], linestyle=:dash, label="DFDC blade twist"
    )

    plot!(pa, cdump.phi * 180.0 / pi, rpc; color=mycolors[2], label="DuctTAPE inflow angle")
    plot!(pa, dfdcphi, dfdcr; color=mycolors[2], linestyle=:dash, label="DFDC inflow angle")

    plot!(
        pa, cdump.alpha * 180.0 / pi, rpc; color=mycolors[3], label="DuctTAPE attack angle"
    )
    plot!(
        pa, dfdcalpha, dfdcr; color=mycolors[3], linestyle=:dash, label="DFDC attack angle"
    )

    savefig(pa, project_dir * "/examples/dfdc_comp/anglescomp.pdf")

    return nothing
end

function plotclcd(cdump, inputs)

    # - extract ducttape radial positions - #
    rpc = inputs.rotor_panel_centers
    cl = cdump.cl
    cd = cdump.cd

    # - Plot Coefficients - #
    pcl = plot(; xlabel=L"c_\ell", ylabel="radial positions")

    plot!(pcl, cl, rpc; color=mycolors[1], label="DuctTAPE")
    plot!(pcl, dfdccl, dfdcr; color=mycolors[1], linestyle=:dash, label="DFDC")

    pcd = plot(; xlabel=L"c_d", ylabel="radial positions")

    plot!(pcd, cd, rpc; color=mycolors[1], label="DuctTAPE")
    plot!(pcd, dfdccd, dfdcr; color=mycolors[1], linestyle=:dash, label="DFDC")

    savefig(pcl, project_dir * "/examples/dfdc_comp/clcomp.pdf")
    savefig(pcd, project_dir * "/examples/dfdc_comp/cdcomp.pdf")

    return nothing
end

function ploths(cdump, inputs)

    # - extract ducttape radial positions - #
    rpc = inputs.rotor_panel_centers

    # - Calculate enthalpy disk jump - #
    Htilde = dt.calculate_enthalpy_jumps(
        cdump.Gamr, inputs.blade_elements[1].Omega, inputs.blade_elements[1].B
    )

    # - Calculate entropy disk jump - #
    Stilde = dt.calculate_entropy_jumps(cdump.sigr, cdump.Wm_rotor)

    ph = plot(; xlabel=L"\widetilde{H}", ylabel="radial positions")
    ps = plot(; xlabel=L"\widetilde{S}", ylabel="radial positions")

    plot!(ph, Htilde, rpc; label="DuctTAPE")
    plot!(ph, deltaH, dfdcr; label="DFDC")

    plot!(ps, Stilde, rpc; label="DuctTAPE")
    plot!(ps, deltaS, dfdcr; label="DFDC")

    savefig(ph, project_dir * "/examples/dfdc_comp/Htildecomp.pdf")
    savefig(ps, project_dir * "/examples/dfdc_comp/Stildecomp.pdf")

    return nothing
end

dump = runandplot()
