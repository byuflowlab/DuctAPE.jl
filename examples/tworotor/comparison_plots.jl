project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/examples/tworotor/"

include(project_dir * "/visualize/plots_default.jl")

using DuctTAPE
const dt = DuctTAPE

using FLOWMath
const fm = FLOWMath

function runandplot(; filename="test.case.jl")
    println("Setting Up")
    include(datapath * filename)

    println("Running Analysis")
    out, converged_states, inputs, initial_states, convergeflag = rundt(
        duct_coordinates,
        hub_coordinates,
        rotor_parameters,
        paneling_constants,
        freestream,
        reference_parameters,
    )

    println("Plotting Geometry")
    plotgeom(inputs, paneling_constants)

    # println("Plotting Rotor/Wake States")
    # plotstates(cdump.Gamr, cdump.sigr, cdump.gamw, inputs, convergeflag)

    # println("Plotting Blade Velocities")
    # plotbladevelocities(cdump, inputs)

    # println("Plotting Blade Angles")
    # plotangles(cdump, inputs)

    # println("Plotting Lift/Drag on Blade")
    # plotclcd(cdump, inputs)

    # println("Plotting Enthalpy and Entropy Disk Jumps")
    # ploths(cdump, inputs)

    println("Plotting Body Aerodynamics")
    plotbodyaero(out, convergeflag, Vref)

    return out, inputs
end

function rundt(
    duct_coordinates,
    hub_coordinates,
    rotor_parameters,
    paneling_constants,
    freestream,
    reference_parameters,
)
    #---------------------------------#
    #           Run Solver            #
    #---------------------------------#

    out, converged_states, inputs, initial_states, convergeflag = dt.analyze_propulsor(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters;
        debug=false,
        verbose=true,
        maximum_linesearch_step_size=1e6,
        iteration_limit=100,
    )

    return out, converged_states, inputs, initial_states, convergeflag
end

function plotgeom(inputs, paneling_constants)

    ##### ----- Plot GEOMETRY ----- #####
    #initialize plot
    pgeom = plot(; aspectratio=1, xlabel="x", ylabel="r")
    nrotor = length(inputs.blade_elements.xrotor)
    for ir in 1:nrotor
        plot!(
            pgeom,
            inputs.blade_elements[ir].xrotor *
            ones(length(inputs.rotor_panel_edges[:, ir])),
            inputs.rotor_panel_edges[:, ir];
            color=mycolors[2],
            linewidth=0.25,
            markersize=0.5,
            markershape=:rect,
            label="Rotor Panel Edges",
        )

        plot!(
            pgeom,
            inputs.rotor_source_panels[ir].panel_center[:, 1],
            inputs.rotor_source_panels[ir].panel_center[:, 2];
            color=mycolors[3],
            seriestype=:scatter,
            markersize=0.75,
            markershape=:circle,
            label="Rotor Panel Centers",
        )
    end

    for iw in 1:(paneling_constants.nwake_sheets)
        if iw == 1
            wakelab = "Wake Panel Edges"
            wakelab2 = "Wake Panel Centers"
        else
            wakelab = ""
            wakelab2 = ""
        end
        plot!(
            pgeom,
            inputs.wakexgrid[:, iw],
            inputs.wakergrid[:, iw];
            linewidth=0.25,
            markersize=0.5,
            markershape=:rect,
            color=:black,
            label=wakelab,
        )

        plot!(
            pgeom,
            inputs.wake_vortex_panels[iw].panel_center[:, 1],
            inputs.wake_vortex_panels[iw].panel_center[:, 2];
            seriestype=:scatter,
            markersize=0.75,
            markershape=:circle,
            color=mycolors[2],
            label=wakelab2,
        )
    end

    savefig(datapath * "precomputed-rotor-wake-geometry.pdf")

    # plot body panels
    plot!(
        pgeom,
        inputs.duct_coordinates[:, 1],
        inputs.duct_coordinates[:, 2];
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        color=mycolors[3],
        label="Body Panel Edges",
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
        if ib == 1
            bodlab = "Body Panel Centers"
        else
            bodlab = ""
        end
        plot!(
            pgeom,
            inputs.body_panels[ib].panel_center[:, 1],
            inputs.body_panels[ib].panel_center[:, 2];
            color=mycolors[1],
            seriestype=:scatter,
            markersize=0.75,
            label=bodlab,
        )
    end

    savefig(pgeom, datapath * "precomputed-full-geometry.pdf")

    return nothing
end

function plotstates(Gamr, sigr, gamw, inputs, convergeflag)

    #TODO: need sigr and gamw from dfdc

    if convergeflag
        convlabel = "DuctTAPE Converged"
    else
        convlabel = "NOT converged"
    end

    ##### ----- Plot rotor circulation distribution ----- #####
    # initialize plot
    pG = plot(; xlabel=L"\Gamma", ylabel="r")

    # plot solution
    plot!(pG, Gamr, inputs.rotor_panel_centers; label=convlabel)
    plot!(pG, dfdcGamr, dfdcr; label="DFDC")

    #save
    savefig(pG, datapath * "rotorcirculation-check.pdf")

    ##### ----- Plot Wake Strengths ----- #####

    #pg = plot(; xlabel=L"\gamma_\theta^{wake}", ylabel="r")

    ## plot solution
    #plot!(pg, gamw, inputs.rotor_panel_edges; label=convlabel)

    ##save
    #savefig(pg, datapath* "wake-strength.pdf")

    ##### ----- Plot Source Strengths ----- #####
    ps = plot(; xlabel=L"\sigma", ylabel="r")

    # plot solution
    plot!(ps, sigr, inputs.rotor_panel_centers; label=convlabel)
    plot!(ps, dfdcsigma, dfdcr; label="DFDC")

    #save
    savefig(ps, datapath * "source-strength-check.pdf")

    return nothing
end

function plotbodyaero(out, convergeflag, Vref)

    # # - Get DFDC Data - #
    # include(datapath * "DFDC_HUB_FORCES.jl")
    # dfdc_hub_cp = dfdc_hub_thrust[:, 3]
    # include(datapath * "DFDC_DUCT_FORCES.jl")
    # dfdc_duct_cp = dfdc_duct_thrust[:, 3]
    # include(datapath * "DFDC_DUCTWAKE_FORCES.jl")
    # dfdc_ductwake_cp = dfdc_ductwake_thrust[:, 3]
    # include(datapath * "DFDC_HUBWAKE_FORCES.jl")
    # dfdc_hubwake_cp = dfdc_hubwake_thrust[:, 2]
    # include(datapath * "DFDC_ELEMENT_CENTERS.jl")
    # dfdc_hubx = elem1[:, 2]
    # dfdc_hubr = elem1[:, 3]
    # dfdc_ductx = elem2[:, 2]
    # dfdc_ductr = elem2[:, 3]
    # dfdc_ductwakex = elem14[:, 2]
    # dfdc_ductwaker = elem14[:, 3]
    # dfdc_hubwakex = elem4[:, 2]
    # dfdc_hubwaker = elem4[:, 3]
    # include(datapath * "DFDC_SURFACE_VELOCITIES.jl")
    # dfdc_duct_vs = dfdc_vs[Int.(elem2[:, 1]), 2]
    # dfdc_ductwake_vs = dfdc_vs[Int.(elem14[:, 1]), 2]
    # dfdc_hub_vs = dfdc_vs[Int.(elem1[:, 1]), 2]
    # dfdc_hubwake_vs = dfdc_vs[Int.(elem4[:, 1]), 2]

    #---------------------------------#
    #              Plot               #
    #---------------------------------#
    # (; gamb, gamw, Gamr, sigr) = out

    pvs = plot(; xlabel="x", ylabel=L"V_s")

    # duct surface
    plot!(
        pvs,
        [out.duct_inner_x; out.duct_outer_x],
        abs.([out.duct_inner_vs; out.duct_outer_vs]);
        color=mycolors[1],
        label="DuctTAPE Duct surface",
    )

    # # dfdc duct surface
    # plot!(
    #     pvs, dfdc_ductx, dfdc_duct_vs; color=mycolors[1], linestyle=:dash, label="DFDC Duct"
    # )

    # hub surface
    plot!(pvs, out.hub_x, abs.(out.hub_vs); color=mycolors[2], label="DuctTAPE Hub")

    # # dfdc hub surface
    # plot!(pvs, dfdc_hubx, dfdc_hub_vs; color=mycolors[2], linestyle=:dash, label="DFDC Hub")

    # # duct wake
    # plot!(pvs, out.ductwake_x, abs.(out.ductwake_vs); color=mycolors[3], label="DuctTAPE Duct Wake")

    # # dfdc duct wake
    # plot!(
    #     pvs,
    #     dfdc_ductwakex,
    #     dfdc_ductwake_vs;
    #     color=mycolors[3],
    #     linestyle=:dash,
    #     label="DFDC Duct Wake",
    # )

    # hub wake
    plot!(
        pvs,
        out.hubwake_x,
        abs.(out.hubwake_vs);
        color=mycolors[4],
        label="DuctTAPE Hub Wake",
    )

    # # dfdc hub wake
    # plot!(
    #     pvs,
    #     dfdc_hubwakex,
    #     dfdc_hubwake_vs;
    #     color=mycolors[4],
    #     linestyle=:dash,
    #     label="DFDC Hub Wake",
    # )

    #save
    savefig(pvs, datapath * "body-velocity.pdf")

    ##### ----- Plot Surface Pressure ----- #####

    pcp = plot(; xlabel="x", ylabel=L"C_p", yflip=true)

    #duct surface
    plot!(
        pcp,
        [out.duct_inner_x; out.duct_outer_x],
        [out.duct_inner_cp; out.duct_outer_cp];
        color=mycolors[1],
        label="DuctTAPE Duct surface",
    )

    ##dfdc duct surface
    #plot!(
    #    pcp, dfdc_ductx, dfdc_duct_cp; color=mycolors[1], linestyle=:dash, label="DFDC Duct"
    #)

    #hub surface
    plot!(pcp, out.hub_x, out.hub_cp; color=mycolors[2], label="DuctTAPE Hub")

    ##dfdc hub surface
    #plot!(pcp, dfdc_hubx, dfdc_hub_cp; color=mycolors[2], linestyle=:dash, label="DFDC Hub")

    #duct wake
    plot!(
        pcp, out.ductwake_x, out.ductwake_cp; color=mycolors[3], label="DuctTAPE Duct Wake"
    )

    # # dfdc duct wake
    # plot!(
    #     pcp,
    #     dfdc_ductwakex,
    #     dfdc_ductwake_cp;
    #     color=mycolors[3],
    #     linestyle=:dash,
    #     label="DFDC Duct Wake",
    # )

    #hub wake
    plot!(pcp, out.hubwake_x, out.hubwake_cp; color=mycolors[4], label="DuctTAPE Hub Wake")

    ##dfdc hub wake
    #plot!(
    #    pcp,
    #    dfdc_hubwakex,
    #    dfdc_hubwake_cp;
    #    color=mycolors[4],
    #    linestyle=:dash,
    #    label="DFDC Hub Wake",
    #)

    #save
    savefig(pcp, datapath * "body-pressure.pdf")

    return nothing
end

function testpieces(; filename="tworotorcase.jl")
    println("Setting Up")
    include(datapath * filename)

    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters;
        debug=false,
    )

    # calculate initial guess for state variables
    initial_states = dt.initialize_states(inputs)
    initials = copy(initial_states)

    p = (; debug=false, verbose=false)

    # - Define closure that allows for parameters - #
    rwrap(r, states) = dt.residual!(r, states, inputs, p)

    res = NLsolve.nlsolve(
        rwrap,
        initial_states;
        autodiff=:forward,
        method=:newton,
        iterations=25,
        show_trace=true,
        linesearch=BackTracking(; maxstep=1e6),
    )

    return inputs, res.zero
end

# out, inputs = runandplot(; filename="test.case.jl")
# out, inputs = runandplot(; filename="tworotorcase.jl")

# inputs, states = testpieces(; filename="tworotorcase.jl")
