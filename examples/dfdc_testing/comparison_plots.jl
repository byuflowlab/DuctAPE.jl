project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/examples/dfdc_testing/"

include(project_dir * "/plots_default.jl")

using DuctTAPE
const dt = DuctTAPE

using FLOWMath
const fm = FLOWMath

function runandplot(; filename="straightwake.case.jl")
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

    return out
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
    plot!(
        pgeom,
        inputs.blade_elements[1].xrotor * ones(length(inputs.rotor_panel_edges)),
        inputs.rotor_panel_edges;
        color=mycolors[2],
        linewidth=0.25,
        markersize=0.5,
        markershape=:rect,
        label="Rotor Panel Edges",
    )

    plot!(
        pgeom,
        inputs.rotor_source_panels[1].panel_center[:, 1],
        inputs.rotor_source_panels[1].panel_center[:, 2];
        color=mycolors[3],
        seriestype=:scatter,
        markersize=0.75,
        markershape=:circle,
        label="Rotor Panel Centers",
    )

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
    # - Get DFDC Data - #
    include(datapath * "DFDC_HUB_FORCES.jl")
    dfdc_hub_cp = dfdc_hub_thrust[:, 3]
    include(datapath * "DFDC_DUCT_FORCES.jl")
    dfdc_duct_cp = dfdc_duct_thrust[:, 3]
    include(datapath * "DFDC_ELEMENT_CENTERS.jl")
    dfdc_hubx = elem1[:, 2]
    dfdc_hubr = elem1[:, 3]
    dfdc_ductx = elem2[:, 2]
    dfdc_ductr = elem2[:, 3]

    if convergeflag
        convlabel = "Converged"
    else
        convlabel = "NOT converged"
    end

    #---------------------------------#
    #    Extract Outputs and Plot     #
    #---------------------------------#
    (; gamb, gamw, Gamr, sigr) = out

    ##### ----- Plot duct surface velocity ----- #####

    # # initialize plot
    # pb = plot(; xlabel="x", ylabel="Q/Qref")

    # plot!(
    #     pb,
    #     dp,
    #     abs.(gamb[1:length(dp)]) ./ Vref;
    #     color=mycolors[1],
    #     label=convlabel * " DuctTAPE Duct",
    # )

    #plot!(
    #    pb,
    #    dfdc_ductx,
    #    #TODO: get surface velocities
    #    color=mycolors[1],
    #    linestyle=:dash,
    #    label="DFDC Duct",
    #)

    # plot!(
    #     pb,
    #     hp,
    #     abs.(gamb[(length(dp) + 1):end]) ./ Vref;
    #     color=mycolors[2],
    #     label=convlabel * " DuctTAPE Hub",
    # )

    #plot!(
    #    pb,
    #    dfdc_hubx,
    #   #TODO: get surface velocities
    #    color=mycolors[2],
    #    linestyle=:dash,
    #    label="DFDC Hub",
    #)

    ##plot rotor location
    #plot!(
    #    pb,
    #    inputs.blade_elements[1].xrotor * ones(2),
    #    abs.([0.0; maximum(abs.(gamb ./ Vref))]);
    #    # linewidth=0.25,
    #    linestyle=:dot,
    #    color=mycolors[3],
    #    label="rotor location",
    #)

    ##save
    #savefig(pb, datapath * "body-velocity.pdf")

    ##### ----- Plot Surface Pressure ----- #####

    pcp = plot(; xlabel="x", ylabel=L"C_p", yflip=true)

    plot!(
        pcp,
        [out.duct_inner_x; out.duct_outer_x],
        [out.duct_inner_cp; out.duct_outer_cp];
        color=mycolors[1],
        label="DuctTAPE Duct surface",
    )

    plot!(
        pcp, dfdc_ductx, dfdc_duct_cp; color=mycolors[1], linestyle=:dash, label="DFDC Duct"
    )

    plot!(pcp, out.hub_x, out.hub_cp; color=mycolors[2], label="DuctTAPE Hub")
    plot!(pcp, dfdc_hubx, dfdc_hub_cp; color=mycolors[2], linestyle=:dash, label="DFDC Hub")

    #save
    savefig(pcp, datapath * "body-pressure.pdf")

    return nothing
end

out = runandplot(; filename="straightwake.case.jl")