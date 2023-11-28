project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/dev_debug_archive/dfdc_testing/"

include(project_dir * "/plots_default.jl")

using DuctAPE
const dt = DuctAPE

using FLOWMath
const fm = FLOWMath

function runandplot(; filename="straightwake.case.jl")
    println("Setting Up")
    include(datapath * filename)

    nb = 20:5:40
    # nb = 1
    # nw = 1:20
    nw = 1
    effs = zeros(length(nb),length(nw))
    rotor_thrusts = zeros(length(nb),length(nw))
    body_thrusts = zeros(length(nb),length(nw))
    total_thrusts = zeros(length(nb),length(nw))
    num_body_panels_conv = zeros(length(nb),length(nw))
    num_wake_x_panels_conv = zeros(length(nb),length(nw))
    nwake_sheets_conv = zeros(length(nb),length(nw))

    brange = nb == 1 ? (4:4) : (nb)
    wrange = nw == 1 ? (4:4) : (nw)
    idb = 1
    for ib in brange
        idw = 1
        for iw in wrange
            nhub_inlet = 10 * ib
            nduct_inlet = 10 * ib
            nwake_sheets = 3 * iw
            npanels = [10; 20] .* ib
            ref_paneling_constants = (;
                npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets
            )

            println("Running Analysis: $idb, $idw of $nb, $nw")
            out, converged_states, inputs, initial_states, convergeflag = rundt(
                duct_coordinates,
                hub_coordinates,
                rotor_parameters,
                ref_paneling_constants,
                freestream,
                reference_parameters,
            )

            println("Saving Outputs: $idb, $idw of $nb, $nw")
            num_body_panels_conv[idb, idw] = inputs.num_body_panels
            num_wake_x_panels_conv[idb, idw] = inputs.num_wake_x_panels
            nwake_sheets_conv[idb, idw] = nwake_sheets
            effs[idb, idw] = out.total_efficiency
            total_thrusts[idb, idw] = out.total_thrust
            body_thrusts[idb, idw] = out.body_thrust
            rotor_thrusts[idb, idw] = out.rotor_total_thrust[1]

            # println("Plotting Geometry: $idb, $idw of $nb, $nw")
            # plotgeom(
            #     inputs,
            #     ref_paneling_constants;
            #     filename="geometry_nbpan-$(inputs.num_body_panels)_nwr-$(nwake_sheets)_nwx-$(inputs.num_wake_x_panels).pdf",
            # )

            idw += 1
        end
        idb += 1
    end

    println("Plotting Convergence Stuff")
    plotconv(
        num_body_panels_conv,
        nwake_sheets_conv,
        num_wake_x_panels_conv,
        effs,
        (; total_thrusts, body_thrusts, rotor_thrusts);
        # fileprefix="nb$nb-nw$nw-",
        fileprefix="more-body-panels-",
    )

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

function plotgeom(inputs, paneling_constants; filename="precomputed_geometry.pdf")

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

    savefig(pgeom, datapath * filename)

    return nothing
end

function plotconv(
    num_body_panels_conv,
    num_wake_x_panels_conv,
    nwake_sheets_conv,
    effs,
    thrusts;
    fileprefix="",
)
    ntot = num_body_panels_conv .+ num_wake_x_panels_conv .+ nwake_sheets_conv
    pe = plot(;
        xlabel="total number of panels",
        ylabel="efficiency",
        title="body and wake x panel convergence",
    )
    pre = plot(;
        xlabel="number of wake sheets",
        ylabel="efficiency",
        title="number of wakes convergence",
    )
    pt = plot(;
        xlabel="total number of panels",
        ylabel="thrust",
        title="body and wake x panel convergence",
    )
    prt = plot(;
        xlabel="number of wake sheets", ylabel="thrust", title="number of wakes convergence"
    )

    plot!(pe, ntot[:, 1], effs[:, 1]; marker=true)
    savefig(pe, datapath * fileprefix * "bodywake-efficiency-convergence.pdf")

    plot!(pt, ntot[:, 1], thrusts.total_thrusts[:, 1]; marker=true, label="Total")
    plot!(pt, ntot[:, 1], thrusts.rotor_thrusts[:, 1]; marker=true, label="Rotor(s)")
    plot!(pt, ntot[:, 1], thrusts.body_thrusts[:, 1]; marker=true, label="Bodies")
    savefig(pt, datapath * fileprefix * "bodywake-thrust-convergence.pdf")

    plot!(pre, nwake_sheets_conv[1, :], effs[1, :]; marker=true)
    savefig(pre, datapath * fileprefix * "nwake-efficiency-convergence.pdf")

    plot!(
        prt,
        nwake_sheets_conv[1, :],
        thrusts.total_thrusts[1, :];
        marker=true,
        label="Total",
    )
    plot!(
        prt,
        nwake_sheets_conv[1, :],
        thrusts.rotor_thrusts[1, :];
        marker=true,
        label="Rotor(s)",
    )
    plot!(
        prt,
        nwake_sheets_conv[1, :],
        thrusts.body_thrusts[1, :];
        marker=true,
        label="Bodies",
    )
    savefig(prt, datapath * fileprefix * "nwake-thrust-convergence.pdf")

    return nothing
end
-
runandplot()
