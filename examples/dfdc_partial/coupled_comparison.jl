project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE
using PrettyTables
const pt = PrettyTables

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath * "outputs/"

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/visualize_flowfield.jl")
include(project_dir * "/visualize/plots_default_new.jl")
# pyplot()

# read in parameters file
include(datapath * "lewis_with_rotor.case.jl")

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 50.0))
pcd = plot(; xlabel="x", ylabel=L"c_p", ylim=(-2.25, 1.0), yflip=true)
pch = plot(; xlabel="x", ylabel=L"c_p", ylim=(-2.25, 1.0), yflip=true)

#---------------------------------#
#               DFDC              #
#---------------------------------#
# include(datapath * "DFDC_VELOCITY_BREAKDOWN.jl")
# include(datapath * "DFDC_SURFACE_VELOCITIES.jl")
# include(datapath * "DFDC_ELEMENT_IDS.jl")
include(datapath * "dfdc_outs/lewis_coupled/DFDC_CPS.jl")
include(datapath * "dfdc_outs/lewis_coupled/DFDC_CIRC.jl")

hubx = dfdc_hub_cp[:, 1]
hubr = dfdc_hub_cp[:, 2]
hubvs = dfdc_hub_cp[:, end - 1]
hubcp = dfdc_hub_cp[:, 4]

ductx = dfdc_duct_cp[:, 1]
ductr = dfdc_duct_cp[:, 2]
ductvs = dfdc_duct_cp[:, end - 1]
ductcp = dfdc_duct_cp[:, 4]

ductwakex = dfdc_wake14_cp[:, 1]
ductwaker = dfdc_wake14_cp[:, 2]
ductwakevs = dfdc_wake14_cp[:, end - 1]
ductwakecp = dfdc_wake14_cp[:, 4]

plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[1], label="DFDC Hub")
plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[1], label="DFDC Duct")
# plot!(pv, hubwakex, hubwakevs; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pv, ductwakex, ductwakevs; linestyle=:dash, label="DFDC Duct Wake")

plot!(pch, hubx, hubcp; linestyle=:dash, color=myred[1], label="DFDC Hub")
plot!(pcd, ductx, ductcp; linestyle=:dash, color=myblue[1], label="DFDC")
# plot!(pcd, hubwakex, hubwakecp; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pcd, ductwakex, ductwakecp; linestyle=:dash, label="DFDC Duct Wake")

pr = plot(; xlabel="Circulation", ylable="r")

dfdc_duct_coordinates = reverse([ductx ductr]; dims=1)
dfdc_hub_coordinates = reverse([hubx hubr]; dims=1)
dfdc_panels = dt.generate_panels([dfdc_duct_coordinates, dfdc_hub_coordinates])

##viz dfdc geometry
#visualize_paneling(;
#    body_panels=dfdc_panels,
#    coordinates=[dfdc_duct_coordinates, dfdc_hub_coordinates],
#    controlpoints=true,
#    nodes=true,
#    normals=true,
#    normal_scaling=0.1,
#    savepath=savepath,
#    filename=["lc-dfdcbodygeometry.pdf"],
#    legendloc=:right,
#)

#---------------------------------#
#             Exp data            #
#---------------------------------#
include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")

# plot!(
#     pcd,
#     pressurexupper,
#     pressureupper;
#     seriestype=:scatter,
#     color=myblue[1],
#     markershape=:utriangle,
#     markersize=3,
#     label="Exp Outer Duct",
# )

# plot!(
#     pcd,
#     pressurexlower,
#     pressurelower;
#     seriestype=:scatter,
#     color=myblue[1],
#     markershape=:dtriangle,
#     markersize=3,
#     label="Exp Inner Duct",
# )

# plot!(
#     pv,
#     Vs_over_Vinf_x,
#     Vs_over_Vinf_vs * Vinf;
#     seriestype=:scatter,
#     color=myred[1],
#     markershape=:utriangle,
#     markersize=3,
#     label="Exp Hub",
# )

###########################################
# savefig(pv, savepath * "lewis-with-rotor-velocity-dfdc.pdf")
# savefig(pcd, savepath * "lewis-with-rotor-pressure-dfdc.pdf")

#---------------------------------#
#             DuctAPE            #
#---------------------------------#
# include(datapath * "lewis_refined_duct.jl")
# include(datapath * "lewis_refined_hub.jl")

# - Generate Panels - #

#overwrite geometry with DFDC input geometry
duct_coordinates = reverse([ductx ductr]; dims=1)
duct_coordinates[1, :] .= duct_coordinates[end, :]
hub_coordinates = reverse([hubx hubr]; dims=1)

# - Use Smooth Duct Coordinates - #
include(datapath * "nacasmoothgeom.jl")
revductcoords = reverse(smoothnormduct; dims=1)
revductcoords[:, 2] .+= duct_coordinates[1, 2]
duct_coordinates = dt.repanel_airfoil(revductcoords; normalize=false, N=200)

#overwite npanels for now (should put conditions for this in auto generation function)

refine = [4; 8; 16; 32; 64; 128]
# refine = [4; 8; 16]

total_thrust_err = zeros(length(refine))
body_thrust_err = zeros(length(refine))
rotor_thrust_err = zeros(length(refine))
total_efficiency_err = zeros(length(refine))
rotor_efficiency_err = zeros(length(refine))
power_err = zeros(length(refine))
torque_err = zeros(length(refine))
bpan = zeros(length(refine))

for (ip, pref) in enumerate(refine)
    global paneling_constants

    paneling_constants = (;
        paneling_constants...,
        nhub_inlet=ceil(Int, pref * 4),
        nduct_inlet=ceil(Int, pref * 4),
        npanels=ceil.(Int, pref .* [3, 2, 2]),
    )

    # # run analyze_propulsor function
    # out, converged_states, inputs, initial_states, convergeflag = @time dt.analyze_propulsor(
    #     duct_coordinates,
    #     hub_coordinates,
    #     paneling_constants,
    #     rotor_parameters,
    #     freestream,
    #     reference_parameters;
    #     debug=false,
    #     verbose=true,
    #     maximum_linesearch_step_size=1e6,
    #     iteration_limit=50,
    # )

    ## -- do things one at a time, step by step -- ##

    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters;
    )

    ptot = inputs.body_doublet_panels.totpanel
    bpan[ip] = ptot
    println("\nTotal body panels = ", ptot, "\n\n")

    # calculate initial guess for state variables
    initial_states = dt.initialize_states(inputs)

    # - Define closure that allows for parameters - #
    p = (;)
    rwrap(r, statesvars) = dt.residual!(r, statesvars, inputs, p)

    @time "Solve" begin
        res = dt.NLsolve.nlsolve(
            rwrap,
            initial_states;
            autodiff=:forward,
            method=:newton,
            iterations=50,
            show_trace=true,
            linesearch=dt.BackTracking(; maxstep=1e6),
        )
    end

    states = res.zero

    out = @time "Post-process" dt.post_process(states, inputs)

    #---------------------------------#
    #           Comparisons           #
    #---------------------------------#
    println("\n\nGENERATING COMPARISON AND VISUALZITION\n")

    ## -- Print Post Processed Comparsion Values -- ##
    include(datapath * "dfdc_outs/lewis_coupled/DFDC_FORCES.jl")

    function perr(x, y)
        return 100.0 * (y - x) / x
    end

    # save errors
    total_thrust_err[ip] = perr(dfdc_total_thrust, out.total_thrust)
    total_efficiency_err[ip] = perr(dfdc_total_efficiency, out.total_efficiency)
    body_thrust_err[ip] = perr(dfdc_body_thrust, out.body_thrust)
    rotor_thrust_err[ip] = perr(dfdc_rotor_thrust[1], out.rotor_thrust[1])
    torque_err[ip] = perr(dfdc_rotor_torque[1], out.rotor_torque[1])
    power_err[ip] = perr(dfdc_rotor_power[1], out.rotor_power[1])
    rotor_efficiency_err[ip] = perr(dfdc_rotor_efficiency[1], out.rotor_efficiency[1])

    # things to print:
    outcomps = [
        "Total Thrust (N):" dfdc_total_thrust out.total_thrust perr(dfdc_total_thrust, out.total_thrust)
        "Total Torque (N-m):" sum(dfdc_rotor_torque) out.total_torque perr(sum(dfdc_rotor_torque), out.total_torque)
        "Total Power (W):" sum(dfdc_rotor_power) out.total_power perr(sum(dfdc_rotor_power), out.total_power)
        "Total Efficiency:" dfdc_total_efficiency out.total_efficiency perr(dfdc_total_efficiency, out.total_efficiency)
        "Body Thrust (N):" dfdc_body_thrust out.body_thrust perr(dfdc_body_thrust, out.body_thrust)
    ]

    outscomprot = [
        [
            "Rotor $ir Thrust (N):" dfdc_rotor_thrust[ir] out.rotor_thrust[ir] perr(dfdc_rotor_thrust[ir], out.rotor_thrust[ir])
            "Rotor $ir Torque (N-m):" dfdc_rotor_torque[ir] out.rotor_torque[ir] perr(dfdc_rotor_torque[ir], out.rotor_torque[ir])
            "Rotor $ir Power (W):" dfdc_rotor_power[ir] out.rotor_power[ir] perr(dfdc_rotor_power[ir], out.rotor_power[ir])
            "Rotor $ir Efficiency:" dfdc_rotor_efficiency[ir] out.rotor_efficiency[ir] perr(dfdc_rotor_efficiency[ir], out.rotor_efficiency[ir])
        ] for ir in 1:(inputs.num_rotors)
    ]

    # print the table
    pt.pretty_table(
        [outcomps; reduce(vcat, outscomprot)];
        header=["Value", "DFDC", "DuctAPE", "% Error"],
        backend=Val(:latex),
    )

    ## -- Plotting -- ##

    @time "Visualize 2D" begin
        println("Visualizing Paneling")
        # check normals

        visualize_paneling(;
            body_panels=inputs.body_doublet_panels,
            rotor_panels=inputs.rotor_source_panels,
            coordinates=[duct_coordinates, hub_coordinates],
            controlpoints=true,
            nodes=true,
            normals=true,
            normal_scaling=0.1,
            savepath=savepath,
            filename=[
                "$(ptot)_panels_lc-bodygeometry.pdf"
                "$(ptot)_panels_lc-bodygeometry.png"
            ],
            legendloc=:right,
        )

        # # everything
        # visualize_paneling(;
        #     body_panels=inputs.body_doublet_panels,
        #     rotor_panels=inputs.rotor_source_panels,
        #     wake_panels=inputs.wake_vortex_panels,
        #     coordinates=[duct_coordinates, hub_coordinates],
        #     controlpoints=true,
        #     nodes=false,
        #     normals=false,
        #     normal_scaling=0.1,
        #     savepath=savepath,
        #     filename=[
        #         "$(ptot)_panels_lc-fullgeometry.pdf"
        #         "$(ptot)_panels_lc-fullgeometry.png"
        #     ],
        #     legendloc=:outerright,
        # )

        # # everything
        # visualize_paneling(;
        #     body_panels=inputs.body_doublet_panels,
        #     rotor_panels=inputs.rotor_source_panels,
        #     wake_panels=inputs.wake_vortex_panels,
        #     wakeinterfaceid=inputs.ductwakeinterfaceid,
        #     coordinates=[duct_coordinates, hub_coordinates],
        #     controlpoints=true,
        #     nodes=true,
        #     TEnodes=false,
        #     normals=false,
        #     normal_scaling=0.1,
        #     savepath=savepath,
        #     filename=["$(ptot)_panels_lc-fullgeometry-zoom.pdf"; "$(ptot)_panels_lc-fullgeometry-zoom.png"],
        #     legendloc=:outerright,
        #     limits=(; ylim=(0.65, 0.8), xlim=(0.4, 0.6)),
        #     zoom=true,
        # )

    end

    # println("Generating VTKs of Velocity Field")
    # # Flowfield VTK generation
    # visualize_flowfield(
    #     inputs.Vinf;
    #     body_panels=inputs.body_doublet_panels,
    #     rotor_panels=inputs.rotor_source_panels,
    #     wake_panels=inputs.wake_vortex_panels,
    #     mub=out.mub,
    #     sigr=out.sigr,
    #     gamw=out.gamw,
    #     Gamr=out.Gamr,
    #     # Pmax=nothing,
    #     # Pmin=nothing,
    #     verbose=true,
    #     run_name="$(ptot)_panels_lc-velocity_fields",
    #     save_path=savepath,
    #     cellsizescale=0.005,
    # )

    # visualize_surfaces(;
    #     body_panels=inputs.body_doublet_panels,
    #     rotor_panels=inputs.rotor_source_panels,
    #     wake_panels=inputs.wake_vortex_panels,
    #     run_name="$(ptot)_panels_lc-surfaces",
    #     save_path=savepath,
    # )

    # - Plot surface Velocity and Pressure - #
    plot!(
        pv,
        out.duct_inner_x,
        out.duct_inner_vs;
        color=myblue[2],
        label="DuctAPE Duct inner",
    )
    plot!(
        pv,
        out.duct_outer_x,
        out.duct_outer_vs;
        color=myblue[3],
        label="DuctAPE Duct outer",
    )
    plot!(pv, out.hub_x, out.hub_vs; color=myred[2], label="DuctAPE Hub")
    plot!(
        pcd,
        out.duct_inner_x,
        out.duct_inner_cp;
        color=ip,
        label="DuctAPE $(ptot) Body Panels",
    )
    plot!(pcd, out.duct_outer_x, out.duct_outer_cp; color=ip, label="")
    plot!(pch, out.hub_x, out.hub_cp; color=ip, label="DuctAPE $(ptot) Body Panels")

    #- plot rotor circulation - #
    for ir in inputs.num_rotors
        if ip == 1
            plot!(
                pr,
                dfdc_circulation / inputs.blade_elements[ir].B,
                inputs.blade_elements[ir].rbe;
                linestyle=:dash,
                color=myblue[1],
                label="DFDC",
            )
        end
        plot!(
            pr,
            out.Gamr[:, ir],
            inputs.blade_elements[ir].rbe;
            color=ip,
            label="DuctAPE $(ptot) Body Panels",
        )
    end

    # # - plot body surface velocity components - #
    # duct_inner_vxr_nogradmu = repeat([inputs.Vinf 0], length(out.duct_inner_x))
    # duct_inner_vxr_nogradmu .+= out.duct_inner_vs_from_body .+ out.duct_inner_vs_from_TE
    # duct_inner_vs_nogradmu = dt.norm.(eachrow(duct_inner_vxr_nogradmu))
    # duct_outer_vxr_nogradmu = repeat([inputs.Vinf 0], length(out.duct_outer_x))
    # duct_outer_vxr_nogradmu .+= out.duct_outer_vs_from_body .+ out.duct_outer_vs_from_TE
    # duct_outer_vs_nogradmu = dt.norm.(eachrow(duct_outer_vxr_nogradmu))
    # hub_vxr_nogradmu = repeat([inputs.Vinf 0], length(out.hub_x))
    # hub_vxr_nogradmu .+= out.hub_vs_from_body .+ out.hub_vs_from_TE
    # hub_vs_nogradmu = dt.norm.(eachrow(hub_vxr_nogradmu))
    # pvngm = plot(; xlabel="x", ylabel=L"V_s")
    # plot!(pvngm, out.duct_inner_x, duct_inner_vs_nogradmu; label="Duct Inner Surface")
    # plot!(pvngm, out.duct_outer_x, duct_outer_vs_nogradmu; label="Duct Outer Surface")
    # # plot!(pvngm, out.hub_x, hub_vs_nogradmu; label="Hub Surface")
    # savefig(pvngm, savepath * "$(ptot)_panels_lc-vs-without-gradmu.pdf")
    # savefig(pvngm, savepath * "$(ptot)_panels_lc-vs-without-gradmu.png")

    savefig(pv, savepath * "$(ptot)_panels_lc-velocity-comp.pdf")
    savefig(pcd, savepath * "$(ptot)_panels_lc-pressure-comp.pdf")
    savefig(pch, savepath * "$(ptot)_panels_lc-pressure-comp-hub.pdf")
    savefig(pr, savepath * "$(ptot)_panels_lc-circulation-comp.pdf")
    savefig(pv, savepath * "$(ptot)_panels_lc-velocity-comp.png")
    savefig(pcd, savepath * "$(ptot)_panels_lc-pressure-comp.png")
    savefig(pch, savepath * "$(ptot)_panels_lc-pressure-comp-hub.png")
    savefig(pr, savepath * "$(ptot)_panels_lc-circulation-comp.png")
end

perr = plot(; xaxis=:log, xlabel="Number of Body Panels", ylabel="Percent Errors")
plot!(perr, bpan, abs.(rotor_thrust_err); label="Rotor Thrust")
plot!(perr, bpan, abs.(torque_err); label="Rotor Torque")
plot!(perr, bpan, abs.(rotor_efficiency_err); label="Rotor Efficiency")
savefig(perr, savepath * "rotor-lc-errors-comp.pdf")
savefig(perr, savepath * "rotor-lc-errors-comp.png")
