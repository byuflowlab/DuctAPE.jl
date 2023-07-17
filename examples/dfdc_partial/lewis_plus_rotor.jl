project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctTAPE
const dt = DuctTAPE
using PrettyTables
const pt = PrettyTables

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/visualize_flowfield.jl")
include(project_dir * "/visualize/plots_default_new.jl")
# pyplot()

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 50.0))
pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-2.25, 1.0), yflip=true)

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

plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
# plot!(pv, hubwakex, hubwakevs; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pv, ductwakex, ductwakevs; linestyle=:dash, label="DFDC Duct Wake")

plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
# plot!(pc, hubwakex, hubwakecp; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pc, ductwakex, ductwakecp; linestyle=:dash, label="DFDC Duct Wake")

dfdc_duct_coordinates = reverse([ductx ductr]; dims=1)
dfdc_hub_coordinates = reverse([hubx hubr]; dims=1)
dfdc_panels = dt.generate_panels([dfdc_duct_coordinates, dfdc_hub_coordinates])

#viz dfdc geometry
visualize_paneling(;
    body_panels=dfdc_panels,
    coordinates=[dfdc_duct_coordinates, dfdc_hub_coordinates],
    controlpoints=true,
    nodes=true,
    normals=true,
    normal_scaling=0.1,
    savepath=savepath,
    filename="lewis-with-rotor-dfdcbodygeometry.pdf",
    legendloc=:right,
)

#---------------------------------#
#             Exp data            #
#---------------------------------#
include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")

# plot!(
#     pc,
#     pressurexupper,
#     pressureupper;
#     seriestype=:scatter,
#     color=myblue[1],
#     markershape=:utriangle,
#     markersize=3,
#     label="Exp Outer Duct",
# )

# plot!(
#     pc,
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
savefig(pv, savepath * "lewis-with-rotor-velocity-dfdc.pdf")
savefig(pc, savepath * "lewis-with-rotor-pressure-dfdc.pdf")

#---------------------------------#
#             DuctTAPE            #
#---------------------------------#
# include(datapath * "lewis_refined_duct.jl")
# include(datapath * "lewis_refined_hub.jl")

# - Generate Panels - #
# TODO: create parameters file to read in
include(datapath * "lewis_with_rotor.case.jl")

#overwrite geometry with DFDC input geometry
duct_coordinates = reverse([ductx ductr]; dims=1)
duct_coordinates[1, :] .= duct_coordinates[end, :]
hub_coordinates = reverse([hubx hubr]; dims=1)

#overwite npanels for now (should put conditions for this in auto generation function)
paneling_constants = (;
    paneling_constants..., nhub_inlet=30, nduct_inlet=30, npanels=[30, 20, 40]
)

# # run inputs first to check things
# # initialize various inputs used in analysis
# inputs = dt.precomputed_inputs(
#     duct_coordinates,
#     hub_coordinates,
#     paneling_constants,
#     rotor_parameters,
#     freestream,
#     reference_parameters;
# )

# # calculate initial guess for state variables
# initial_states = dt.initialize_states(inputs)

# # - Define closure that allows for parameters - #
# p = (;)
# rwrap(r, statesvars) = dt.residual!(r, statesvars, inputs, p)

# res = dt.NLsolve.nlsolve(
#     rwrap,
#     initial_states;
#     autodiff=:forward,
#     method=:newton,
#     iterations=50,
#     show_trace=true,
#     linesearch=dt.BackTracking(; maxstep=1e6),
# )

# states = res.zero

# out = dt.post_process(states, inputs)

# # run analyze_propulsor function
# out, converged_states, inputs, initial_states, convergeflag = dt.analyze_propulsor(
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

#---------------------------------#
#           Comparisons           #
#---------------------------------#
## -- Print Post Processed Comparsion Values -- ##
include(datapath * "dfdc_outs/lewis_coupled/DFDC_FORCES.jl")

function perr(x, y)
    return 100.0 * (y - x) / x
end

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
    [outcomps; reduce(vcat, outscomprot)]; header=["Value", "DFDC", "DuctTAPE", "% Error"]
)

## -- Plotting -- ##

println("Visualizing Paneling")
# check normals
visualize_paneling(;
    body_panels=inputs.body_doublet_panels,
    coordinates=[duct_coordinates, hub_coordinates],
    controlpoints=true,
    nodes=true,
    normals=true,
    normal_scaling=0.1,
    savepath=savepath,
    filename="lewis-with-rotor-bodygeometry.pdf",
    legendloc=:right,
)

# everything
visualize_paneling(;
    body_panels=inputs.body_doublet_panels,
    rotor_panels=inputs.rotor_source_panels,
    wake_panels=inputs.wake_vortex_panels,
    coordinates=[duct_coordinates, hub_coordinates],
    controlpoints=true,
    nodes=true,
    normals=false,
    normal_scaling=0.1,
    savepath=savepath,
    filename="lewis-with-rotor-fullgeometry.pdf",
    legendloc=:outerright,
)

# everything
visualize_paneling(;
    body_panels=inputs.body_doublet_panels,
    rotor_panels=inputs.rotor_source_panels,
    wake_panels=inputs.wake_vortex_panels,
    wakeinterfaceid=inputs.ductwakeinterfaceid,
    coordinates=[duct_coordinates, hub_coordinates],
    controlpoints=true,
    nodes=true,
    TEnodes=false,
    normals=false,
    normal_scaling=0.1,
    savepath=savepath,
    filename="lewis-with-rotor-fullgeometry-zoom.pdf",
    legendloc=:outerright,
    limits=(; ylim=(0.65, 0.8), xlim=(0.4, 0.6)),
    zoom=true,
)

println("Generating VTKs of Velocity Field")
# Flowfield VTK generation
visualize_flowfield(
    inputs.Vinf;
    body_panels=inputs.body_doublet_panels,
    rotor_panels=inputs.rotor_source_panels,
    wake_panels=inputs.wake_vortex_panels,
    mub=out.mub,
    sigr=out.sigr,
    gamw=out.gamw,
    # Pmax=nothing,
    # Pmin=nothing,
    verbose=true,
    run_name="velocity_fields",
    save_path=savepath,
    cellsizescale=0.005,
)

visualize_surfaces(
    inputs.Vinf;
    body_panels=inputs.body_doublet_panels,
    rotor_panels=inputs.rotor_source_panels,
    wake_panels=inputs.wake_vortex_panels,
    mub=out.mub,
    sigr=out.sigr,
    gamw=out.gamw,
    run_name="velocity_field",
    save_path=savepath,
)

# - Plot surface Velocity and Pressure - #
plot!(pv, out.duct_inner_x, out.duct_inner_vs; color=myblue[2], label="DuctTAPE Duct inner")
plot!(pv, out.duct_outer_x, out.duct_outer_vs; color=myblue[3], label="DuctTAPE Duct outer")
plot!(pv, out.hub_x, out.hub_vs; color=myred[2], label="DuctTAPE Hub")
plot!(pc, out.duct_inner_x, out.duct_inner_cp; color=myblue[2], label="DuctTAPE Duct inner")
plot!(pc, out.duct_outer_x, out.duct_outer_cp; color=myblue[3], label="DuctTAPE Duct outer")
plot!(pc, out.hub_x, out.hub_cp; color=myred[2], label="DuctTAPE Hub")

#- plot rotor circulation - #
pr = plot(; xlabel="Circulation", ylable="r")
for ir in inputs.num_rotors
    plot!(
        pr,
        dfdc_circulation / inputs.blade_elements[ir].B,
        inputs.blade_elements[ir].rbe;
        linestyle=:dash,
        color=myblue[ir],
        label="DFDC Rotor $ir",
    )
    plot!(
        pr,
        out.Gamr[:, ir],
        inputs.blade_elements[ir].rbe;
        color=myblue[ir],
        label="DuctTAPE Rotor $ir",
    )
end

savefig(pv, savepath * "lewis-with-rotor-velocity-comp.pdf")
savefig(pc, savepath * "lewis-with-rotor-pressure-comp.pdf")
savefig(pr, savepath * "lewis-with-rotor-circulation-comp.pdf")
savefig(pv, savepath * "lewis-with-rotor-velocity-comp.png")
savefig(pc, savepath * "lewis-with-rotor-pressure-comp.png")
