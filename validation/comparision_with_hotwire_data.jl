#---------------------------------#
#              SETUP              #
#---------------------------------#

project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

savepath = joinpath(project_dir, "validation_scripts/", "outputs/")

using DelimitedFiles
using Statistics
using FLOWMath
const fm = FLOWMath
using DuctTAPE
const dt = DuctTAPE
using LineSearches
using ForwardDiff
using NLsolve

# plotting defaults
include("../plots_default_new.jl")
pyplot()

i2m = 0.0254
ft2m = 12.0 * i2m

#---------------------------------#
#              RUN DT             #
#---------------------------------#
# - Read in Parameters - #
include("ducttape_parameters.jl")

# optionally all together
out, converged_strengths, inputs, initial_strengths, convergeflag = @time dt.analyze_propulsor(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotorstator_parameters,
    freestream,
    reference_parameters;
    verbose=true, #prints out solve progress
    maximum_linesearch_step_size=1e6, # used in Newton Solve
    iteration_limit=15, # used in Newton Solve
)

#---------------------------------#
#            Get Angles           #
#---------------------------------#
mub0, gamw0, Gamr0, sigr0 = dt.extract_state_variables(converged_strengths, inputs)

overwriterpm = range(1000, 7808; step=25)
overwriteomega = overwriterpm .* pi / 30
for (io, omega) in enumerate(overwriteomega)
    println("Running at RPM = ", overwriterpm[io])

    global initial_strengths, converged_strengths, inputs

    overwrite_params = [(; rotor_parameters..., Omega=omega), stator_parameters]
    # repanels bodies and rotors, generates wake "grid", precomputes influence matrices, etc.
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        overwrite_params,
        freestream,
        reference_parameters;
    )

    # calculate initial guess for state variables
    initial_strengths = converged_strengths
    initials = copy(initial_strengths)

    # - Define closure that allows for parameters - #
    p = []
    rwrap(r, states) = dt.residual!(r, states, inputs, p)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        initial_strengths;
        autodiff=:forward,
        method=:newton,
        iterations=25,
        show_trace=true,
        linesearch=BackTracking(; maxstep=1e6),
        ftol=1e-8,
    )

    converged_strengths = res.zero
end

##---------------------------------#
##            Get Angles           #
##---------------------------------#
#mub0, gamw0, Gamr0, sigr0 = dt.extract_state_variables(converged_strengths, inputs)

## Get velocities at rotor
#vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = dt.calculate_rotor_velocities(
#    Gamr0, gamw0, sigr0, mub0, inputs
#)

## get the inflow and attack angles
#rotor_inflow, rotor_aoa = dt.calculate_inflow_angles(
#    Wm_rotor[:, 1], Wtheta_rotor[:, 1], inputs.blade_elements[1].twists
#)

#stator_inflow, stator_aoa = dt.calculate_inflow_angles(
#    Wm_rotor[:, 2], Wtheta_rotor[:, 2], inputs.blade_elements[2].twists
#)

##---------------------------------#
##              PLOT               #
##---------------------------------#
#rotor_r = inputs.rotor_panel_centers[:, 1] ./ inputs.reference_parameters.Rref
#stator_r = inputs.rotor_panel_centers[:, 2] ./ inputs.reference_parameters.Rref

#pr = plot(;
#    title="Rotor Angles, rotor RPM=$RPM_rotor",
#    xlabel="Angles (degrees)",
#    ylabel="Normalized Radial Position",
#)
#plot!(pr, inputs.blade_elements[1].twists * 180 / pi, rotor_r; label="Twist")
#plot!(pr, rotor_inflow * 180 / pi, rotor_r; label="Inflow")
#plot!(pr, rotor_aoa * 180 / pi, rotor_r; label="Attack")
#savefig(pr, savepath * "rotor_initial_angles$RPM_rotor.pdf")
#savefig(pr, savepath * "rotor_initial_angles$RPM_rotor.png")

#ps = plot(;
#    title="Stator Angles, rotor RPM=$RPM_rotor",
#    xlabel="Angles (degrees)",
#    ylabel="Normalized Radial Position",
#)
#plot!(ps, inputs.blade_elements[2].twists * 180 / pi, stator_r; label="Twist")
#plot!(ps, stator_inflow * 180 / pi, stator_r; label="Inflow")
#plot!(ps, stator_aoa * 180 / pi, stator_r; label="Attack")
#savefig(ps, savepath * "stator_initial_angles$RPM_rotor.pdf")
#savefig(ps, savepath * "stator_initial_angles$RPM_rotor.png")
## plot re-paneled geometry
#include("../../DuctTAPE.jl/visualize/visualize_geometry.jl")
#visualize_paneling(;
#    body_panels=inputs.body_doublet_panels,
#    rotor_panels=inputs.rotor_source_panels,
#    wake_panels=inputs.wake_vortex_panels,
#    coordinates=[duct_coordinates, hub_coordinates],
#    controlpoints=true,
#    nodes=true,
#    wakeinterfaceid=[],
#    prescribedpanels=nothing,
#    TEnodes=false,
#    normals=false,
#    normal_scaling=0.05,
#    savepath=savepath,
#    filename=["NASA_geometry_repanel.pdf"],
#    legendloc=false,
#    zoom=false,
#    limits=nothing,
#    nodemarkersize=1,
#    cpmarkersize=1,
#)

# # - Generate DFDC Case File - #
# include("../../DuctTAPE.jl/convenience_functions/generate_dfdc_case.jl")
# gen_dfdc_case(
#     filename,
#     op_data,
#     wake_data,
#     airfoil_data,
#     rotor_data,
#     reverse(inputs.hub_coordinates; dims=1),
#     reverse(inputs.duct_coordinates; dims=1);
#     savepath="temp_debug/",
# )

# Define Probe Locations
Nprobes = 25
# x location = 4 inches aft of rotor stacking axis, which is at x=0
probe_poses_x = 4.0 * 0.0254 * ones(Nprobes) #4 inches behind rotor stacking axis
# interpolate to get radial location of hub and duct
probetip = FLOWMath.akima(inner_duct[:, 1], inner_duct[:, 2], probe_poses_x[1])
probehub = FLOWMath.akima(
    hub_coordinates_closed_TE[:, 1], hub_coordinates_closed_TE[:, 2], probe_poses_x[1]
)
# r locations evenly spaced
probe_poses_r = range(probehub, probetip, Nprobes)

# get probed velocities
dt_vx_init, dt_vr_init, dt_vtheta_init = dt.probe_velocity_field(
    [probe_poses_x probe_poses_r], inputs, initial_strengths; debug=false
)

dt_vx, dt_vr, dt_vtheta = dt.probe_velocity_field(
    [probe_poses_x probe_poses_r], inputs, converged_strengths
)

#---------------------------------#
#         EXP HOT-WIRE DATA       #
#---------------------------------#

# - Get inner duct radius for non-dimensionalization - #
# use inner_duct coordinates extracted in parameters file
Rduct = FLOWMath.akima(inner_duct[:, 1], inner_duct[:, 2], 4.0 * i2m)

# - Load Data Files - #
rexp = readdlm("../Hotwire_RMS/Immersion_Station.data")[:, 1] .* i2m ./ Rduct
vxexp = readdlm("../Hotwire_RMS/U_1.data")
vrexp = readdlm("../Hotwire_RMS/W_1.data")
vtexp = readdlm("../Hotwire_RMS/V_1.data")

# - Vx for plotting - #
vxexpbar = zeros(length(rexp))
vxexpmax = zeros(length(rexp))
vxexpmin = zeros(length(rexp))
for (ir, r) in enumerate(rexp)
    vxexpbar[ir] = Statistics.mean(vxexp[:, ir]) * ft2m ./ Vinf
    vxexpmax[ir] = maximum(vxexp[:, ir]) * ft2m ./ Vinf
    vxexpmin[ir] = minimum(vxexp[:, ir]) * ft2m ./ Vinf
end
vxerr = (vxexpbar .- vxexpmin, vxexpmax .- vxexpbar)

# - Vr for plotting - #
vrexpbar = zeros(length(rexp))
vrexpmax = zeros(length(rexp))
vrexpmin = zeros(length(rexp))
for (ir, r) in enumerate(rexp)
    vrexpbar[ir] = Statistics.mean(vrexp[:, ir]) * ft2m ./ Vinf
    vrexpmax[ir] = maximum(vrexp[:, ir]) * ft2m ./ Vinf
    vrexpmin[ir] = minimum(vrexp[:, ir]) * ft2m ./ Vinf
end
vrerr = (vrexpbar .- vrexpmin, vrexpmax .- vrexpbar)

# - Vtheta for plotting - #
vtexpbar = zeros(length(rexp))
vtexpmax = zeros(length(rexp))
vtexpmin = zeros(length(rexp))
for (ir, r) in enumerate(rexp)
    vtexpbar[ir] = Statistics.mean(vtexp[:, ir]) * ft2m ./ Vinf
    vtexpmax[ir] = maximum(vtexp[:, ir]) * ft2m ./ Vinf
    vtexpmin[ir] = minimum(vtexp[:, ir]) * ft2m ./ Vinf
end
vterr = (vtexpbar .- vtexpmin, vtexpmax .- vtexpbar)

#---------------------------------#
#              PLOTS              #
#---------------------------------#

# - Axial Velocity Comparison - #
px = plot(;
    xlabel=L"\mathrm{Axial~Velocity}~\left(\frac{V_x}{V_\infty}\right)",
    ylabel=L"\mathrm{Radial~Position}~\left(\frac{r}{R_\mathrm{duct}}\right)",
)

# hot-wire axial velocity

plot!(
    px,
    vxexpbar,
    rexp;
    color=myblue[1],
    xerror=vxerr,
    markerstrokecolor=mygray[3],
    label="Experimental",
)

# DuctTAPE axial velocity
# initial
plot!(
    px,
    dt_vx_init ./ inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE Initial States",
)

# converged
plot!(
    px,
    dt_vx ./ inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE Converged States",
)

# - Radial Velocity Comparison - #
pr = plot(;
    xlabel=L"\mathrm{Radial~Velocity}~\left(\frac{V_r}{V_\infty}\right)",
    ylabel=L"\mathrm{Radial~Position}~\left(\frac{r}{R_\mathrm{duct}}\right)",
)

# hot-wire radial velocity
plot!(
    pr,
    vrexpbar,
    rexp;
    color=myblue[1],
    xerror=vrerr,
    markerstrokecolor=mygray[3],
    label="Experimental",
)

# DuctTAPE radial velocity
# initial
plot!(
    pr,
    dt_vr_init .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE Initial States",
)

# converged
plot!(
    pr,
    dt_vr .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE",
)

# - Tangential Velocity Comparison - #
pt = plot(;
    xlabel=L"\mathrm{Tangential~Velocity}~\left(\frac{V_\theta}{V_\infty}\right)",
    ylabel=L"\mathrm{Radial~Position}~\left(\frac{r}{R_\mathrm{duct}}\right)",
)

# hot-wire tangential velocity
plot!(
    pt,
    vtexpbar,
    rexp;
    color=myblue[1],
    xerror=vterr,
    markerstrokecolor=mygray[3],
    label="Experimental",
)

# DuctTAPE tangential velocity
# initial
plot!(
    pt,
    dt_vtheta_init .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE Initial States",
)

# converged
plot!(
    pt,
    dt_vtheta .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctTAPE",
)

savefig(px, "figures/axial_velocity_comparison.pdf")
savefig(pr, "figures/radial_velocity_comparison.pdf")
savefig(pt, "figures/tangential_velocity_comparison.pdf")

savefig(px, "figures/axial_velocity_comparison.png")
savefig(pr, "figures/radial_velocity_comparison.png")
savefig(pt, "figures/tangential_velocity_comparison.png")

