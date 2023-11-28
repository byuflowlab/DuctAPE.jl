######################################################################
#                                                                    #
#                        How to run DuctAPE                         #
#                                                                    #
######################################################################
# - Load Package - #
using DuctAPE
const dt = DuctAPE

# - Read in Parameters - #
# most of the setup happens here.
include("ducttape_parameters.jl")

# - Run DuctAPE - #
out, converged_strengths, inputs, initial_strengths, convergeflag = @time dt.analyze_propulsor(
    duct_coordinates, # [x,r] duct coordinates, clockwise starting at TE
    hub_coordinates, # [x,r] center body coordinates, clockwise starting at LE
    paneling_constants, # named tuple including values used in repaneling geometry
    rotorstator_parameters, # vector of named tuples containing parameters for rotor and stator geometry and operating conditions
    freestream, # freestream values, velocity, density, viscosity, speed of sound
    reference_parameters; # reference values used in post-process, velocity, fan radius
    verbose=true, #prints out solve progress
    maximum_linesearch_step_size=1e6, # used in Newton Solve
    iteration_limit=15, # used in Newton Solve
)

#---------------------------------#
#          Probe Velocity         #
#---------------------------------#
# this is an extra post-processing step to get the velocity at the hot-wire locations

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

# get probed velocities for converged values
dt_vx, dt_vr, dt_vtheta = dt.probe_velocity_field(
    [probe_poses_x probe_poses_r], inputs, converged_strengths
)

######################################################################
#                                                                    #
#    Everything after this is just loading data and plotting         #
#                                                                    #
######################################################################

#---------------------------------#
#              SETUP              #
#---------------------------------#

# this just identifies the directory one up as the project directory, which helps in loading in files from various locations.
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# define save path
savepath = joinpath(project_dir, "validation_scripts/", "outputs/")

# additional Packages to Include
using DelimitedFiles
using Statistics
using FLOWMath
const fm = FLOWMath

# plotting defaults
include("../plots_default_new.jl")
pyplot()

i2m = 0.0254
ft2m = 12.0 * i2m

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

# DuctAPE axial velocity
plot!(
    px,
    dt_vx ./ inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctAPE",
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

# DuctAPE radial velocity
plot!(
    pr,
    dt_vr .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctAPE",
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

# DuctAPE tangential velocity
plot!(
    pt,
    dt_vtheta .* inputs.reference_parameters.Vref,
    probe_poses_r ./ Rduct;
    seriestype=:scatter,
    markershape=:utriangle,
    markersize=5,
    label="DuctAPE",
)

savefig(px, "figures/axial_velocity_comparison.pdf")
savefig(pr, "figures/radial_velocity_comparison.pdf")
savefig(pt, "figures/tangential_velocity_comparison.pdf")
