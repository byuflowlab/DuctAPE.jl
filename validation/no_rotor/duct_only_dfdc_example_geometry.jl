#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/no_rotor/figs/"
dispath =
    project_dir * "/../../Writing/dissertation/src/ductsolvercontents/ductsolverfigures/"

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

include(project_dir * "/validation/no_rotor/data/dfdc_example_geometry_vcp.jl")

ductcoords = reverse(dfdc_example_duct_nodes; dims=1)
npan = size(ductcoords, 1) - 1

# npans = [21, 31, 41, 51, 61, 71, 81, 91, 101, 151, 161, 201]
# for (i, npan) in enumerate(npans)
# println("N Panels = ", npan - 1)
# repanel = dt.repanel_airfoil(ductcoords; N=npan, normalize=false)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
# panels = dt.generate_panels([repanel])
panels = dt.generate_panels([ductcoords])

xn = panels.node[:, 1]
xcp = panels.controlpoint[:, 1]

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# Define freestream on panels
Vinf = 20.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, size(panels.controlpoint, 1)) # need velocity on each panel

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
AICn, AICt = dt.vortex_aic_boundary_on_boundary(
    panels.controlpoint,
    panels.normal,
    panels.tangent,
    panels.node,
    panels.nodemap,
    panels.influence_length,
)

# # - Add Trailing Edge Gap Panel Influences - #
# dt.add_te_gap_aic!(
#     AICn,
#     AICt,
#     panels.controlpoint,
#     panels.normal,
#     panels.tangent,
#     panels.tenode,
#     panels.teinfluence_length,
#     panels.tendotn,
#     panels.tencrossn,
#     panels.teadjnodeidxs,
# )

kids = [
    size(AICn)[1]+1 1
    size(AICn)[1]+1 size(AICn)[2]
]

LHS = zeros(size(AICn)[1] + 1, size(AICn)[2])

dt.add_kutta!(LHS, AICn, kids)

RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
push!(RHS, 0.0)

#---------------------------------#
#             Solving             #
#---------------------------------#
gamb = LHS \ RHS

pg = plot(xn, gamb; xlabel="x", ylabel="node strengths", label="")
savefig(pg, savepath * "duct-gammas.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

## --- Velocity Contributions --- ###

# get tangent
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

# add in body induced tangent velocity
Vtan .+= AICt * gamb

# add in jump term
jump = (gamb[1:(end - 1)] + gamb[2:end]) / 2
Vtan .-= jump / 2.0

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (Vtan / Vinf) .^ 2

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pp = plot(;
    xlabel="x",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pp, xcp, cp; label="DuctAPE", color=1)

savefig(savepath * "duct-pressure-comp-$(npan-1)-panels.pdf")
