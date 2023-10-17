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

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

# - load DFDC outputs - #
include(project_dir * "/validation/no_rotor/data/dfdc_lewis_geometry_vcp.jl")
hubx = dfdc_lewis_hub[:, 1]
hubcp = dfdc_lewis_hub[:, 3]
hubvs = dfdc_lewis_hub[:, end]
ductx = dfdc_lewis_duct[:, 1]
ductcp = dfdc_lewis_duct[:, 3]
ductvs = dfdc_lewis_duct[:, end]

hubcoords = reverse(dfdc_lewis_hub_nodes[2:end, :]; dims=1)
ductcoords = reverse(dfdc_lewis_duct_nodes; dims=1)
npanduct = size(ductcoords)[1]
npanhub = size(hubcoords)[1]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
# panels = dt.generate_panels([repanel];body=true)
panels = dt.generate_panels([ductcoords, hubcoords])

# ##### ----- Visualize to Check ----- #####
# visualize_paneling(;
#     body_panels=panels,
#     coordinates=[coordinates],
#     controlpoints=true,
#     nodes=true,
#     normals=true,
#     savepath=savepath,
#     filename=["duct-geometry.pdf"],
# )

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
AICn, AICt = dt.vortex_panel_influence_matrices(
    panels.controlpoint,
    panels.normal,
    panels.tangent,
    panels.node,
    panels.nodemap,
    panels.influence_length,
)

kids = [
    size(AICn)[1]+1 1
    size(AICn)[1]+1 npanduct
    size(AICn)[1]+2 npanduct+1
]

LHS = zeros(size(AICn)[1] + 2, size(AICn)[2])

dt.add_kutta!(LHS, AICn, kids)

RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
push!(RHS, 0.0)
push!(RHS, 0.0)

#---------------------------------#
#             Solving             #
#---------------------------------#
gamb = LHS \ RHS

# pg = plot(xn, gamb; xlabel="x", ylabel="node strengths", label="")
# savefig(pg, savepath * "duct-gammas.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

## --- Velocity Contributions --- ###

# get tangent
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

# add in body induced tangent velocity
Vtan .+= AICt * gamb

# add in jump term
jumpduct = (gamb[1:(npanduct - 1)] + gamb[2:npanduct]) / 2
jumphub = (gamb[(npanduct + 1):(end - 1)] + gamb[(npanduct + 2):end]) / 2
jump = [jumpduct; jumphub]
Vtan .-= jump / 2.0

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (Vtan / Vinf) .^ 2

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pcp = plot(;
    xlabel="x",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pcp, xcp[1:(npanduct - 1)], cp[1:(npanduct - 1)]; label="DuctAPE", color=1)
plot!(pcp, ductx, ductcp; linestyle=:dash, color=2, label="DFDC")

savefig(
    pcp,
    savepath *
    "dfdclewis-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
)
savefig(
    pcp,
    savepath *
    "dfdclewis-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
)

pvs = plot(;
    xlabel="x",
    ylabel=L"\frac{V_s}{V_\infty}",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)

plot!(pvs, xcp[(npanduct):end], Vtan[(npanduct):end] ./ Vinf; color=1, label="DuctAPE")

plot!(pvs, hubx, hubvs ./ Vinf; linestyle=:dash, color=2, label="DFDC")

savefig(
    pvs,
    savepath *
    "dfdclewis-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
)
savefig(
    pvs,
    savepath *
    "dfdclewis-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
)
