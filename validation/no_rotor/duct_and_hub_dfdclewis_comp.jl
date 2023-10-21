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

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")

# - load DFDC outputs - #
include(project_dir * "/validation/no_rotor/data/dfdc_lewis_geometry_vcp.jl")
hubx = dfdc_lewis_hub[:, 1]
hubcp = dfdc_lewis_hub[:, 3]
hubvs = dfdc_lewis_hub[:, end]
ductx = dfdc_lewis_duct[:, 1]
ductcp = dfdc_lewis_duct[:, 3]
ductvs = dfdc_lewis_duct[:, end]

hubcoords = reverse(dfdc_lewis_hub_nodes; dims=1)
ductcoords = reverse(dfdc_lewis_duct_nodes; dims=1)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
# panels = dt.generate_panels([repanel];body=true)
panels = dt.generate_panels([ductcoords, hubcoords])

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

# -- Assemble LHS Matrix -- #
# - Boundary on boundary influence coefficients - #
AICn, AICt = dt.vortex_aic_boundary_on_boundary(
    panels.controlpoint,
    panels.normal,
    panels.tangent,
    panels.node,
    panels.nodemap,
    panels.influence_length,
)

# - Boundary on internal psuedo control point influence coefficients - #
AICpcp, _ = dt.vortex_aic_boundary_on_field(
    panels.itcontrolpoint,
    panels.itnormal,
    panels.ittangent,
    panels.node,
    panels.nodemap,
    panels.influence_length,
)

# - Freetream influence for RHS vector - #
vdnb = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel))
vdnpcp = dt.freestream_influence_vector(
    panels.itnormal, repeat(Vs, size(panels.itcontrolpoint, 1))
)

LHS = dt.assemble_lhs_matrix(AICn, AICpcp, panels; dummyval=1.0)
RHS = dt.assemble_rhs_matrix(vdnb, vdnpcp, panels)
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
Vtan .+= AICt * gamb[1:size(AICt, 2)]

# add in jump term
jumpduct = (gamb[1:(panels.npanel[1])] + gamb[2:panels.nnode[1]]) / 2
jumphub = (gamb[(panels.nnode[1] + 1):(panels.totpanel+1)]+gamb[(panels.nnode[1]+2):(panels.totnode)]) / 2
jump = [jumpduct; jumphub]
Vtan .-= jump / 2.0

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (Vtan / Vinf) .^ 2

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pcp = plot(;
    size=(225, 169),
    xlabel="x",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pcp, xcp[1:panels.npanel[1]], cp[1:panels.npanel[1]]; label="", color=1)
plot!(pcp, ductx, ductcp; linewidth=2, linestyle=:dash, color=2, label="")

savefig(
    pcp,
    savepath *
    "dfdclewis-pressure-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.pdf",
)
savefig(
    pcp,
    dispath *
    "dfdclewis-pressure-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.tikz",
)

pvs = plot(;
    size=(225, 169),
    xlabel="x",
    ylabel=L"\frac{V_s}{V_\infty}",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)

plot!(
    pvs,
    xcp[panels.nnode[1]:(panels.totpanel)],
    Vtan[panels.nnode[1]:(panels.totpanel)] ./ Vinf;
    color=1,
    label="",
)

plot!(pvs, hubx, hubvs ./ Vinf; linewidth=2, linestyle=:dash, color=2, label="")

savefig(
    pvs,
    savepath *
    "dfdclewis-velocity-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.pdf",
)
savefig(
    pvs,
    dispath *
    "dfdclewis-velocity-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.tikz",
)

