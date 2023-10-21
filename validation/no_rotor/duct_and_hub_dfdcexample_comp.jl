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
include(project_dir * "/validation/no_rotor/data/dfdc_example_geometry_vcp.jl")
hubx = dfdc_example_hub[:, 1]
hubcp = dfdc_example_hub[:, 3]
hubvs = dfdc_example_hub[:, end]
ductx = dfdc_example_duct[:, 1]
ductcp = dfdc_example_duct[:, 3]
ductvs = dfdc_example_duct[:, end]

hubcoords = reverse(dfdc_example_hub_nodes; dims=1)
ductcoords = reverse(dfdc_example_duct_nodes; dims=1)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
coordinates = [ductcoords, hubcoords]
panels = dt.generate_panels(coordinates)

# - get raw DFDC outputs - #
include(project_dir * "/validation/no_rotor/data/dfdc_example_LHS.jl")
include(project_dir * "/validation/no_rotor/data/dfdc_example_RHS.jl")
gamdfdc = dfdc_example_lhs \ dfdc_example_rhs
dfdc_hub_gammas = gamdfdc[1:panels.nnode[2]]
dfdc_duct_gammas = gamdfdc[(panels.nnode[2] + 1):(end - 2)]
pgh = plot(; xlabel="x", ylabel="hub node strengths")
plot!(
    pgh,
    dfdc_example_hub_nodes[:, 1],
    dfdc_hub_gammas;
    linestyle=:dash,
    linewidth=2,
    label="DFDC",
    color=myred,
)
pgd = plot(; xlabel="x", ylabel="duct node strengths")
plot!(
    pgd,
    dfdc_example_duct_nodes[:, 1],
    dfdc_duct_gammas;
    linestyle=:dash,
    linewidth=2,
    label="DFDC",
    color=myred,
)

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

ductgam = gamb[1:panels.nnode[1]]
hubgam = gamb[(panels.nnode[1] + 1):(panels.totnode)]
plot!(pgd, ductcoords[:, 1], ductgam; label="DuctAPE", color=myblue)
plot!(pgh, hubcoords[:, 1], hubgam; label="DuctAPE", color=myblue)
savefig(pgd, savepath * "raw_gamma_comp_dfdc_example_duct.pdf")
savefig(pgh, savepath * "raw_gamma_comp_dfdc_example_hub.pdf")

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

## --- Velocity Contributions --- ###

# get tangent
Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

# add in body induced tangent velocity
Vtan .+= AICt * gamb[1:size(AICt, 2)]

# add in jump term
jumpduct = (gamb[1:(panels.npanel[1])] + gamb[2:(panels.nnode[1])]) / 2
jumphub =
    (
        gamb[(panels.nnode[1]):(panels.totpanel)] +
        gamb[(panels.nnode[1] + 1):(panels.totpanel + 1)]
    ) / 2
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
    "dfdcexample-pressure-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.pdf",
)
savefig(
    pcp,
    dispath *
    "dfdcexample-pressure-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.tikz",
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
    "dfdcexample-velocity-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.pdf",
)
savefig(
    pvs,
    dispath *
    "dfdcexample-velocity-comp-$(panels.npanel[1])-duct-panels-$(panels.npanel[2])-hub-panels.tikz",
)
