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

Nnodeduct = size(ductcoords)[1]
Nnodehub = size(hubcoords)[1]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels([ductcoords, hubcoords])

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
    size(AICn)[1]+1 Nnodeduct
    size(AICn)[1]+2 Nnodeduct+1
]

LHS = zeros(size(AICn)[1] + 2, size(AICn)[2])

dt.add_kutta!(LHS, AICn, kids)

RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
push!(RHS, 0.0)
push!(RHS, 0.0)

# - Duct On Duct - #
LHSductduct = LHS[1:Nnodeduct, 1:Nnodeduct]
LHSductduct[end, :] .= LHS[end - 1, 1:Nnodeduct]
LHSductduct

include(project_dir * "/validation/no_rotor/data/dfdc_lhs_duct-on-duct.jl")
LHSductductdfdc = dfdclhsductduct[:, 1:(end - 1)]
maxdderr = maximum(LHSductduct .- LHSductductdfdc)
println("\nMaximum LHS Duct-Duct Error: ", maxdderr)

# - Duct On Hub - #
LHSducthub = LHS[Nnodeduct:(end - 2), 1:Nnodeduct]
LHSducthub

include(project_dir * "/validation/no_rotor/data/dfdc_lhs_duct-on-hub.jl")
LHSducthubdfdc = dfdclhsductonhub[2:end-1,1:end-1]
maxdherr = maximum(LHSducthub .- LHSducthubdfdc)
println("\nMaximum LHS Duct-on-Hug Error: ", maxdherr)

# - Hub On Hub - #
LHShubhub = LHS[Nnodeduct:(end - 2), (Nnodeduct + 1):end]
LHShubhubflipv = reverse(LHShubhub,dims=1)
LHShubhubfliph = reverse(LHShubhubflipv,dims=2)
LHShubhubfliph
LHShubhub = [LHShubhubfliph; reverse(LHS[end, (Nnodeduct + 1):end])']

include(project_dir * "/validation/no_rotor/data/dfdc_lhs_hub-on-hub.jl")
LHShubhubdfdc = dfdclhshubhub[2:end,1:end-1]
maxhherr = findmax(LHShubhub[2:end-1,:] .- LHShubhubdfdc[2:end-1,:])
println("\nMaximum LHS Duct-Duct Error: ", maxhherr)

# - Hub On Duct - #
# TODO


# - RHS - #
include(project_dir * "/validation/no_rotor/data/dfdc_rhs.jl")
dtrhs = reverse(RHS[1:end-2])
insert!(dtrhs,1,0.0)
insert!(dtrhs,80,0.0)
insert!(dtrhs,188,0.0)
maxrhserr = findmax(dtrhs .+ dfdc_rhs[1:end-1])
println("\nMaximum LHS Duct-Duct Error: ", maxrhserr)


