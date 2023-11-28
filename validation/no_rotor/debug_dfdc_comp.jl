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

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

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

plot(hubcoords[:,1], hubcoords[:,2], aspectratio=1, color=1)
plot!(ductcoords[:,1], ductcoords[:,2], aspectratio=1, color=1)

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
AICn, AICt = dt.vortex_aic_boundary_on_boundary(
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
#TODO: these aren't going to match, you integrate in different directions.
LHSductduct = LHS[1:Nnodeduct, 1:Nnodeduct]
LHSductduct[end, :] .= LHS[end - 1, 1:Nnodeduct]
LHSductduct

include(project_dir * "/validation/no_rotor/data/dfdc_example_lhs_duct-on-duct.jl")
LHSductductdfdc = dfdclhsductduct
maxdderr = maximum(LHSductduct .- LHSductductdfdc)
println("\nMaximum LHS Duct-Duct Error: ", maxdderr)

# - Duct On Hub - #
#TODO: these aren't going to match, you integrate in different directions.
LHSducthub = LHS[Nnodeduct:(end - 2), 1:Nnodeduct]
LHSducthub

include(project_dir * "/validation/no_rotor/data/dfdc_example_lhs_duct-on-hub.jl")
LHSducthubdfdc = dfdclhsductonhub[2:end-1,1:end-1]
maxdherr = maximum(LHSducthub .- LHSducthubdfdc)
println("\nMaximum LHS Duct-on-Hub Error: ", maxdherr)

# - Hub On Hub - #
#TODO: these aren't going to match, you integrate in different directions.
#TODO: flipping the matrix is not the same as integrating the other direction...
LHShubhub = LHS[Nnodeduct:(end - 2), (Nnodeduct + 1):end]
LHShubhubflipv = reverse(LHShubhub,dims=1)
LHShubhubfliph = reverse(LHShubhubflipv,dims=2)
LHShubhubfliph
LHShubhub = [LHShubhubfliph; reverse(LHS[end, (Nnodeduct + 1):end])']

include(project_dir * "/validation/no_rotor/data/dfdc_example_lhs_hub-on-hub.jl")
LHShubhubdfdc = dfdclhshubhub[2:end,1:end-1]
maxhherr = findmax(LHShubhub[2:end-1,:] .- LHShubhubdfdc[2:end-1,:])
println("\nMaximum LHS Duct-Duct Error: ", maxhherr)

# - Hub On Duct - #
#TODO: these aren't going to match, you integrate in different directions.
# TODO


# - RHS - #
include(project_dir * "/validation/no_rotor/data/dfdc_example_rhs.jl")
dtrhs = reverse(RHS[1:end-2])
insert!(dtrhs,62,0.0)
push!(dtrhs,0.0)
maxrhserr1 = findmax(dtrhs[1:62] .+ dfdc_rhs[1:62])
maxrhserr2 = findmax(dtrhs[63:end] .+ dfdc_rhs[64:end-1])
println("\nMaximum RHS Error: ", max(maxrhserr1,maxrhserr2))
