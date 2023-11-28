#=
Isolated Duct Validation using data from Lewis annular airofils section.
Geometry is a NACA 66-015, with coordinates generated from OpenVSP, with the repeated LE point manually removed.
In addition, geometry was interpolated using a cosine spaced scheme
=#

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

include(project_dir * "/visualize/plots_default.jl")

# - load experimental data - #
include(project_dir * "/test/data/naca_662-015.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
coordinates = reverse(duct_coordinates; dims=1)

npan = 161

# interpolate geometry with cosine spacing and given number of panels
repanel = dt.repanel_airfoil(coordinates; N=npan, normalize=false)

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels([repanel])

# rename axial coordinates for convenience in later plotting
# NOTE: probably should use views
xn = panels.node[:, 1]
xcp = panels.controlpoint[:, 1]

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#

# Define freestream on panels
Vinf = 1.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, size(panels.controlpoint, 1)) # need velocity on each panel

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# NOTE: these two functions allocate more than they need to.

# - Initial System Matrices - #
@time AICn, AICt = dt.vortex_aic_boundary_on_boundary(
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

## -- Manually Assemble Linear System -- ##
# need to manually assemble since we don't have a way of automating the isolated case right now
# initialize LHS Matrix
LHS = zeros(size(AICn)[2] + 1, size(AICn)[2] + 1)
# fill LHS matrix with standard AIC terms
LHS[1:size(AICn, 1), 1:size(AICn, 2)] .= AICn
# add on pseudo control point influence terms
LHS[size(AICn, 2), 1:size(AICn, 2)] .= AICpcp'
# add in dummy variable influence terms
LHS[1:size(AICn, 1), size(AICn, 2) + 1] .= 1.0
# add kutta condition
LHS[size(AICn, 2) + 1, 1] = LHS[size(AICn, 2) + 1, size(AICn, 2)] = 1.0

# - assemble RHS - #
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
push!(RHS, -1.0)
push!(RHS, 0.0)

#---------------------------------#
#             Solving             #
#---------------------------------#
# TODO: use LU decomp here in real application
# standard linear solve
gamb = LHS \ RHS
# don't include dummy variable associated with internal panel
gamb = gamb[1:(end - 1)]

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
