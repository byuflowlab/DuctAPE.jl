#TODO: update for constant vortex case
#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/examples/body_only_doublet/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default_new.jl")
include(project_dir * "/visualize/visualize_geometry.jl")

# - load geometry - #
# read data file
# TODO: load smooth geometry
include(project_dir * "/test/data/naca_662-015.jl")
# put coordinates together
coordinates = [x_duct r_duct]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
panels = dt.generate_panels(coordinates)

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    body_panels=panels,
    coordinates=[coordinates],
    controlpoints=true,
    nodes=true,
    TEnodes=true,
    normals=true,
    savepath=savepath,
    filename="duct-geometry.pdf",
)

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#
# Define freestream on panels
Vinf = 1.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
gamb = dt.freestream_influence_vector(panels.normal, Vsmat)

# - Adding Kutta Condition - #
dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

# - LU decomp - #
LHSLU = la.lu!(LHS)

#---------------------------------#
#             Solving             #
#---------------------------------#
# - Use LinearAlgebra's ldiv! to solve linear system using factorized matrix
# note that gamb is overwritten to be the RHS vector and then overwritten to be the body strengths
la.ldiv!(LHSLU, gamb)

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
Vtot = dt.vfromvortexpanels(panels.controlpoint, panels.controlpoints, gamb)

### --- Steady Surface Pressure --- ###
cp = 1.0 .- (dt.norm.(eachrow(Vtot)) / Vinf) .^ 2

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pp = plot(; xlabel="x", ylabel=L"c_p", yflip=true)
plot!(
    pp,
    pressurexupper,
    pressureupper;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    label="exp outer",
)
plot!(
    pp,
    pressurexlower,
    pressurelower;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:dtriangle,
    label="exp inner",
)
xs = panels.controlpoint[:, 1]
plot!(pp, xs, cp; label="DuctTAPE")

savefig(savepath * "duct-pressure-comp.pdf")
