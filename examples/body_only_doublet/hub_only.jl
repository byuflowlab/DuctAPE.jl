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
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
coordinates = [x_hub r_hub]

#---------------------------------#
#             Paneling            #
#---------------------------------#
##### ----- Generate Panels ----- #####
#= use new Neumann paneling functions in panel.jl
generate_panels(coordinates::Vector{Matrix{TF}}) where {TF}
note that input must be a vector of matrices, so even if only modeling one body, still
needs to be a vector.
Also note, multiple dispatch enabled to allow for single body, simply calls the function
placing the coordinate matrix in a vector as an input.
=#
panels = dt.generate_panels(coordinates)

##### ----- Visualize to Check ----- #####
visualize_paneling(
    panels;
    coordinates=coordinates,
    controlpoints=true,
    nodes=true,
    TEnodes=true,
    normals=true,
    savepath=savepath,
    filename="hub-geometry.pdf",
)

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#
# Define freestream on panels
Vinf = 30.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel

# prescribe a panel for the least squares solve:
# choose the first panel to be prescirbed to zero (assumes first panel is not hub leading/traling edge).
prescribedpanels = [(1, 0.0)]

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
LHS = dt.init_body_lhs(panels)
RHS = dt.gen_body_rhs(panels.normal, Vsmat)

# - Prepping for Least Sqaures Solve - #
LHSlsq, RHSlsq = prep_leastsquares(LHS, RHS, prescribedpanels)

#---------------------------------#
#             Solving             #
#---------------------------------#

mured = LHSlsq \ RHSlsq

# - Solving With Kutta - #
mu = dt.mured2mu(mured, prescribedpanels)

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu)

# - "Wake"-induced Surface Velocity - #
Vb_wake = dt.vfromTE(panels.controlpoint, panels.TEnodes, panels.TEidxs, mu)

# - ∇μ/2 surface velocity - #
Vb_gradmu = dt.vfromgradmu(panels, mu)

# - Total Velocity - #
V_nogradmu = Vb .+ Vb_wake
Vtot = Vsmat .+ Vb .+ Vb_wake .+ Vb_gradmu

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
pp = plot(; xlabel="x", ylabel=L"Vs/V_\infty")
plot!(
    pp,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    color=blue[1],
    markershape=:utriangle,
    label="experimental",
)
xs = panels.controlpoint[:, 1]
plot!(pp, xs, dt.norm.(eachrow(Vtot)) ./ Vinf; label="DuctTAPE")

savefig(savepath * "hub-vel-comp.pdf")
# TODO: need to coordinate with Eduardo on who is adding what where for this stuff
