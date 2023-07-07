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

# # - load geometry - #
# # read data file
include(project_dir * "/test/data/naca_662-015.jl")
# # put coordinates together
# coordinates = [x_duct r_duct]

import GeometricTools as gt
using CSV
using DataFrames

filename = "/home/juddmehr/code-development/julia_devs/FLOWPanel.jl/examples/data/naca662015.csv"
contour = CSV.read(filename, DataFrame)
aspectratio = 0.6                       # Duct trailing edge aspect ratio l/d
d = 2 * 0.835                   # (m) duct diameter
n_rfl = 10                        # This controls the number of chordwise panels

NDIVS_rfl_up = [                            # Discretization of airfoil upper surface
    # 0 to 0.25 of the airfoil has `n_rfl` panels at a geometric expansion of 10 that is not central
    (0.25, n_rfl, 10.0, false),
    # 0.25 to 0.75 of the airfoil has `n_rfl` panels evenly spaced
    (0.50, n_rfl, 1.0, true),
    # 0.75 to 1.00 of the airfoil has `n_rfl` panels at a geometric expansion of 0.1 that is not central
    (0.25, n_rfl, 0.1, false),
]

NDIVS_rfl_lo = NDIVS_rfl_up                 # Discretization of airfoil lower surface
# Re-discretize the contour of the body of revolution according to NDIVS
xs, ys = gt.rediscretize_airfoil(
    contour[:, 1], contour[:, 2], NDIVS_rfl_up, NDIVS_rfl_lo; verify_spline=false
)
ys[end] = ys[1]

# Scale contour by duct length
xs *= d * aspectratio
ys *= d * aspectratio

# Move contour to the radial position
ys .+= d / 2

# Revert points to make the normals point outwards
reverse!(xs)
reverse!(ys)

coordinates = [xs ys]

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
visualize_paneling(;
    body_panels=panels,
    coordinates=coordinates,
    controlpoints=true,
    nodes=true,
    endpoints=true,
    normals=true,
    savepath=savepath,
    filename="duct-geometry.pdf",
)

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#
# Define freestream on panels
Vinf = 30.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel
# println("Uinfs: ")
# display(Vsmat)

# prescribe a panel for the least squares solve:
# choose the first panel to be prescirbed to zero (assumes first panel is not hub leading/traling edge).
prescribedpanels = [(1, 0.0)]

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
# LHS = dt.init_body_lhs(panels)
LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
# LHSnokutta = deepcopy(LHS)
# RHSnokutta = deepcopy(RHS)

# - Adding Kutta Condition - #
dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

# - Prepping for Least Sqaures Solve - #
# without kutta condition
# LHSlsq_nokutta, RHSlsq_nokutta = dt.prep_leastsquares(
#     LHSnokutta, RHSnokutta, prescribedpanels
# )
# with kutta condition
LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

#---------------------------------#
#             Solving             #
#---------------------------------#

mured = LHSlsq \ RHSlsq
# mured_nokutta = LHSlsq_nokutta \ RHSlsq_nokutta

# - Solving Without Kutta - #
# mu_nokutta = LHSlsq_nokutta\RHSlsq_nokutta
# mu_nokutta = dt.mured2mu(mured_nokutta, prescribedpanels)

# - Solving With Kutta - #
mu = dt.mured2mu(mured, prescribedpanels)

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
# Vb_nokutta = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu_nokutta)
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu)

# - "Wake"-induced Surface Velocity - #
# Vb_nokutta_wake = dt.vfromTE(
# panels.controlpoint, panels.endpoints, panels.endpointidxs, mu_nokutta
# )
Vb_wake = dt.vfromTE(panels.controlpoint, panels.TEnodes, mu; verbose=true)

# - ∇μ/2 surface velocity - #
Vb_gradmu = dt.vfromgradmu(panels, mu)

# - Total Velocity - #
# V_nokutta = Vb_nokutta .+ Vb_nokutta_wake
V_nogradmu = Vsmat .+ Vb .+ Vb_wake
Vtot = V_nogradmu .+ Vb_gradmu

# ### --- Velocity Tangent to Surface --- ###
# Vtan_nokutta = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nokutta), eachrow(panels.tangent))]
# Vtan_nogradmu = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nogradmu), eachrow(panels.tangent))]
vnngm = [dt.dot(v, n) for (v, n) in zip(eachrow(V_nogradmu), eachrow(panels.normal))]
vtngm = [dt.dot(v, t) for (v, t) in zip(eachrow(V_nogradmu), eachrow(panels.tangent))]
vmagngm = [dt.norm(v) for v in eachrow(V_nogradmu)]
vtcheck = abs.(vtngm) .- vmagngm

# vtotmag = dt.norm.(eachrow(Vtot))

### --- Steady Surface Pressure --- ###
# cp_nokutta = 1.0 .- (dt.norm.(eachrow(V_nokutta)) / Vinf) .^ 2
cp_nogradmu = 1.0 .- (dt.norm.(eachrow(V_nogradmu)) / Vinf) .^ 2
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
# plot!(xs,cp_nokutta,label="no Kutta")
# plot!(xs,cp_nogradmu,label=L"no~ \nabla\mu")
ncut = 4
plot!(pp, xs[ncut:(end - ncut)], cp[ncut:(end - ncut)]; label="DuctTAPE")

savefig(savepath * "duct-pressure-comp.pdf")
