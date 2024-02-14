#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/dev_debug_archive/body_only_doublet/"

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")

# - load geometry - #
# read duct data file
include(project_dir * "/test/data/naca_662-015.jl")
# put coordinates together
duct_coordinates = [x_duct r_duct]

# read hub data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
hub_coordinates = [x_hub r_hub]

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
panels = dt.generate_panels([duct_coordinates, hub_coordinates])

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    body_panels=panels,
    coordinates=duct_coordinates,
    controlpoints=true,
    nodes=true,
    normals=true,
    savepath=savepath,
    filename="body-geometry.pdf",
)

#---------------------------------#
#       Operating Conditions      #
#---------------------------------#
# Define freestream on panels
Vinf = 30.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.totpanel) # need velocity on each panel

# prescribe a panel for the least squares solve:
# choose the first panel to be prescirbed to zero (assumes first panel is not hub leading/traling edge).
prescribedpanels = [(1, 0.0); (panels.totpanel, 0.0)]
# prescribedpanels = [(1, 0.0)]

#---------------------------------#
#        Induced Velocities       #
#---------------------------------#

# - Initial System Matrices - #
LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
LHSnokutta = deepcopy(LHS)
RHSnokutta = deepcopy(RHS)

# - Adding Kutta Condition - #
dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

# - Prepping for Least Sqaures Solve - #
# without kutta condition
LHSlsq_nokutta, RHSlsq_nokutta = dt.prep_leastsquares(
    LHSnokutta, RHSnokutta, prescribedpanels
)
# with kutta condition
LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

#---------------------------------#
#             Solving             #
#---------------------------------#

mured = LHSlsq \ RHSlsq
mured_nokutta = LHSlsq_nokutta \ RHSlsq_nokutta

# - Solving Without Kutta - #
# mu_nokutta = LHSlsq_nokutta\RHSlsq_nokutta
mu_nokutta = dt.mured2mu(mured_nokutta, prescribedpanels)

# - Solving With Kutta - #
mu = dt.mured2mu(mured, prescribedpanels)

#---------------------------------#
#         Post-Processing         #
#---------------------------------#

### --- Velocity Contributions --- ###
# - Body-induced Surface Velocity - #
Vb_nokutta = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu_nokutta)
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu)

# - "Wake"-induced Surface Velocity - #
Vb_nokutta_wake = dt.vfromTE(panels.controlpoint, panels.TEnodes, mu_nokutta)
Vb_wake = dt.vfromTE(panels.controlpoint, panels.TEnodes, mu)

# - ∇μ/2 surface velocity - #
Vb_gradmu = dt.vfromgradmu(panels, mu)

# - Total Velocity - #
V_nokutta = Vb_nokutta .+ Vb_nokutta_wake
V_nogradmu = Vb .+ Vb_wake
Vtot = Vsmat .+ Vb .+ Vb_wake .+ Vb_gradmu

# ### --- Velocity Tangent to Surface --- ###
# Vtan_nokutta = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nokutta), eachrow(panels.tangent))]
# Vtan_nogradmu = [dt.dot(v,t) for (v,t) in zip(eachrow(V_nogradmu), eachrow(panels.tangent))]
# Vtan = [dt.dot(v,t) for (v,t) in zip(eachrow(Vtot), eachrow(panels.tangent))]

### --- Steady Surface Pressure --- ###
cp_nokutta = 1.0 .- (dt.norm.(eachrow(V_nokutta)) / Vinf) .^ 2
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
ncut = 2
plot!(pp, xs[ncut:(40 - ncut)], cp[ncut:(40 - ncut)]; label="DuctAPE")

savefig(savepath * "body-pressure-comp.pdf")

pv = plot(; xlabel="x", ylabel=L"Vs/V_\infty")
plot!(
    pv,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    label="experimental",
)
xs = panels.controlpoint[:, 1]
plot!(pv, xs[41:end], dt.norm.(eachrow(Vtot[41:end, :])) ./ Vinf; label="DuctAPE")

savefig(pv, savepath * "body-vel-comp.pdf")

ps = plot(; xlabel="x", ylabel="panel strength")
plot!(ps, xs, mu; seriestype=:scatter)
savefig(ps, savepath * "body-panel-strengths.pdf")
