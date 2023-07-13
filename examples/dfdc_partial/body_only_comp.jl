project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctTAPE
const dt = DuctTAPE

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")
pyplot()

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 30.0))
pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-1.25, 1.0), yflip=true)

#---------------------------------#
#               DFDC              #
#---------------------------------#
include(datapath * "DFDC_VELOCITY_BREAKDOWN.jl")
include(datapath * "DFDC_SURFACE_VELOCITIES.jl")
include(datapath * "DFDC_ELEMENT_IDS.jl")

hubid = elids[1, 2]:elids[1, 3]
ductid = elids[2, 2]:elids[2, 3]
hwid = elids[4, 2]:elids[4, 3]
dwid = elids[14, 2]:elids[14, 3]

hubx = dfdc_velocities[hubid, 10]
hubr = dfdc_velocities[hubid, 11]
hubvs = dfdc_vs[hubid, 2]
ductx = dfdc_velocities[ductid, 10]
ductr = dfdc_velocities[ductid, 11]
ductvs = dfdc_vs[ductid, 2]
# hubwakex = dfdc_velocities[hwid, 10]
# hubwakevs = dfdc_vs[hwid, 2]
ductwakex = dfdc_velocities[dwid, 10]
ductwakevs = dfdc_vs[dwid, 2]

Vinf = 20.0
hubcp = 1.0 .- (hubvs ./ Vinf) .^ 2
ductcp = 1.0 .- (ductvs ./ Vinf) .^ 2
# hubwakecp = 1.0 .- (hubwakevs ./ Vinf) .^ 2
# ductwakecp = 1.0 .- (ductwakevs ./ Vinf) .^ 2

plot!(pv, hubx, hubvs; linestyle=:dash, label="DFDC Hub")
plot!(pv, ductx, ductvs; linestyle=:dash, label="DFDC Duct")
# plot!(pv, hubwakex, hubwakevs; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pv, ductwakex, ductwakevs; linestyle=:dash, label="DFDC Duct Wake")

plot!(pc, hubx, hubcp; linestyle=:dash, label="DFDC Hub")
plot!(pc, ductx, ductcp; linestyle=:dash, label="DFDC Duct")
# plot!(pc, hubwakex, hubwakecp; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pc, ductwakex, ductwakecp; linestyle=:dash, label="DFDC Duct Wake")

#---------------------------------#
#             DuctTAPE            #
#---------------------------------#
# include(datapath * "lewis_refined_duct.jl")
# include(datapath * "lewis_refined_hub.jl")

# - Generate Panels - #
# use DFDC repaneling for comparison
duct_coordinates = reverse([ductx ductr]; dims=1)
hub_coordinates = reverse([hubx hubr]; dims=1)

write_coordinates(
    duct_coordinates;
    filename="dfdc_lewis_duct.jl",
    savepath=savepath,
    varname="duct_coordinates",
)
write_coordinates(
    hub_coordinates;
    filename="dfdc_lewis_hub.jl",
    savepath=savepath,
    varname="hub_coordinates",
)

panels = dt.generate_panels([duct_coordinates, hub_coordinates])
# panels = dt.generate_panels([hub_coordinates])
# panels = dt.generate_panels([duct_coordinates])

##### ----- Visualize to Check ----- #####
visualize_paneling(;
    body_panels=panels,
    coordinates=[duct_coordinates, hub_coordinates],
    controlpoints=true,
    nodes=true,
    normals=true,
    normal_scaling=0.1,
    savepath=savepath,
    filename="body-geometry.pdf",
    legendloc=:right,
)

# - Operating Conditions - #
# Define freestream on panels
Vinf = 20.0 #magnitude doesn't matter yet.
Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
Vsmat = repeat(Vs, panels.npanels) # need velocity on each panel

prescribedpanels = [(1, 0.0); (panels.npanels, 0.0)]
# prescribedpanels = [(1, 0.0)]
# prescribedpanels = []

## -- Induced Velocities -- ##

# - Initial System Matrices - #
LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
RHS = dt.freestream_influence_vector(panels.normal, Vsmat)

# - Adding Kutta Condition - #
dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

# - Prepping for Least Sqaures Solve - #
LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

# - Solving - #
mured = LHSlsq \ RHSlsq
mu = dt.mured2mu(mured, prescribedpanels)
# mu = LHS \ RHS

# - Post-processing - #
# - Body-induced Surface Velocity - #
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mu)

# - "Wake"-induced Surface Velocity - #
# Vb_wake = dt.vfromTE(panels.controlpoint, panels.TEnodes, mu)

# - ∇μ/2 surface velocity - #
Vb_gradmu = dt.vfromgradmu(panels, mu)

# - Total Velocity - #
# Vtot = Vsmat .+ Vb .+ Vb_wake .+ Vb_gradmu
Vtot = Vsmat .+ Vb .+ Vb_gradmu
Vs = dt.norm.(eachrow(Vtot))
### --- Steady Surface Pressure --- ###
bodycp = 1.0 .- (Vs / Vinf) .^ 2

# - Plotting - #
xs = panels.controlpoint[:, 1]
ncut = 2
plot!(pv, xs, Vs; label="DuctTAPE")
plot!(pc, xs, bodycp; label="DuctTAPE")

#---------------------------------#
#             Exp data            #
#---------------------------------#
include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
plot!(
    pc,
    pressurexupper,
    pressureupper;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    markersize=3,
    label="exp outer",
)
plot!(
    pc,
    pressurexlower,
    pressurelower;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:dtriangle,
    markersize=3,
    label="exp inner",
)

plot!(
    pv,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs * Vinf;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    markersize=3,
    label="experimental",
)
###########################################
savefig(pv, savepath * "velocity-comp.pdf")
savefig(pc, savepath * "pressure-comp.pdf")
