project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")
# pyplot()
gr()

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 30.0))
pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-1.25, 1.0), yflip=true)

#---------------------------------#
#               DFDC              #
#---------------------------------#
# include(datapath * "DFDC_VELOCITY_BREAKDOWN.jl")
# include(datapath * "DFDC_SURFACE_VELOCITIES.jl")
# include(datapath * "DFDC_ELEMENT_IDS.jl")
include(datapath * "DFDC_CPS.jl")

# hubid = elids[1, 2]:elids[1, 3]
# ductid = elids[2, 2]:elids[2, 3]
# hwid = elids[4, 2]:elids[4, 3]
# dwid = elids[14, 2]:elids[14, 3]

hubx = dfdc_hub_cp[:, 1]
hubr = dfdc_hub_cp[:, 2]
hubvs = dfdc_hub_cp[:, end - 1]
hubcp = dfdc_hub_cp[:, 4]

ductx = dfdc_duct_cp[:, 1]
ductr = dfdc_duct_cp[:, 2]
ductvs = dfdc_duct_cp[:, end - 1]
ductcp = dfdc_duct_cp[:, 4]

ductwakex = dfdc_wake14_cp[:, 1]
ductwaker = dfdc_wake14_cp[:, 2]
ductwakevs = dfdc_wake14_cp[:, end - 1]
ductwakecp = dfdc_wake14_cp[:, 4]

plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
# plot!(pv, hubwakex, hubwakevs; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pv, ductwakex, ductwakevs; linestyle=:dash, label="DFDC Duct Wake")

plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
# plot!(pc, hubwakex, hubwakecp; linestyle=:dash, label="DFDC Hub Wake")
# plot!(pc, ductwakex, ductwakecp; linestyle=:dash, label="DFDC Duct Wake")

dfdc_duct_coordinates = reverse([ductx ductr]; dims=1)
dfdc_hub_coordinates = reverse([hubx hubr]; dims=1)
dfdc_panels = dt.generate_panels([dfdc_duct_coordinates, dfdc_hub_coordinates])

#---------------------------------#
#             DuctAPE            #
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

# panels = dt.generate_panels([duct_coordinates, hub_coordinates])
# panels = dt.generate_panels([hub_coordinates])
panels = dt.generate_panels([duct_coordinates])


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
Vsmat = repeat(Vs, panels.totpanel) # need velocity on each panel

_, leid = findmin(duct_coordinates[:, 1])
# prescribedpanels = [(1, 0.0); (panels.totpanel, 0.0)]
prescribedpanels = [(1, 0.0)]
# prescribedpanels = [(leid, 0.0)]
# prescribedpanels = []

## -- Induced Velocities -- ##

# - Initial System Matrices - #
A_bb = dt.doublet_panel_influence_matrix(panels.nodes, panels)
# - add internal panel stuff - #
LHS = dt.doublet_panel_influence_on_internal_panels(A_bb, panels, panels)
b_bf = dt.freestream_influence_vector(panels.normal, Vsmat)
# - add internal panel stuff - #
RHS = dt.freestream_influence_on_internal_panels(b_bf, panels, Vs)

# - Adding Kutta Condition - #
dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

mu = dt.solve_body_strengths(LHS,RHS, prescribedpanels, panels.nbodies)

# # - Prepping for Least Sqaures Solve - #
# LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)
# # - Solving - #
# mured = LHSlsq \ RHSlsq
# mu = dt.mured2mu(mured, prescribedpanels, panels.nbodies)
# # mu = LHS \ RHS

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
plot!(pv, xs, Vs; label="DuctAPE")
plot!(pc, xs, bodycp; label="DuctAPE")

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
    label="exp duct outer",
)
plot!(
    pc,
    pressurexlower,
    pressurelower;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:dtriangle,
    markersize=3,
    label="exp duct inner",
)

plot!(
    pv,
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs * Vinf;
    seriestype=:scatter,
    color=myblue[1],
    markershape=:utriangle,
    markersize=3,
    label="hub experimental",
)
###########################################
savefig(pv, savepath * "velocity-comp.png")
savefig(pc, savepath * "pressure-comp.png")
savefig(pv, savepath * "velocity-comp.pdf")
savefig(pc, savepath * "pressure-comp.pdf")
