#=

Example of running duct/hub alone (without rotor) in DuctTAPE.
Note that currently (as of April 2023) this is best done with FLOWFoil.jl, which specializes in 2D panel methods.

Authors: Judd Mehr,

=#
using DuctTAPE
const dt = DuctTAPE
using FLOWFoil
const ff = FLOWFoil

include("../plots_default.jl")

######################################################################
#                                                                    #
#                        DuctTAPE Analysis                           #
#                                                                    #
######################################################################
#---------------------------------#
#        Define Coordinates       #
#---------------------------------#

# - Duct Coordinates - #
# use duct coordinates from FLOWFoil validation cases
#=
The file containing the duct coordinates contains the following geometry items:
- x_duct : x-coordinates of duct geometry defined from trailing edge to trailing edge clockwise
- r_duct : r-coordinates of duct geometry defined from trailing edge to trailing edge clockwise
note in this case, we include the radial offset of the duct since we have no rotor tip radius to define the duct radial location.
=#
include("../test/data/naca_662-015.jl")
duct_coordinates = [x_duct r_duct .+ 10.0]

##### ----- Panels ----- #####
body_panels = dt.generate_body_panels(duct_coordinates, nothing)

#---------------------------------#
#  Compute Coupling Coefficients  #
#---------------------------------#

##### ----- Mesh ----- #####
# mesh_bb = dt.generate_body_mesh(body_panels)
mesh_bb = dt.generate_one_way_mesh(body_panels, body_panels)

kutta_idxs = dt.get_kutta_indices([false], mesh_bb)

##### ----- Linear System ----- #####
#=
Note: the system returned by FLOWFoil already includes the Kutta condtion, which in the axisymmetric case uses a subtractive method. This means that the system is actually an N-1xN-1 system where N is the total number of panels (where we apply the Kutta condition to the duct but not the hub).
We will recover the missing duct trailing edge panel stregth in the post-processing phase
panel_strengths = dt.solve_body_only(coupling_coefficient_matrix, boundary_conditions)
=#
#returns the right-hand-side coupling coefficient matrix (A in Ax=b)
#and returns the left-hand-side boundary condition vector (b in Ax=b)
# A, b = dt.generate_body_linear_system(body_panels, mesh_bb)
A_bb = dt.assemble_induced_velocity_on_body_matrix(
    mesh_bb, body_panels, body_panels; singularity="vortex"
)

# apply back-diagonal correction to duct portions of coefficient matrix
dt.apply_back_diagonal_correction!(
    A_bb, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
)

# - freestream to body - #
# b_bfff = ff.assemble_ring_boundary_conditions_raw(
#     ff.Dirichlet(), [false; true], body_panels, body_meshff
# )
Vinf = 5.0
b_bf = Vinf .* dt.assemble_body_freestream_boundary_conditions(body_panels, mesh_bb)

#---------------------------------#
#              SOLVE              #
#---------------------------------#
# This solve function uses ImplicitAD's implicit_linear function for convenience in automatic differentiation.
# panel_strengths = dt.solve_body_system(A, b)
panel_strengths = dt.solve_body_system(A_bb, b_bf, kutta_idxs) # get circulation strengths from solving body to body problem

#---------------------------------#
#          Post-Process           #
#---------------------------------#
#=
For now, post-processing returns the extracted panel strenghts (synonymous with the surface velocities) and the associated pressure coefficients.
=#
vs_duct, cp_duct = dt.states_to_outputs_body_only(
    panel_strengths, body_panels, mesh_bb; Vinf=Vinf
)

#---------------------------------#
#              PLOTS              #
#---------------------------------#
#note that we already have the x-locations associated with the outputs as the panel center locations of the geometry.

# - Plot Surface Pressure Coefficients - #
plot(; xlabel="x/c", ylabel=L"c_p", yflip=true)
plot!(body_panels[1].panel_center[:, 1], cp_duct; label="DuctTAPE")
# plot!(pressurexupper, pressureupper; seriestype=:scatter, label="Experimental Outer")
# plot!(pressurexlower, pressurelower; seriestype=:scatter, label="Experimental Inner")

######################################################################
#                                                                    #
#                        FLOWFoil Analysis                           #
#                                                                    #
######################################################################

method = ff.AxisymmetricProblem(ff.Vortex(ff.Constant()), ff.Dirichlet(), [false])
flow_angle = [0.0]
reynolds = [-1.0]
mach = [-1.0]
coordinates = duct_coordinates

# Generate Problem Object
problem = ff.define_problem(method, coordinates, flow_angle, reynolds, mach)

# Generate Panel Geometry
panels = ff.generate_panels(method, coordinates)

# Generate Influence Mesh
mesh = ff.generate_mesh(method, panels)

# Assemble Linear System
system = ff.generate_inviscid_system(method, panels, mesh)

# Solve Linear System
solution = ff.solve(system)

# Post Process Solution
polar = ff.post_process(method, problem, panels, mesh, solution)

# panel method solution
plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_pressure[1, :, 1];
    linestyle=:dash,
    linewidth=2,
    label="FLOWFoil Axisym",
)


ppolar = FLOWFoil.solve(duct_coordinates, flow_angle, ff.PlanarProblem(ff.Vortex(ff.Linear()), ff.Dirichlet()))
# planar method solution
plot!(
    ppolar.xsmooth[1, :, 1],
    ppolar.surface_pressure[1, :, 1];
    linestyle=:dash,
    linewidth=2,
    label="FLOWFoil Planar",
)



savefig("examples/duct-shift-surface-pressure.pdf")