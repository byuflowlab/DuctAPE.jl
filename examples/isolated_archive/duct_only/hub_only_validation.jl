#=

Example of running hub/hub alone (without rotor) in hubTAPE.
Note that currently (as of April 2023) this is best done with FLOWFoil.jl, which specializes in 2D panel methods.

Authors: Judd Mehr,

=#
using DuctAPE
const dt = DuctAPE

include("../plots_default.jl")

#---------------------------------#
#        Define Coordinates       #
#---------------------------------#

# - Hub Coordinates - #
# use hub coordinates from FLOWFoil validation cases
#=
Similarly, the hub coordinates here contain x_hub and r_hub, which contain the x and r coordinates from the leading edge to the trailing edge of the center body (hub), thus the coordinates are also effectively clockwise.
=#
include("../test/data/bodyofrevolutioncoords.jl")
hub_coordinates = [x_hub r_hub]

##### ----- Panels ----- #####
body_panels = dt.generate_body_panels(nothing, hub_coordinates)

#---------------------------------#
#  Compute Coupling Coefficients  #
#---------------------------------#

##### ----- Mesh ----- #####
# mesh_bb = dt.generate_body_mesh(body_panels)
mesh_bb = dt.generate_one_way_mesh(body_panels, body_panels)

kutta_idxs = dt.get_kutta_indices([true], mesh_bb)

##### ----- Linear System ----- #####
#=
Note: the system returned by FLOWFoil already includes the Kutta condtion, which in the axisymmetric case uses a subtractive method. This means that the system is actually an N-1xN-1 system where N is the total number of panels (where we apply the Kutta condition to the hub but not the hub).
We will recover the missing hub trailing edge panel stregth in the post-processing phase
panel_strengths = dt.solve_body_only(coupling_coefficient_matrix, boundary_conditions)
=#
#returns the right-hand-side coupling coefficient matrix (A in Ax=b)
#and returns the left-hand-side boundary condition vector (b in Ax=b)
# A, b = dt.generate_body_linear_system(body_panels, mesh_bb)
A_bb = dt.assemble_induced_velocity_on_body_matrix(
    mesh_bb, body_panels, body_panels; singularity="vortex"
)

# # apply back-diagonal correction to hub portions of coefficient matrix
# dt.apply_back_diagonal_correction!(
#     A_bb, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
# )

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
vs_hub, cp_hub = dt.states_to_outputs_body_only(
    panel_strengths, body_panels, mesh_bb; Vinf=Vinf
)

#---------------------------------#
#              PLOTS              #
#---------------------------------#
#note that we already have the x-locations associated with the outputs as the panel center locations of the geometry.

# - Plot Surface Pressure Coefficients - #
plot(; xlabel="x/c", ylabel=L"V/V_\infty")
plot!(body_panels[1].panel_center[:, 1], vs_hub; label="hub")
plot!(Vs_over_Vinf_x, Vs_over_Vinf_vs; seriestype=:scatter, label="Experimental")
savefig("examples/hub-only-surface-velocity.pdf")
