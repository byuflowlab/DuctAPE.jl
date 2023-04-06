#=

Example of running duct/hub alone (without rotor) in DuctTAPE.
Note that currently (as of April 2023) this is best done with FLOWFoil.jl, which specializes in 2D panel methods.

Authors: Judd Mehr,

=#
using DuctTAPE
const dt = DuctTAPE

include("../plots_default.jl")

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
duct_coordinates = [x_duct r_duct]

# - Hub Coordinates - #
# use hub coordinates from FLOWFoil validation cases
#=
Similarly, the hub coordinates here contain x_hub and r_hub, which contain the x and r coordinates from the leading edge to the trailing edge of the center body (hub), thus the coordinates are also effectively clockwise.
=#
include("../test/data/bodyofrevolutioncoords.jl")
hub_coordinates = [x_hub r_hub]

##### ----- Panels ----- #####
body_panels = dt.generate_body_panels(duct_coordinates, hub_coordinates)

#---------------------------------#
#  Compute Coupling Coefficients  #
#---------------------------------#

##### ----- Mesh ----- #####
body_mesh = dt.generate_body_mesh(body_panels)

##### ----- Linear System ----- #####
#=
Note: the system returned by FLOWFoil already includes the Kutta condtion, which in the axisymmetric case uses a subtractive method. This means that the system is actually an N-1xN-1 system where N is the total number of panels (where we apply the Kutta condition to the duct but not the hub).
We will recover the missing duct trailing edge panel stregth in the post-processing phase
panel_strengths = dt.solve_body_only(coupling_coefficient_matrix, boundary_conditions)
=#
#returns the right-hand-side coupling coefficient matrix (A in Ax=b)
#and returns the left-hand-side boundary condition vector (b in Ax=b)
A, b = dt.generate_body_linear_system(body_panels, body_mesh)

#---------------------------------#
#              SOLVE              #
#---------------------------------#
# This solve function uses ImplicitAD's implicit_linear function for convenience in automatic differentiation.
panel_strengths = dt.solve_body_system(A, b)

#---------------------------------#
#          Post-Process           #
#---------------------------------#
#=
For now, post-processing returns the extracted panel strenghts (synonymous with the surface velocities) and the associated pressure coefficients.
=#
vs_duct, vs_hub, cp_duct, cp_hub = dt.states_to_outputs_body_only(
    panel_strengths, body_panels, body_mesh
)

#---------------------------------#
#              PLOTS              #
#---------------------------------#
#note that we already have the x-locations associated with the outputs as the panel center locations of the geometry.

# - Plot Surface Velocity - #
plot(; xlabel="x/c", ylabel=L"v_s/V_\infty")
plot!(body_panels[1].panel_center[:, 1], vs_duct; label="Duct")
plot!(body_panels[2].panel_center[:, 1], vs_hub; label="Hub")
savefig("examples/duct-only-surface-velocity.pdf")

# - Plot Surface Pressure Coefficients - #
plot(; xlabel="x/c", ylabel=L"c_p", yflip=true)
plot!(body_panels[1].panel_center[:, 1], cp_duct; label="Duct")
plot!(body_panels[2].panel_center[:, 1], cp_hub; label="Hub")
savefig("examples/duct-only-surface-pressure.pdf")
