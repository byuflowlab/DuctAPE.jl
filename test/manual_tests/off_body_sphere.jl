using FLOWFoil
const ff = FLOWFoil
using DuctTAPE
const dt = DuctTAPE
include("../../plots_default.jl")

Vinf = 1.0
theta = range(180.0, 0.0; length=181)
radius = 0.5
x = radius * cosd.(theta) #.+ 0.5
r = radius * sind.(theta)
geom = plot(x, r; aspectratio=:equal, xlabel="x", ylabel="r", label="sphere surface")
coordinates = [x r]

method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [true])
problem = ff.define_problem(method, coordinates, 0.0, -1.0, -1.0)
panels = ff.generate_panels(method, coordinates)
mesh = ff.generate_mesh(method, panels)
system = ff.generate_inviscid_system(method, panels, mesh)
solution = ff.solve(system)

#---------------------------------#
#  Sanity Check: Surface Values   #
#---------------------------------#

# - Get Panel Solution - #
surface_velocity = solution.x
surface_pressure = 1.0 .- surface_velocity .^ 2
xcp = panels.panel_center[:, 1]

# - Analytic Solution - #
function flow_about_sphere(Vinf, theta, r, sphere_radius=0.5)
    cpr1 = -Vinf .* cos.(theta) .* (1.0 .- (sphere_radius ./ r) .^ 3)
    utheta = Vinf .* sin.(theta) .* (1.0 .+ 0.5 * (sphere_radius ./ r) .^ 3)
    return cpr1, utheta
end

analytic_radial_surface_velocity, analytic_tangential_surface_velocity = flow_about_sphere(
    Vinf, theta * pi / 180.0, radius, radius
)
cpa = 1.0 .- analytic_tangential_surface_velocity .^ 2 # radial velocity is zero on surface

# - Visual Comparison - #
plot(; xlabel="x/c", ylabel="surface velocity")
plot!(xcp, surface_velocity; label="Panel Method")
plot!(
    x,
    analytic_tangential_surface_velocity;
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)
savefig("test/manual_tests/sphere_surface_velocity.pdf")
plot(; xlabel="x/c", ylabel="surface pressure")
plot!(xcp, surface_pressure; label="Panel Method")
plot!(x, cpa; linestyle=:dash, linewidth=2, label="Analytic Solution")
savefig("test/manual_tests/sphere_surface_pressure.pdf")

#---------------------------------#
#       Off-body Velocities       #
#---------------------------------#

# Create points at which to sample velocity, choose some straight up from middle, and some off at an angle.
x1 = zeros(50)
r1 = range(0.5, 3.0; length=50)
coords1 = [x1 r1]

x2 = range(-0.5, 3.0; length=50)
r2 = range(0.5, 3.0; length=50)
coords2 = [x2 r2]

# Get "panels" and "meshes"
panels1 = ff.generate_panels(method, coords1)
cpx1 = panels1.panel_center[:, 1]
cpr1 = panels1.panel_center[:, 2]
plot!(geom, cpx1, cpr1; label="line of interest 1")
mesh1 = dt.generate_one_way_mesh(panels, panels1)

panels2 = ff.generate_panels(method, coords2)
cpx2 = panels2.panel_center[:, 1]
cpr2 = panels2.panel_center[:, 2]
plot!(geom, cpx2, cpr2; label="line of interest 2")
mesh2 = dt.generate_one_way_mesh(panels, panels2)

# Calculate Coefficients
unit_induced_axial_velocity1, unit_induced_radial_velocity1 = dt.assemble_induced_velocity_matrices(
    mesh1, panels, panels1
)

unit_induced_axial_velocity2, unit_induced_radial_velocity2 = dt.assemble_induced_velocity_matrices(
    mesh2, panels, panels2
)

# Compute Velocity

axial_velocity1 = Vinf .+ unit_induced_axial_velocity1 * surface_velocity
radial_velocity1 = unit_induced_radial_velocity1 * surface_velocity
Vmag1 = sqrt.(axial_velocity1 .^ 2 .+ radial_velocity1 .^ 2)

axial_velocity2 = Vinf .+ unit_induced_axial_velocity2 * surface_velocity
radial_velocity2 = unit_induced_radial_velocity2 * surface_velocity
Vmag2 = sqrt.(axial_velocity2 .^ 2 .+ radial_velocity2 .^ 2)

# Analytic solution uses same function
# need to find r and theta from x and r definitions
arr = sqrt.(cpx2 .^ 2 .+ cpr2 .^ 2)
at = atan.(cpr2, cpx2)
plot!(geom, arr .* cos.(at), arr .* sin.(at); linestyle=:dash, label="double check line 2")

ur1, ut1 = flow_about_sphere(1.0, pi / 2.0, panels1.panel_center[:, 2])
analytic_vmag1 = sqrt.(ur1 .^ 2 .+ ut1 .^ 2)

ur2, ut2 = flow_about_sphere(1.0, at, arr)
analytic_vmag2 = sqrt.(ur2 .^ 2 .+ ut2 .^ 2)

# - Visually Compare - #

# xvup = plot(; xlabel="x/c", ylabel="velocity")
# plot!(xvup, cpx1, Vmag1; label="Panel Method")
# plot!(xvup, cpx1, analytic_vmag1; linestyle=:dash, linewidth=2, label="Analytic Solution")

yvup = plot(; xlabel="r/c", ylabel="Case 1 Velocity Magnitude")
plot!(yvup, cpr1, Vmag1; label="Panel Method")

plot!(yvup, cpr1, analytic_vmag1; linestyle=:dash, linewidth=2, label="Analytic Solution")

# plot(xvup, yvup; layout=(1, 2))
savefig("test/manual_tests/off-body-case1.pdf")

yvangle = plot(; xlabel="r/c", ylabel="Case 2 Velocity Magnitude")
plot!(yvangle, cpr2, Vmag2; label="Panel Method")
plot!(
    yvangle, cpr2, analytic_vmag2; linestyle=:dash, linewidth=2, label="Analytic Solution"
)
savefig("test/manual_tests/off-body-case2.pdf")

savefig(geom, "test/manual_tests/geometry.pdf")
