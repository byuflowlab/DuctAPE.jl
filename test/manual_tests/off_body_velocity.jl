using FLOWFoil
const ff = FLOWFoil
using DuctTAPE
const dt = DuctTAPE
include("../../plots_default.jl")

theta = range(180.0, 0.0; length=181)
radius = 0.5
x = radius * cosd.(theta) #.+ 0.5
r = radius * sind.(theta)
geom = plot(x, r; aspectratio=:equal, label="sphere surface")
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
    ur = -Vinf .* cos.(theta) .* (1.0 .- (sphere_radius ./ r) .^ 3)
    utheta = Vinf .* sin.(theta) .* (1.0 .+ 0.5 * (sphere_radius ./ r) .^ 3)
    return ur, utheta
end

ura, uta = flow_about_sphere(1.0, theta * pi / 180.0, radius, radius)
cpa = 1.0 .- uta .^ 2 # radial velocity is zero on surface

# - Visual Comparison - #
plot(; xlabel="x/c", ylabel="surface velocity")
plot!(xcp, surface_velocity; label="Panel Method")
plot!(x, uta; linestyle=:dash, linewidth=2, label="Analytic Solution")
savefig("test/manual_tests/sphere_surface_velocity.pdf")
plot(; xlabel="x/c", ylabel="surface pressure")
plot!(xcp, surface_pressure; label="Panel Method")
plot!(x, cpa; linestyle=:dash, linewidth=2, label="Analytic Solution")
savefig("test/manual_tests/sphere_surface_pressure.pdf")

#---------------------------------#
#       Off-body Velocities       #
#---------------------------------#

# Create points at which to sample velocity, choose some straight up from middle, and some off at an angle.
xup = zeros(50)
rup = range(0.5, 3.0; length=50)
ucoordinates = [xup rup]

xangle = range(-0.5, 3.0; length=50)
rangle = range(0.5, 3.0; length=50)
acoordinates = [xangle rangle]

# Get "panels" and "meshes"
upanels = ff.generate_panels(method, ucoordinates)
ux = upanels.panel_center[:, 1]
ur = upanels.panel_center[:, 2]
plot!(geom, ux, ur; label="line of interest 1")
umesh = dt.generate_one_way_mesh(panels, upanels)

apanels = ff.generate_panels(method, acoordinates)
ax = apanels.panel_center[:, 1]
ar = apanels.panel_center[:, 2]
plot!(geom, ax, ar; label="line of interest 2")
amesh = dt.generate_one_way_mesh(panels, apanels)

# Calculate Coefficients
ucoeff = dt.assemble_one_way_coefficient_matrix(umesh, panels, upanels)
acoeff = dt.assemble_one_way_coefficient_matrix(amesh, panels, apanels)

# Compute Velocity
vupmag = ucoeff * surface_velocity
vanglemag = acoeff * surface_velocity

# Analytic solution uses same function
# need to find r and theta from x and r definitions
arr = sqrt.(ax .^ 2 .+ ar .^ 2)
at = atan.(ar, ax)
plot!(geom, arr .* cos.(at), arr .* sin.(at); linestyle=:dash, label="double check line 2")

urup, utup = flow_about_sphere(1.0, pi / 2.0, upanels.panel_center[:, 2])
umag = sqrt.(urup .^ 2 .+ utup .^ 2)

urangle, utangle = flow_about_sphere(1.0, at, arr)
amag = sqrt.(urangle .^ 2 .+ utangle .^ 2)

# - Visually Compare - #

# xvup = plot(; xlabel="x/c", ylabel="velocity")
# plot!(xvup, ux, vupmag; label="Panel Method")
# plot!(xvup, ux, umag; linestyle=:dash, linewidth=2, label="Analytic Solution")

yvup = plot(; xlabel="r/c", ylabel="velocity")
plot!(yvup, ur, vupmag; label="Panel Method")
plot!(yvup, ur, umag; linestyle=:dash, linewidth=2, label="Analytic Solution")

# plot(xvup, yvup; layout=(1, 2))
savefig("test/manual_tests/off-body-case1.pdf")

yvangle = plot(; xlabel="r/c", ylabel="velocity")
plot!(yvangle, ar, vanglemag; label="Panel Method")
plot!(yvangle, ar, amag; linestyle=:dash, linewidth=2, label="Analytic Solution")
savefig("test/manual_tests/off-body-case2.pdf")

savefig(geom, "test/manual_tests/geometry.pdf")
