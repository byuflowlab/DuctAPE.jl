using FLOWFoil
const ff = FLOWFoil
using DuctTAPE
const dt = DuctTAPE
include("../../plots_default.jl")

Vinf = 1.0
x, r = ff.naca4(0.0, 0.0, 12.0)
r .+= 100.0
coordinates = [x r]
geom = plot(x, r; aspectratio=:equal, xlabel="x", ylabel="r", label="Airfoil Surface")

method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [false])
problem = ff.define_problem(method, coordinates, 0.0, -1.0, -1.0)
panels = ff.generate_panels(method, coordinates)
mesh = ff.generate_mesh(method, panels)
system = ff.generate_inviscid_system(method, panels, mesh)
solution = ff.solve(system)

surface_velocity = solution.x[1:(end - 1)]
surface_pressure = 1.0 .- surface_velocity .^ 2

sv = plot(
    panels.panel_center[:, 1],
    surface_velocity;
    label="Panel Method",
    xlabel="x/c",
    ylabel=L"\frac{V_s}{V_\infty}",
)
savefig("test/manual_tests/0012_surface_velocity.pdf")

sp = plot(
    panels.panel_center[:, 1],
    surface_pressure;
    label="Panel Method",
    xlabel="x/c",
    ylabel=L"c_p",
    yflip=true,
)
savefig("test/manual_tests/0012_surface_pressure.pdf")

## -- Set Up Test Case 1 -- ##
minr = findfirst(x -> x < 0.5, panels.panel_center[:, 1]) #findmin(panels.panel_center[:, 2])
x1 = ones(50) * panels.panel_center[minr, 1]
r1 = range(panels.panel_center[minr, 2], panels.panel_center[minr, 2] - 15.0; length=50)
coords1 = [x1 r1]
panels1 = ff.generate_panels(method, coords1)
cpr1 = panels1.panel_center[:, 2]
mesh1 = dt.generate_one_way_mesh(panels, panels1)

plot!(geom, x1, r1; label="Line of Interest 1")

# Calculate Coefficients
unit_induced_axial_velocity1, unit_induced_radial_velocity1 = dt.assemble_induced_velocity_matrices(
    mesh1, panels, panels1
)

axial_velocity1 = Vinf .+ unit_induced_axial_velocity1 * surface_velocity
radial_velocity1 = unit_induced_radial_velocity1 * surface_velocity
Vmag1 = sqrt.(axial_velocity1 .^ 2 .+ radial_velocity1 .^ 2)

rotor_plane = plot(; xlabel="r", ylabel="Case 1 Velocity Magnitude")
plot!(rotor_plane, cpr1, Vmag1; label="Panel Method")
plot!(
    rotor_plane,
    cpr1,
    ones(length(cpr1)) * abs(surface_velocity[minr]);
    linestyle=:dash,
    label="Surface Velocity Magnitude",
)
savefig("test/manual_tests/0012_off_body_case_1.pdf")

savefig(geom, "test/manual_tests/naca0012_off_body_geometry.pdf")
