# using FLOWFoil
# const ff = FLOWFoil
# using DuctTAPE
# const dt = DuctTAPE
# include("../../plots_default.jl")

# Vinf = 1.0
# prescribed_strength = 2.0
dist = [
    1e-10
    1e-9
    1e-8
    1e-7
    1e-6
    1e-5
    1e-4
    1e-3
    1e-2
    1e-1
    0.2
    0.3
    0.4
    0.5
]
N = length(dist)

# x = [-0.5; 0.5]
# r = [1.0; 1.0]
# coordinates = [x r]

# method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
# problem = ff.define_problem(method, coordinates, 0.0, -1.0, -1.0)
# panels = ff.generate_panels(method, coordinates)

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

method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
problem = ff.define_problem(method, coordinates, 0.0, -1.0, -1.0)
panels = ff.generate_panels(method, coordinates)
mesh = ff.generate_mesh(method, panels)
system = ff.generate_inviscid_system(method, panels, mesh)
solution = ff.solve(system)

paneln = 45

test_x = [-0.1; 0.1] .+ panels.panel_center[paneln, 1]
test_r = [0.0; 0.0] .+ panels.panel_center[paneln, 2]
# test_x = [0.0; 0.0] .+ panels.panel_center[paneln, 1]
# test_r = [-0.1; 0.1] .+ panels.panel_center[paneln, 2]
vx = zeros(N)
vr = zeros(N)
for i in 1:N
    # test_x .+= 10.0^(-i)
    test_r .+= dist[i]

    test_coords = [test_x test_r]
    test_panel = ff.generate_panels(method, test_coords)
    # println("panel length: ", test_panel.panel_length[1])
    test_mesh = dt.generate_one_way_mesh(panels, test_panel)
    println("mesh x: ", test_mesh.x[paneln])
    println("mesh r: ", test_mesh.r[paneln])

    unit_axial_vel, unit_radial_vel = dt.assemble_induced_velocity_matrices(
        test_mesh, panels, test_panel
    )

    vx[i] = (Vinf .+ unit_axial_vel * solution.x)[1]
    vr[i] = (unit_radial_vel * solution.x)[1]
    # test_x .-= 10.0^(-i)
    test_r .-= dist[i]
end

plot(
    dist,
    vx;
    xaxis=:log,
    # yaxis=:log,
    xlabel=L"Distance: ",
    ylabel="axial induced velocity + freestream",
)
savefig("test/manual_tests/axial_very_close.pdf")
plot(dist, vr; xaxis=:log, xlabel=L"Distance: ", ylabel="radial induced velocity")
savefig("test/manual_tests/radial_very_close.pdf")
