using FLOWFoil
const ff = FLOWFoil
using DuctAPE
const dt = DuctAPE

include("../../plots_default.jl")

# Duct Geometry
include("../data/marin_19a_duct.jl")
include("../data/marin_19a_duct_alone_surface_pressure.jl")

plot(; aspectratio=:equal)
plot!(full_x, full_r; markershape=:circle, markersize=2, label="Raw Geometry")

# - FLOWFoil Solution for Duct Alone - #

coordinates = ff.repanel_airfoil([full_x full_r]; N=100)

plot(; aspectratio=:equal)
plot!(
    coordinates[:, 1],
    coordinates[:, 2];
    seriestype=:scatter,
    markershape=:circle,
    markersize=1,
    label="Smooth Duct Geometry",
)
savefig("duct_geometry.pdf")

# - Check that geometry yields smooth results - #

method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false])
problem = ff.define_problem(method, coordinates, 0.0, -1.0, -1.0)
body_panels = ff.generate_panels(method, coordinates)
mesh_body_to_body = ff.generate_mesh(method, body_panels)

body_system = ff.generate_inviscid_system(method, body_panels, mesh_body_to_body)

gammas = body_system.A \ body_system.b
surface_velocity = gammas[1:(end - 1)]
surface_pressure = 1.0 .- surface_velocity .^ 2

plot(
    body_panels.panel_center[:, 1],
    surface_velocity;
    xlabel="x/c",
    ylabel=L"Vs/V_\infty",
    label="",
)
savefig("marin_19a_duct_surface_velocity.pdf")

plot(
    body_panels.panel_center[:, 1],
    surface_pressure;
    xlabel="x/c",
    ylabel=L"c_p",
    yflip=true,
    label="flowfoil",
)
plot!(inner_cp_vs_x[:, 1], inner_cp_vs_x[:, 2]; seriestype=:scatter, label="exp inner")
plot!(outer_cp_vs_x[:, 1], outer_cp_vs_x[:, 2]; seriestype=:scatter, label="exp outer")
savefig("marin_19a_duct_surface_pressure.pdf")
