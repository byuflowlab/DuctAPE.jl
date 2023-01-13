#=
Basic Ducted Propeller Example

Author: Judd Mehr
=#

using FLOWFoil
using DuctTAPE
using FLOWMath
using PyPlot

### --- Create Duct Geometry and Run FLOWFoil

## -- Geometry

#use a symmetric naca airfoil for the hub
_, xhub, _, rhub = FLOWFoil.naca4(0.0, 0.0, 10.0; N=80, split=true)

#use a flipped naca airfoil for the duct
xduct, rduct = FLOWFoil.naca4(; N=80)
coords = [xduct rduct]
scale = 1.0
angle = 8.0
location = [0.0; 0.5]
xduct, rduct = FLOWFoil.position_coordinates(coords, scale, angle, location; flipped=true)

# rhub = 1e-8 * ones(length(rhub))
# xduct = [reverse(xhub); xhub[2:end]]
# rduct = [0.5 * ones(length(rhub)); rhub[2:end] .+ 0.5]

# generate mesh objects
duct = FLOWFoil.generate_axisym_mesh(xduct, rduct; bodyofrevolution=false)
hub = FLOWFoil.generate_axisym_mesh(xhub, rhub; bodyofrevolution=true)
meshes = [duct; hub]

## -- Set up FLOWFoil Problem
problem = FLOWFoil.Problem(meshes; axisymmetric=true, viscous=false)

## -- Solve FLOWFoil System
ff_solution = FLOWFoil.solve(problem)

####create dummy rotor for now...
struct Rotor
    radial_stations
    xlocation
end
rotor = Rotor(range(0.0, 1.0; length=100), 0.5)

figure(1; figsize=(6, 4))
clf()

plot(xduct, rduct)
plot(xhub, rhub)
axis("equal")

### --- Sample Velocity Field
# Need to get wall mesh geometry in order to define rotor blade dimensions
xduct, yduct, xhub, yhub = DuctTAPE.extract_ff_geom(ff_solution)

# duct wall geometry needs to be split before being splined
xductinner, _, yductinner, _ = DuctTAPE.split_wall(xduct, yduct)

#spline duct and hub geometries
ductspline = FLOWMath.Akima(reverse(xductinner), reverse(yductinner))
hubspline = FLOWMath.Akima(xhub, yhub)

#find Rhub from point on hub spline corresponding with rotor xlocation
rhub = hubspline(rotor.xlocation)

#find Rtip from point on duct spline corresponding with rotor xlocation
rtip = ductspline(rotor.xlocation)

#get dimensional rotor radial station locations
rotor_rdims = DuctTAPE.lintran(
    rhub, #range a start
    rtip, #range a end
    rotor.radial_stations[1], #range b start
    rotor.radial_stations[end], #range b end
    rotor.radial_stations, #range to transform
)

#assemble field points
field_points = [[rotor.xlocation; rotor_rdims[i]] for i in 1:length(rotor_rdims)]

freestream = DuctTAPE.Freestream(10.0, 1.225, 1.81e-5, 343.0)

rotor_vels = DuctTAPE.probe_ff_velocity(
    ff_solution,
    field_points,
    freestream.vinf;
    rho=1.225,
    mu=1.81e-5,
    treat_singularity=true,
)

total_vels_x = getindex.(rotor_vels, 1)

plot(
    [rotor.xlocation; rotor.xlocation],
    [rotor_rdims[1]; rotor_rdims[end]];
    linewidth=3,
    label="rotor position",
)
legend()
savefig("test_geometry.pdf"; bbox_inches="tight")

figure(2)
# clf()

plot(total_vels_x, rotor_rdims; linewidth=2, label="x induced velocity + Vinf")

plot(getindex.(rotor_vels, 2), rotor_rdims; linewidth=2, label="r induced velocity")

plot(
    [freestream.vinf; freestream.vinf],
    [rotor_rdims[1]; rotor_rdims[end]],
    "--k";
    label="Vinf",
)
xlabel("v/Vinf")
ylabel("r")
legend()
savefig("velocity_probe_test.pdf"; bbox_inches="tight")
