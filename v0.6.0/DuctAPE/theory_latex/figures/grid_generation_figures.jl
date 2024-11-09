#=
Figures for Duct Grid Generation section of dissertation
=#

include("../../plots_default.jl")

using FLOWFoil
using DuctTAPE

# Get duct wall geometry
xwall, ywall = naca4(4, 4, 20)
wallcoords = [xwall -ywall]
wallangle = 10.0
walllocation = [0.0; 0.75]
FLOWFoil.position_coordinates(wallcoords, 1.0, wallangle, walllocation)
wallx = wallcoords[:, 1]
wallr = wallcoords[:, 2]

# Get Duct Hub Geometry
_, hubxcoordinates, _, hubrcoordinates = naca4(2, 4, 20; split=true)
hubcoordinates = [hubxcoordinates hubrcoordinates]
position_coordinates(hubcoordinates, 0.67, 0.0, [0.25; 0.0])
hubx = hubcoordinates[:, 1]
hubr = hubcoordinates[:, 2]

# Get Boundaries
xinlet = -0.5 * maximum(wallx)
xoutlet = 2.0 * maximum(wallx)
rmax = 1.5 * maximum(wallr)

#set up plot axis
p = plot(; size=(500, 300), ticks=false, showaxis=false, aspect_ratio=:equal)
#plot wall
plot!(wallx, wallr)
annotate!(0.5, 0.825, text("Duct Wall", 8; color=mycolors[1]))
#plot hub
plot!(hubx, hubr)
plot!(hubx, zeros(length(hubx)); linestyle=:dot, linecolor=2)
annotate!(0.5, 0.15, text("Hub", 8; color=mycolors[2]))
#plot inlet bound
plot!([xinlet; xinlet], [0.0; rmax]; linestyle=:dash, linecolor=3)
annotate!(xinlet + 0.01, rmax / 2.0, text("Inlet", 8, :left; color=mycolors[3]))
#plot outlet bound
plot!([xoutlet; xoutlet], [0.0; rmax]; linestyle=:dash, linecolor=3)
annotate!(xoutlet - 0.01, rmax / 2.0, text("Outlet", 8, :right; color=mycolors[3]))
#plot rmax bound
plot!([xinlet; xoutlet], [rmax; rmax]; linestyle=:dash, linecolor=3)
annotate!(
    xoutlet * 3 / 4,
    rmax - 0.1,
    text("Maximum Radial Boundary", 8, :right; color=mycolors[3]),
)
#plot centerline bounds
plot!([xinlet; minimum(hubx)], [0.0; 0.0]; linestyle=:dash, linecolor=3)
plot!([maximum(hubx); xoutlet], [0.0; 0.0]; linestyle=:dash, linecolor=3)
annotate!(xoutlet * 7 / 8, 0.1, text("Axis of Rotation", 8, :right; color=mycolors[3]))

savetightplot(p, "ductboundaries.tikz")

# function plot_test_grid(
#     wallloc=0.0, wallscale=1.0, hubloc=0.0, hubscale=1.0, filename="testgrid.pdf"
# )
#     # --- WALL GEOMETRY DEFINITION
#     xwall, ywall = naca4(4, 4, 20)
#     wallcoords = [xwall -ywall]
#     wallangle = 8.0
#     walllocation = [wallloc; 0.75]
#     FLOWFoil.position_coordinates(wallcoords, wallscale, wallangle, walllocation)
#     wallx = wallcoords[:, 1]
#     wallr = wallcoords[:, 2]
#     outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(wallx, wallr)

#     # --- HUB GEOMETRY DEFINITION

#     _, hubxcoordinates, _, hubrcoordinates = naca4(2, 4, 10; split=true)
#     hubcoordinates = [hubxcoordinates hubrcoordinates]
#     position_coordinates(hubcoordinates, hubscale, 0.0, [hubloc; 0.0])
#     hubx = hubcoordinates[:, 1]
#     hubr = hubcoordinates[:, 2]

#     # --- GRID POINTS DEFINITION
#     rotors = []
#     grid_options = DuctTAPE.GridOptions(15, 35, 35, 35)
#     duct = DuctTAPE.Duct(
#         innerwallx,
#         innerwallr,
#         outerwallx,
#         outerwallr,
#         hubcoordinates[:, 1],
#         hubcoordinates[:, 2],
#     )

#     x_grid_points, r_grid_points = DuctTAPE.generate_grid_points(
#         duct, rotors, grid_options; debug=false
#     )

#     # PLOTTING
#     figure(1; figsize=(10, 3))
#     clf()

#     plot(innerwallx, innerwallr, "C0"; linewidth=2, label="inner duct wall", zorder=2)
#     plot(outerwallx, outerwallr, "--C0"; linewidth=2, label="outer duct wall")

#     plot(hubx, hubr, "C1"; linewidth=2, label="hub", zorder=2)

#     plot(x_grid_points, r_grid_points, "C2"; zorder=1)
#     plot(x_grid_points', r_grid_points', "C2"; zorder=1)

#     axis("equal")

#     legend(; loc=1)

#     savefig(filename; bbox_inches="tight")

#     return nothing
# end
