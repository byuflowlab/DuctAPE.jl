
using DuctTAPE
using Plots

include("../plots_default.jl")
include("../data/dfdc/dstestr2_case.jl");

"""
"""
# function initialize_geometry()
plot(; xlabel="x", ylabel="r")
# - Split Wall Coordinates
outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(ductx, ductr)

# - Plot Geometry
plot!(innerwallx, innerwallr; aspectratio=:equal)
# plot!(outerwallx, outerwallr; linestyle=:dash, color=1)
plot!(hubx, hubr; color=2)

ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
    innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
)

# --- GENERATE ROTOR OBJECT ARRAY

#generate rotor object
rotor1 = DuctTAPE.RotorGeometry(
    xdisk1, nblade1, rnondim1, 0.0, chord1, beta1, nothing, nothing, nothing, nothing, rpm
)

#generate stator object (rpm is zero for stator)
rotor2 = DuctTAPE.RotorGeometry(
    xdisk2, nblade2, rnondim2, 0.0, chord2, beta2, nothing, nothing, nothing, nothing, 0.0
)

#assemble array
rotors = [rotor1; rotor2]

num_radial_stations = 10
grid_options = DuctTAPE.defineGridOptions(num_radial_stations)

xg, rg, nx, nr = DuctTAPE.generate_grid_points(
    ductgeometry, ductsplines, rotors, grid_options
)

plot!(xg, rg; color=3, linewidth=0.5)
plot!(xg', rg'; color=3, linewidth=0.5)

# return nothing
# end

# initialize_geometry()

#function initialize_geometry(;
#    wallloc=0.0,
#    wallscale=1.0,
#    hubloc=0.0,
#    hubscale=1.0,
#    nrad=9,
#    filename="testgrid.pdf",
#    saveplot=false,
#)
#    # --- WALL GEOMETRY DEFINITION
#    # use naca 4420
#    xwall, ywall = FLOWFoil.naca4(4, 4, 20; N=20)

#    #flip the airfoil y coordinates
#    wallcoords = [xwall -ywall]

#    #rotate the airfoile
#    wallangle = 8.0

#    #translate the airfoil up
#    walllocation = [wallloc; 0.75]

#    #apply modifications to airfoil
#    FLOWFoil.position_coordinates(wallcoords, wallscale, wallangle, walllocation)

#    #extract coordinates
#    wallx = wallcoords[:, 1]
#    wallr = wallcoords[:, 2]

#    #split coordinates
#    outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(wallx, wallr)

#    # --- HUB GEOMETRY DEFINITION

#    # use naca 2410, get only top half
#    _, hubxcoordinates, _, hubrcoordinates = naca4(2, 4, 10; N=20, split=true)

#    # combine coordinates
#    hubcoordinates = [hubxcoordinates hubrcoordinates]

#    #adjust coordinates as needed
#    position_coordinates(hubcoordinates, hubscale, 0.0, [hubloc; 0.0])

#    #extract coordinates
#    hubx = hubcoordinates[:, 1]
#    hubr = hubcoordinates[:, 2]

#    # --- GRID POINTS DEFINITION
#    # define rotor radial stations
#    radstash = collect(range(0.0, 1.0; length=nrad))

#    #define rotor objects
#    rotors = [
#        DuctTAPE.Rotor(
#            0.25,
#            5,
#            radstash,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#        )
#        DuctTAPE.Rotor(
#            0.5,
#            7,
#            radstash,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#            nothing,
#        )
#    ]

#    #define grid options
#    grid_options = DuctTAPE.defineGridOptions(nrad)

#    #define geometry objects
#    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
#        innerwallx,
#        innerwallr,
#        outerwallx,
#        outerwallr,
#        hubcoordinates[:, 1],
#        hubcoordinates[:, 2],
#    )

#    # Generate initialized grid points
#    grid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

#    #extract grid points
#    xg = grid.x_grid_points
#    rg = grid.r_grid_points
#    nx = grid.nx
#    nr = grid.nr

#    # PLOT GEOMETRY
#    figure(1; figsize=(10, 3))
#    clf()

#    plot(
#        ductgeometry.wallinnerxcoordinates,
#        ductgeometry.wallinnerrcoordinates,
#        "C0";
#        linewidth=2,
#        label="inner duct wall",
#        zorder=2,
#    )
#    plot(
#        ductgeometry.wallouterxcoordinates,
#        ductgeometry.wallouterrcoordinates,
#        "--C0";
#        linewidth=2,
#        label="outer duct wall",
#    )

#    plot(
#        ductgeometry.hubxcoordinates,
#        ductgeometry.hubrcoordinates,
#        "C1";
#        linewidth=2,
#        label="hub",
#        zorder=2,
#    )

#    plot(xg, rg, "C2"; zorder=1)
#    plot(xg', rg', "C2"; zorder=1)
#    plot(xg[:, 1] .* NaN, rg[:, 1] .* NaN, "C2"; zorder=1, label="Initial Grid")

#    axis("equal")

#    legend(; loc=1)

#    axis("off")

#    if saveplot
#        savefig("./examples/"*filename; bbox_inches="tight")
#    end

#    return ductgeometry, ductsplines, rotors, grid
#end

#initialize_geometry()
