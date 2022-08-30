"""
"""
function setup_geometry(; plotgeometry=false)

    # - Split Wall Coordinates
    outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(ductx, ductr)

    # -- DEFINE DUCT OBJECT
    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    if plotgeometry
        # PLOT GEOMETRY
        plot(; aspectratio=:equal, xlabel="x", ylabel="r", legend=true, label="")

        plot!(
            ductgeometry.wallinnerxcoordinates,
            ductgeometry.wallinnerrcoordinates;
            color=1,
            linewidth=2,
            label="inner duct wall",
        )
        plot!(
            ductgeometry.wallouterxcoordinates,
            ductgeometry.wallouterrcoordinates;
            color=1,
            linestyle=:dot,
            linewidth=2,
            label="outer duct wall",
        )

        plot!(
            ductgeometry.hubxcoordinates,
            ductgeometry.hubrcoordinates;
            color=2,
            linewidth=2,
            label="hub",
        )

        savefig("./examples/test_geometry.pdf")
    end

    return ductgeometry, ductsplines
end

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
