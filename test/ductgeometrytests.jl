@testset "Grid Geometry" begin

    # probably want to put in basic unit tests to make sure that the grid geometry is being added as it should.

end

@testset "Grid Relaxation Stability" begin
    wallloc = 0.0
    wallscale = 1.0
    hubloc = 0.0
    hubscale = 1.0
    filename = "testgrid.pdf"

    # --- WALL GEOMETRY DEFINITION
    xwall, ywall = FLOWFoil.naca4(4, 4, 20)
    wallcoords = [xwall -ywall]
    wallangle = 8.0
    walllocation = [wallloc; 0.75]
    FLOWFoil.position_coordinates(wallcoords, wallscale, wallangle, walllocation)
    wallx = wallcoords[:, 1]
    wallr = wallcoords[:, 2]
    outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(wallx, wallr)

    # --- HUB GEOMETRY DEFINITION

    _, hubxcoordinates, _, hubrcoordinates = FLOWFoil.naca4(2, 4, 10; split=true)
    hubcoordinates = [hubxcoordinates hubrcoordinates]
    position_coordinates(hubcoordinates, hubscale, 0.0, [hubloc; 0.0])
    hubx = hubcoordinates[:, 1]
    hubr = hubcoordinates[:, 2]

    # --- GRID POINTS DEFINITION
    grid_options = DuctTAPE.GridOptions(10)
    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        hubcoordinates[:, 1],
        hubcoordinates[:, 2],
    )

    #Test that nothing breaks no matter where the rotor is placed relative to the duct
    rotors = [DuctTAPE.Rotor(0.0, nothing, nothing, nothing, nothing, nothing)]
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)
    rotors = [DuctTAPE.Rotor(0.25, nothing, nothing, nothing, nothing, nothing)]
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)
    rotors = [DuctTAPE.Rotor(0.5, nothing, nothing, nothing, nothing, nothing)]
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)
    rotors = [DuctTAPE.Rotor(0.75, nothing, nothing, nothing, nothing, nothing)]
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)
    rotors = [DuctTAPE.Rotor(1.0, nothing, nothing, nothing, nothing, nothing)]
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    @test true == true

    #Test that nothing breaks depending on how big the duct and hub are relative to each other

    rotors = [DuctTAPE.Rotor(0.25, nothing, nothing, nothing, nothing, nothing)]
    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        1.5 .* hubcoordinates[:, 1],
        1.5 .* hubcoordinates[:, 2],
    )
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        0.5 .* hubcoordinates[:, 1],
        0.5 .* hubcoordinates[:, 2],
    )
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    @test true == true

    #Test that nothing breaks when duct or hub are shifted back and forward

    rotors = [DuctTAPE.Rotor(0.5, nothing, nothing, nothing, nothing, nothing)]
    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        0.25 .+ hubcoordinates[:, 1],
        hubcoordinates[:, 2],
    )
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    rotors = [DuctTAPE.Rotor(0.25, nothing, nothing, nothing, nothing, nothing)]
    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        0.5 .+ hubcoordinates[:, 1],
        hubcoordinates[:, 2],
    )
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    duct = DuctTAPE.defineDuctGeometry(
        innerwallx,
        innerwallr,
        outerwallx,
        outerwallr,
        -0.25 .+ hubcoordinates[:, 1],
        hubcoordinates[:, 2],
    )
    grid = DuctTAPE.initialize_grid(duct, rotors, grid_options)

    @test true == true
end
