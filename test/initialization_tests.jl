include("../data/dfdc/dstestr2_case.jl");

@testset "Grid Geometry" begin

    # - Split Wall Coordinates
    outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(ductx, ductr)

    # -- DEFINE DUCT OBJECT
    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    grid_options = DuctTAPE.defineGridOptions(10)

    grid = DuctTAPE.generate_grid_points(ductgeometry, ductsplines, [], grid_options)

    # Test to make sure there's no compilation errors for now.
    @test grid == grid

    #TODO; probably want to put together more basic unit tests to check ins and outs and such.
    #TODO: instead of going directly for full dfdc example, just define a simple square volume with flat walls.  You could probably just make it 1x1 or 2x2 or something small so you can put in exactly what the values for everything should be.

end

@testset "Grid Relaxation Stability" begin

    # - Split Wall Coordinates
    outerwallx, innerwallx, outerwallr, innerwallr = DuctTAPE.split_wall(ductx, ductr)

    # -- DEFINE DUCT OBJECT
    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    # -- GENERATE ROTOR OBJECT ARRAY

    #generate rotor object
    rotor1 = DuctTAPE.RotorGeometry(
        xdisk1,
        nblade1,
        rnondim1,
        0.0,
        chord1,
        beta1,
        nothing,
        nothing,
        nothing,
        nothing,
        rpm,
    )

    #generate stator object (rpm is zero for stator)
    rotor2 = DuctTAPE.RotorGeometry(
        xdisk2,
        nblade2,
        rnondim2,
        0.0,
        chord2,
        beta2,
        nothing,
        nothing,
        nothing,
        nothing,
        0.0,
    )

    #assemble array
    rotors = [rotor1; rotor2]

    #TODO: add test about rotor order, making sure xlocations are in increasing order.
    #TODO: will need to add this method somewhere (probably where the check is for this condition)

    #set up grid options
    num_radial_stations = length(rnondim1)
    grid_options = DuctTAPE.defineGridOptions(num_radial_stations)

    grid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

    # Test to make sure there's no compilation errors for now.
    @test grid == grid

    #TODO; probably want to put together more basic unit tests to check ins and outs and such.
end
