include("../data/dfdc/dstestr2_case.jl");

@testset "Example Duct Geometry" begin

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
end

@testset "Example Grid Geometry" begin

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

    grid = DuctTAPE.initialize_wakegrid(ductgeometry, ductsplines, rotors, grid_options)

    # Test to make sure there's no compilation errors for now.
    @test grid == grid
end

@testset "Example Airfoils" begin

    #test that airfoil functions work
    datapath = "../data/dfdc/airfoils/"

    # generate airfoils for rotor
    af1 = DuctTAPE.ccb.AlphaReAF([
        datapath * "disk1_re5e5.dat",
        datapath * "disk1_re1e6.dat",
        datapath * "disk1_re1.5e6.dat",
        datapath * "disk1_re2e6.dat",
    ])

    cl, cd = DuctTAPE.get_clcd(af1, 5.0, 1.5e6, 0.0, -1.0)
    cldata = 0.5465618348626059
    cddata = 0.02635818567589299

    @test isapprox(cl, cldata)
    @test isapprox(cd, cddata)
end

@testset "Example Rotors and Blades" begin

    #test that airfoil functions work
    datapath = "../data/dfdc/airfoils/"

    # generate airfoils for rotor
    af1 = DuctTAPE.ccb.AlphaReAF([
        datapath * "disk1_re5e5.dat",
        datapath * "disk1_re1e6.dat",
        datapath * "disk1_re1.5e6.dat",
        datapath * "disk1_re2e6.dat",
    ])

    rs = [0.0; 0.25; 0.5; 0.75; 1.0]
    chords = [0.1; 0.1; 0.1; 0.1; 0.1]
    twists = [90.0; 80.0; 75.0; 72.5; 71.25]                 #test rotor and blade definitions.
    rotor = DuctTAPE.RotorGeometry(
        0.5, 5, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 5), nothing, 2000
    )

    # nothing really to test, this is a direct definition...
    @test rotor == rotor

    outerwallx = [0; 0.5; 1]
    innerwallx = outerwallx
    outerwallr = [1; 1.5; 1]
    innerwallr = [1.0; 1; 1]
    hubx = outerwallx
    hubr = [0; 0; 0.0]

    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    bladedims = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor)

    @test bladedims.hubr == 0.0
    @test bladedims.tipr == 1.0
    @test bladedims.rdim == rs
    @test bladedims.tdim == twists
    @test isapprox(bladedims.sweptannulus, pi)
    @test isapprox(bladedims.sweptarea, pi)
end

@testset "Example Paneling Tests:" begin

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

    wakegrid = DuctTAPE.initialize_wakegrid(ductgeometry, ductsplines, rotors, grid_options)

    # get paneling of various objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = DuctTAPE.generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    # Test to make sure there's no compilation errors for now.
    @test true

    #TODO: add unit tests for the various pieces of the paneling procedure.  May need to break it up into several functions.
end
