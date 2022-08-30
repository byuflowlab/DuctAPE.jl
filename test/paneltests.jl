@testset "Paneling Tests:" begin

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

    wakegrid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

    # get paneling of various objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = DuctTAPE.generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    # Test to make sure there's no compilation errors for now.
    @test true

    #TODO: add unit tests for the various pieces of the paneling procedure.  May need to break it up into several functions.
end
