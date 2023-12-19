@testset "Initialization Checks" begin
    #---------------------------------#
    #      Duct and Hub Geometry      #
    #---------------------------------#
    # duct geometry
    duct_x = [1.0; 0.5; 0.0; 0.5; 1.0]
    duct_r = [0.0; -0.25; 0.0; 0.25; 0.0]
    duct_coordinates = [duct_x duct_r]

    # hub geometry
    hub_x = [0.0; 0.5; 1.0]
    hub_r = [0.0; 0.25; 0.0]
    hub_coordinates = [hub_x hub_r]

    #---------------------------------#
    #             ROTORS              #
    #---------------------------------#
    # Required rotor information for duct and wake generation
    rotorzlocs = [0.5; 0.75]

    # rotor parameters
    rotor1_parameters = (;
        Rtip=1.5,
        B=2,
        Omega=50,
        rotorzloc=rotorzlocs[1],
        r=[0.0, 1.0],
        chords=[0.5, 0.25],
        twists=[50.0, 10.0],
        airfoils=[nothing, nothing],
        tip_gap=0.0,
    )

    # stator parameters
    rotor2_parameters = (; rotor1_parameters..., rotorzloc=rotorzlocs[2])

    # array with rotor and stator parameters
    rotor_parameters = [rotor1_parameters, rotor2_parameters]

    nwake_sheets = 10
    nhub_inlet = 4
    nduct_inlet = 5
    npanels = [5; 5; 20] #this is a vector of number of panels between discrete points after the first rotor, e.g. between the rotor and the duct trailing edge and between the duct trailing edge and the end of the wake
    wake_length = 1.0

    paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

    constants = DuctAPE.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters, #vector of named tuples
        (Vinf=5.0,);#freestream;
        finterp=fm.linear,
    )

    # Test that the rotor tips and hubs ended up in the right places.
    for irotor in 1:(constants.num_rotors)
        @test constants.blade_elements[irotor].Rhub ==
            constants.rotor_panel_edges[1, irotor]
        @test constants.blade_elements[irotor].Rtip ==
            constants.rotor_panel_edges[end, irotor]
    end

    @test constants.blade_elements[1].Rhub == 0.25
    @test constants.blade_elements[1].Rtip == 1.5
    @test constants.blade_elements[2].Rhub == 0.125
    @test constants.blade_elements[2].Rtip == 1.625
end
