@testset "I/O Dimension Tests" begin

    # - Start with Initialization Function - #
    duct_coordinates = [1.0 1.0; 0.5 0.75; 0.0 1.0; 0.5 1.25; 1.0 1.0]
    hub_coordinates = [0.0 0.0; 0.5 0.25; 1.0 0.0]
    freestream = dt.Freestream(10.0)
    af = ccb.AlphaAF("data/naca4412.dat")
    airfoils = fill(af, 3)
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
    )

    # stator parameters
    rotor2_parameters = (; rotor1_parameters..., rotorzloc=rotorzlocs[2])

    # array with rotor and stator parameters
    rotor_parameters = [rotor1_parameters, rotor2_parameters]

    constants = DuctAPE.initialize_constants(
        duct_coordinates,
        hub_coordinates,
        rotor_parameters, #vector of named tuples
        (Vinf=5.0,);#freestream;
        wake_length=1.0,
        nwake_sheets=10,
        finterp=fm.linear,
        nhub_inlet=4,
        nduct_inlet=5,
        npanels=[5; 5; 20], #this is a vector of number of panels between discrete points after the first rotor, e.g. between the rotor and the duct trailing edge and between the duct trailing edge and the end of the wake
    )

    @test size(constants.converged) == (1,)

    @test size(constants.b_bf) == (44,)
    @test size(constants.blade_elements) == (2,)
    @test size(constants.body_panels) == (2,)
    @test size(constants.rotor_source_panels) == (2,)
    @test size(constants.wake_vortex_panels) == (10,)
    @test size(constants.A_bb) == (44, 44)
    @test size(constants.vx_rb) == (2, 1)
    @test size(constants.vr_rb) == (2, 1)
    @test size(constants.A_br) == (1, 2)
    @test size(constants.vx_rr) == (2, 2)
    @test size(constants.vr_rr) == (2, 2)
    @test size(constants.A_bw) == (1,10)
    @test size(constants.vx_rw) == (2, 10)
    @test size(constants.vr_rw) == (2, 10)


    # TODO: need to add body and source stuff
    # # - Rotor Circulation and Panel Source Strengths - #
    # rotor_circulation_strengths, rotor_panel_source_strengths = dt.calculate_gamma_sigma(
    #     constants.blade_elements, constants.freestream.Vinf
    # )

    # @test size(rotor_circulation_strengths) == (3, 2)
    # @test size(rotor_panel_source_strengths) == (3, 2)

    # - Variable Extraction - #
    # TODO: update

    # state_variables = [
    #     body_vortex_strengths[1:(end - 1)] # Don't include the bound circulation value used in the kutta condition.
    #     reduce(vcat, wake_vortex_strengths') # wake_vortex_strengths comes out as a matrix, one ROW for each wake, need to make it an array.
    #     reduce(vcat, rotor_circulation_strengths) # Gamma_init will be defined as a matrix, one COLUMN for each rotor, need to reduce to a single vector
    #     # reduce(vcat,blade_element_source_strengths) # Sigma_init will be defined as a matrix, one COLUMN for each rotor, need to reduce to a single vector
    # ]

    # body_vortex_strengths, wake_vortex_strengths, rotor_circulation_strengths, blade_element_source_strengths = dt.extract_state_variables(
    #     state_variables, constants
    # )

    # @test size(body_vortex_strengths) == (6,)
    # @test size(rotor_circulation_strengths) == (3, 2)
    # @test size(blade_element_source_strengths) == (3, 2)
    # @test size(wake_vortex_strengths) == (2, 5)

    # rotor_panel_source_strengths =
    #     (
    #         blade_element_source_strengths[1:(end - 1), :] .+
    #         blade_element_source_strengths[2:end, :]
    #     ) / 2.0

    # @test size(rotor_panel_source_strengths) == (2, 2)

    # - Jumps - #

    # H_tilde = dt.calculate_enthalpy_jumps(
    #     rotor_circulation_strengths, constants.blade_elements
    # )

    # BGamma, Gamma_tilde = dt.calculate_net_circulation(
    #     rotor_circulation_strengths, constants.blade_elements
    # )

    # @test size(H_tilde) == (3, 2)
    # @test size(BGamma) == (3, 2)
    # @test size(Gamma_tilde) == (3, 2)

    # # - Rotor Velocities - #

    # vi_rotor = dt.calculate_induced_velocities_on_rotors(
    #     BGamma,
    #     Gamma_tilde,
    #     constants.blade_elements,
    #     constants.A_body_to_rotor,
    #     body_vortex_strengths,
    #     constants.A_wake_to_rotor,
    #     wake_vortex_strengths,
    #     # constants.A_rotor_to_rotor,
    #     # rotor_panel_source_strengths, # these are the averages
    # )

    # @test size(vi_rotor.vm) == (3, 2)
    # @test size(vi_rotor.vtheta) == (3, 2)

    # # - wake velocities - #

    # wake_velocities = dt.calculate_wake_velocities(
    #     constants.A_body_to_wake,
    #     body_vortex_strengths,
    #     constants.A_wake_to_wake,
    #     wake_vortex_strengths,
    #     # constants.A_rotor_to_wake,
    #     # rotor_panel_source_strengths, # these are the averages
    # )

    # @test size(wake_velocities) == (2, 5)
end
