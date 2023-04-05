@testset "I/O Dimension Tests" begin

    # - Start with Initialization Function - #
    duct_coordinates = [1.0 1.0; 0.5 0.75; 0.0 1.0; 0.5 1.25; 1.0 1.0]
    hub_coordinates = [0.0 0.0; 0.5 0.25; 1.0 0.0]
    freestream = dt.Freestream(10.0)
    af = ccb.AlphaAF("./data/naca_4412_extrapolated_rotated_APCshifted.dat")
    airfoils = fill(af, 3)
    rotor_parameters = [
        (
            rotor_x_position=0.5,
            radial_positions=[0.0; 0.5; 1.0],
            chords=ones(3),
            twists=ones(3),
            airfoils=airfoils,
            num_radial_stations=3,
            num_blades=2,
            Omega=100.0,
        )
        (
            rotor_x_position=0.75,
            radial_positions=[0.0; 0.5; 1.0],
            chords=ones(3),
            twists=ones(3),
            airfoils=airfoils,
            num_radial_stations=3,
            num_blades=2,
            Omega=0.0,
        )
    ]

    params = dt.initialize_parameters(
        duct_coordinates, hub_coordinates, rotor_parameters, freestream
    )

    @test size(params.converged) == (1,)

    @test size(params.bc_freestream_to_body) == (7,)
    @test size(params.blade_elements) == (2,)
    @test size(params.rotoridxs) == (2,)
    @test size(params.body_panels) == 2
    @test size(params.rotor_panels) == (2,)
    @test size(params.wake_panels) == (2,)
    @test size(params.A_body_to_body) == (7, 7)
    @test size(params.A_body_to_rotor) == (2,)
    @test size(params.A_body_to_wake) == (2,)
    @test size(params.A_rotor_to_body) == (2,)
    @test size(params.A_rotor_to_rotor) == (2, 2)
    @test size(params.A_rotor_to_wake) == (2, 2)
    @test size(params.A_wake_to_body) == (2,)
    @test size(params.A_wake_to_rotor) == (2, 2)
    @test size(params.A_wake_to_wake) == (2, 2)

    # - No Rotor Solution - #

    body_vortex_strengths = ImplicitAD.implicit_linear(
        params.A_body_to_body, params.bc_freestream_to_body
    )

    @test size(body_vortex_strengths) == (7,)

    # - Rotor Circulation and Panel Source Strengths - #
    rotor_circulation_strengths, rotor_panel_source_strengths = dt.calculate_gamma_sigma(
        params.blade_elements, params.freestream.Vinf
    )

    @test size(rotor_circulation_strengths) == (3, 2)
    @test size(rotor_panel_source_strengths) == (3, 2)

    # - Wake Vortex Strengths - #
    wake_vortex_strengths = dt.initialize_wake_vortex_strengths(
        rotor_circulation_strengths, params
    )

    @test size(wake_vortex_strengths) == (2, 5)

    # - Variable Extraction - #

    state_variables = [
        body_vortex_strengths[1:(end - 1)] # Don't include the bound circulation value used in the kutta condition.
        reduce(vcat, wake_vortex_strengths') # wake_vortex_strengths comes out as a matrix, one ROW for each wake, need to make it an array.
        reduce(vcat, rotor_circulation_strengths) # Gamma_init will be defined as a matrix, one COLUMN for each rotor, need to reduce to a single vector
        # reduce(vcat,blade_element_source_strengths) # Sigma_init will be defined as a matrix, one COLUMN for each rotor, need to reduce to a single vector
    ]

    body_vortex_strengths, wake_vortex_strengths, rotor_circulation_strengths, blade_element_source_strengths = dt.extract_state_variables(
        state_variables, params
    )

    @test size(body_vortex_strengths) == (6,)
    @test size(rotor_circulation_strengths) == (3, 2)
    @test size(blade_element_source_strengths) == (3, 2)
    @test size(wake_vortex_strengths) == (2, 5)

    rotor_panel_source_strengths =
        (
            blade_element_source_strengths[1:(end - 1), :] .+
            blade_element_source_strengths[2:end, :]
        ) / 2.0

    @test size(rotor_panel_source_strengths) == (2, 2)

    # - Jumps - #

    H_tilde = dt.calculate_enthalpy_jumps(
        rotor_circulation_strengths, params.blade_elements
    )

    BGamma, Gamma_tilde = dt.calculate_net_circulation(
        rotor_circulation_strengths, params.blade_elements
    )

    @test size(H_tilde) == (3, 2)
    @test size(BGamma) == (3, 2)
    @test size(Gamma_tilde) == (3, 2)

    # - Rotor Velocities - #

    vi_rotor = dt.calculate_induced_velocities_on_rotors(
        BGamma,
        Gamma_tilde,
        params.blade_elements,
        params.A_body_to_rotor,
        body_vortex_strengths,
        params.A_wake_to_rotor,
        wake_vortex_strengths,
        # params.A_rotor_to_rotor,
        # rotor_panel_source_strengths, # these are the averages
    )

    @test size(vi_rotor.vm) == (3, 2)
    @test size(vi_rotor.vtheta) == (3, 2)

    # - wake velocities - #

    wake_velocities = dt.calculate_wake_velocities(
        params.A_body_to_wake,
        body_vortex_strengths,
        params.A_wake_to_wake,
        wake_vortex_strengths,
        # params.A_rotor_to_wake,
        # rotor_panel_source_strengths, # these are the averages
    )

    @test size(wake_velocities) == (2, 5)
end
