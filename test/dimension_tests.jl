#TODO; will likely want to update these to make sure all the indexing is being used correctly throughout.

@testset "I/O Dimension Tests" begin


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
