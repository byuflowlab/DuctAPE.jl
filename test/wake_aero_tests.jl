@testset "Wake Aerodynamic Tests" begin
    @testset "Enthalpy Jump Calculations" begin

        # - Single Rotor - #
        Gammas = ones(2)
        blade_elements = [(omega=[1.0], num_blades=[2])]
        @test ones(2) / pi == dt.calculate_enthalpy_jumps(Gammas, blade_elements)

        # - Multi Rotor - #
        Gammas = ones(2, 2)
        blade_elements = [(omega=[1.0], num_blades=[2]) for i in 1:2]
        @test cumsum(ones(2, 2) / pi; dims=2) ==
            dt.calculate_enthalpy_jumps(Gammas, blade_elements)
    end

    @testset "Net Circulation Calculations" begin
        # - Single Rotor - #
        Gammas = ones(2)
        blade_elements = [(omega=[1.0], num_blades=[2])]
        @test ([2.0, 2.0], [2.0, 2.0]) ==
            dt.calculate_net_circulation(Gammas, blade_elements)

        # - Multi Rotor - #
        Gammas = ones(2, 2)
        blade_elements = [(omega=[1.0], num_blades=[2]) for i in 1:2]
        @test (2.0 * ones(2, 2), [2.0 4.0; 2.0 4.0]) ==
            dt.calculate_net_circulation(Gammas, blade_elements)
    end

    # @testset "Duct Surface Velocity Test" begin
    #     body_vortex_strengths = -ones(5)
    #     duct_panels = (panel_center=[2.0; 1.0; 0.0; 1.0; 2.0],)
    #     wake_panels = [(panel_center=[0.5; 1.5],)]
    #     @test ones(2) ==
    #         dt.get_surface_velocity(body_vortex_strengths, duct_panels, wake_panels)
    # end

    #TODO: add in rotor stuff when function is ready for it
    @testset "Wake Velocity Calculations" begin
        A_wake_to_wake = [[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] for i in 1:2, j in 1:2]
        gamma_wake = [1.0 2.0 3.0; 1.0 2.0 3.0]
        A_body_to_wake = [[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] for i in 1:2]
        gamma_body = [1.0; 2.0; 3.0]

        vm = dt.calculate_wake_velocities(
            A_body_to_wake, #vector of matrices size 1 x num wakes
            gamma_body,
            # A_rotor_to_wake, matrix of matrices size num rotors x num wakes
            # rotor_source_strengths,
            A_wake_to_wake, #matrix of matrices size num wakes x num wakes
            gamma_wake,
        )

        @test vm[1, :] ==
            A_wake_to_wake[1, 1] * gamma_wake[1, :] +
              A_wake_to_wake[2, 1] * gamma_wake[2, :] +
              A_body_to_wake[1] * gamma_body
        @test vm[2, :] ==
            A_wake_to_wake[1, 2] * gamma_wake[1, :] +
              A_wake_to_wake[2, 2] * gamma_wake[2, :] +
              A_body_to_wake[2] * gamma_body
    end

    @testset "Wake Vorticity Calculations" begin
        Vm = [1.0 1.0; 2.0 2.0; 3.0 3.0]
        Gamma_tilde = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0]
        H_tilde = [1.0 2.0; 2.0 3.0; 3.0 4.0; 4.0 5.0]
        rotoridxs = [1; 2]
        Vinf = 1.0
        wake_panels = [
            (panel_center=[0.0 1.0; 1.0 1.0],)
            (panel_center=[0.0 2.0; 1.0 2.0],)
            (panel_center=[0.0 3.0; 1.0 3.0],)
        ]

        wake_vortex_strengths = zeros(size(Vm))

        dt.calculate_wake_vortex_strengths!(
            wake_vortex_strengths, wake_panels, Vm, Vinf, Gamma_tilde, H_tilde, rotoridxs
        )

        @test wake_vortex_strengths[1, 1] == 1.0 / 2.0 * (-1.0 / (2.0 * pi) * 3.0 + 2.0)
        @test wake_vortex_strengths[2, 1] == 1.0 / 4.0 * (-1.0 / (4.0 * pi) * 5.0 + 2.0)
        @test wake_vortex_strengths[3, 1] == 1.0 / 6.0 * (-1.0 / (6.0 * pi) * 7.0 + 2.0)
        @test wake_vortex_strengths[1, 2] == 1.0 / 2.0 * (-1.0 / (2.0 * pi) * 5.0 + 2.0)
        @test wake_vortex_strengths[2, 2] == 1.0 / 4.0 * (-1.0 / (4.0 * pi) * 7.0 + 2.0)
        @test wake_vortex_strengths[3, 2] == 1.0 / 6.0 * (-1.0 / (6.0 * pi) * 9.0 + 2.0)
    end

    @testset "Wake Aero Initialization" begin

        # - Vm from Thrust - #
        freestream = (Vinf=1.0, rho=1.0)

        blade_elements = [
            (
                num_blades=1.0,
                omega=1.0,
                num_radial_stations=3,
                radial_positions=[1.0; 2.0; 3.0],
            ) for i in 1:2
        ]

        BGamma = ones(3, 2)
        Gamma_tilde = ones(3, 2)

        rotoridxs = [1, 2]

        wake_panels = [
            (panel_center=[1.0 1.5; 2.0 1.5; 3.0 1.5],)
            (panel_center=[1.0 2.5; 2.0 2.5; 3.0 2.5],)
        ]

        Vm = dt.vm_from_thrust(freestream, blade_elements, BGamma, wake_panels, rotoridxs)

        vin(vi) = sqrt(vi^2 + 1.0 / pi) - 0.5 * vi

        @test all(Vm[:, 1] .== vin(1.0) + 1.0)
        @test all(Vm[:, 2] .== vin(vin(1.0)) + 1.0)
        @test all(Vm[:, 3] .== Vm[:, 2] .* 0.99)

        # - initial vortex strengths - #
        params = (
            rotoridxs=rotoridxs,
            blade_elements=blade_elements,
            wake_panels=wake_panels,
            freestream=freestream,
        )
        gw = dt.initialize_wake_vortex_strengths(Gamma_tilde, params)
        @test gw == zeros(2, 3)
    end
end
