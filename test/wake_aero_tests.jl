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

    @testset "Duct Surface Velocity Test" begin
        body_vortex_strengths = -ones(5)
        duct_panels = (panel_center=[2.0; 1.0; 0.0; 1.0; 2.0],)
        wake_panels = [(panel_center=[0.5; 1.5],)]
        @test ones(2) ==
            dt.get_surface_velocity(body_vortex_strengths, duct_panels, wake_panels)
    end

    @testset "Wake Velocity Calculations" begin

        # - Single Rotor - #
        x_edge = [1.0; 2.0; 3.0; 4.0; 5.0]
        edge_velocity = ones(5)
        rotoridxs = [1]
        Gamma_tilde = [1.0; 2.0]
        H_tilde = [1.0; 2.0]
        blade_elements = [(radial_positions=ones(2),)]

        vm = dt.calculate_wake_velocities(
            x_edge, edge_velocity, rotoridxs, Gamma_tilde, H_tilde, blade_elements
        )

        @test vm ==
            [sqrt(1.0 + (1.0 / (2.0 * pi))^2 * 3.0 + 2.0) * ones(5)'; edge_velocity']

        # - Multi Rotor - #
        x_edge = [1.0; 2.0; 3.0; 4.0; 5.0]
        edge_velocity = ones(5)
        rotoridxs = [1, 3]
        Gamma_tilde = [1.0 1.0; 2.0 2.0]
        H_tilde = [1.0 1.0; 2.0 2.0]
        blade_elements = [(radial_positions=ones(2),) for i in 1:2]

        vm = dt.calculate_wake_velocities(
            x_edge, edge_velocity, rotoridxs, Gamma_tilde, H_tilde, blade_elements
        )

        @test vm ==
            [sqrt(1.0 + (1.0 / (2.0 * pi))^2 * 3.0 + 2.0) * ones(5)'; edge_velocity']
    end

    @testset "Wake Vorticity Calculations" begin
        vm = ones(2, 2)
        @test [1.0 1.0; 0.0 0.0] == dt.calculate_wake_vorticity(vm)
    end
end
