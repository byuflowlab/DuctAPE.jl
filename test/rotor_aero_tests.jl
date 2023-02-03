@testset "Rotor Aerodynamics" begin
    @testset "Angle of Attack Calculation" begin
        twist_angle = pi / 2.0
        meridional_flow_magnitude = 1.0
        tangential_flow_magnitude = -1.0

        @test pi / 4.0 == dt.calculate_angle_of_attack(
            twist_angle, meridional_flow_magnitude, tangential_flow_magnitude
        )
    end

    @testset "Blade Element Inflow Velocity Calculation" begin

        # - Single Blade Element - #
        Vinf = 1.0
        vm = [1.0]
        vtheta = [-1.0]
        blade_elements = [(radial_positions=[1.0], omega=[1.0])]
        inflow = dt.calculate_inflow_velocities(blade_elements, Vinf, vm, vtheta)
        @test inflow.Wm[1] == 2.0
        @test inflow.Wtheta[1] == -2.0
        @test inflow.Wmag[1] == sqrt(8.0)

        # - Single Rotor - #
        Vinf = 1.0
        vm = ones(2)
        vtheta = -ones(2)
        blade_elements = [(radial_positions=ones(2), omega=ones(2))]
        inflow = dt.calculate_inflow_velocities(blade_elements, Vinf, vm, vtheta)
        @test all(inflow.Wm[:] .== 2.0)
        @test all(inflow.Wtheta[:] .== -2.0)
        @test all(inflow.Wmag[:] .== sqrt(8.0))

        # - Multiple Rotors - #
        Vinf = 1.0
        vm = ones(2, 2)
        vtheta = -ones(2, 2)
        blade_elements = [(radial_positions=ones(2), omega=ones(2)) for i in 1:2]
        inflow = dt.calculate_inflow_velocities(blade_elements, Vinf, vm, vtheta)
        @test all(inflow.Wm[:] .== 2.0)
        @test all(inflow.Wtheta[:] .== -2.0)
        @test all(inflow.Wmag[:] .== sqrt(8.0))
    end

    @testset "Airfoil Polar Search" begin
        alpha = [-pi; 0.0; pi]
        cl = 2.0 * ones(length(alpha))
        cd = ones(length(alpha))
        af = ccb.AlphaAF(alpha, cl, cd, "TEST")
        @test dt.search_polars(af, pi / 2.0) == (2.0, 1.0)
    end

    @testset "Blade Element Circulation and Source Strength Calculation" begin
        alpha = [-pi; 0.0; pi]
        cl = 2.0 * ones(length(alpha))
        cd = ones(length(alpha))
        af = AlphaAF(alpha, cl, cd, "TEST")

        Vinf = 1.0
        vm = ones(2, 2)
        vtheta = -ones(2, 2)
        blade_elements = [
            (
                num_radial_stations=[2],
                chords=ones(2),
                twists=pi / 2.0 * ones(2),
                radial_positions=ones(2),
                omega=ones(2),
                airfoils=fill(af, 2),
            ) for i in 1:2
        ]

        Gamma = similar(vm)
        Sigma = similar(vm)

        dt.calculate_gamma_sigma!(Gamma, Sigma, blade_elements, Vinf, vm, vtheta)

        @test all(Gamma .== sqrt(8))
        @test all(Sigma .== 0.5 * sqrt(8))

        Gamma, Sigma = dt.calculate_gamma_sigma(blade_elements, Vinf, vm, vtheta)
        @test all(Gamma .== sqrt(8))
        @test all(Sigma .== 0.5 * sqrt(8))
    end

    @testset "Induced Velocity at Blade Elements Calculation" begin

        # - Single Rotor - #

        BGamma = ones(2)
        Gamma_tilde = ones(2)
        blade_elements = [(
            num_radial_stations=[2],
            chords=ones(2),
            twists=pi / 2.0 * ones(2),
            radial_positions=ones(2),
            omega=ones(2),
        )]
        A_bodies_to_rotor = [1 0; 0 1]
        gamma_bodies = ones(2)
        A_wake_to_rotor = [[1 0; 0 1] for i in 1:2]
        gamma_wake = ones(2, 2)
        # A_rotor_to_rotor =[1 0; 0 1]
        # Sigma =ones(2)

        vi = dt.calculate_induced_velocities(
            BGamma,
            Gamma_tilde,
            blade_elements,
            A_bodies_to_rotor,
            gamma_bodies,
            A_wake_to_rotor,
            gamma_wake,
            # A_rotor_to_rotor,
            # Sigma,
        )

        @test vi.vtheta == 1.5 / (2.0 * pi) * ones(2)
        @test vi.vm == 3.0 * ones(2)

        # - Multiple Rotors - #
        BGamma = [1.0 for i in 1:2, j in 1:2]
        Gamma_tilde = [1.0 for i in 1:2, j in 1:2]
        blade_elements = [
            (
                num_radial_stations=[2],
                chords=ones(2),
                twists=pi / 2.0 * ones(2),
                radial_positions=ones(2),
                omega=ones(2),
            ) for i in 1:2
        ]
        A_bodies_to_rotor = [[1 0; 0 1] for i in 1:2]
        gamma_bodies = ones(2)
        A_wake_to_rotor = [[1 0; 0 1] for i in 1:2, j in 1:2]
        gamma_wake = ones(2, 2)
        # A_rotor_to_rotor =[1 0; 0 1]
        # Sigma =ones(2)

        vi = dt.calculate_induced_velocities(
            BGamma,
            Gamma_tilde,
            blade_elements,
            A_bodies_to_rotor,
            gamma_bodies,
            A_wake_to_rotor,
            gamma_wake,
            # A_rotor_to_rotor,
            # Sigma,
        )

        @test vi.vtheta == 1.5 / (2.0 * pi) * ones(2, 2)
        @test vi.vm == 3.0 * ones(2, 2)
    end
end

