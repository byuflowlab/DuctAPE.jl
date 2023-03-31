@testset "Rotor Aerodynamics" begin
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
        af = ccb.AlphaAF(alpha, cl, cd, "TEST")

        Vinf = 1.0
        vm = ones(2, 2)
        vtheta = -ones(2, 2)
        blade_elements = [
            (
                num_radial_stations=[2],
                num_blades=1,
                chords=ones(2),
                twists=pi / 2.0 * ones(2),
                radial_positions=ones(2),
                omega=ones(2),
                inner_airfoil=fill(af, 2),
                outer_airfoil=fill(af, 2),
                inner_fraction=[1.0; 1.0],
            ) for i in 1:2
        ]

        Gamma = similar(vm) .= 0
        Sigma = similar(vm) .= 0

        W = sqrt.(vm .^ 2 .+ vtheta .^ 2)

        dt.calculate_gamma_sigma!(Gamma, Sigma, blade_elements, vm, vtheta)

        @test all(Gamma .== W)
        @test all(Sigma .== 1.0 / (4.0 * pi) * W)

        Gamma, Sigma = dt.calculate_gamma_sigma(blade_elements, vm, vtheta)
        @test all(Gamma .== W)
        @test all(Sigma .== 1.0 / (4.0 * pi) * W)
    end

    @testset "Induced Velocity at Blade Elements Calculation" begin

        # - Single Rotor - #

        gamma_rotor = ones(2, 1)
        blade_elements = [(
            num_radial_stations=[2],
            num_blades=1,
            chords=ones(2),
            twists=pi / 2.0 * ones(2),
            radial_positions=ones(2),
            omega=ones(2),
        )]
        vx_rw = [[1 0; 0 1] for i in 1:1, j in 1:2]
        vr_rw = [[1 0; 0 1] for i in 1:1, j in 1:2]
        gamma_wake = ones(2, 2)
        # A_rotor_to_rotor =[1 0; 0 1]
        # Sigma =ones(2)

        vx, vr, vt = dt.calculate_induced_velocities_on_rotors(
            blade_elements, gamma_rotor, vx_rw, vr_rw, gamma_wake
        )

        @test all(vt .== 1.0 / (4.0 * pi) * gamma_rotor)
        @test all(vx .== 2.0 * ones(2))
        @test all(vr .== 2.0 * ones(2))

        # - Multiple Rotors - #
        gamma_rotor = ones(2, 2)
        blade_elements = [
            (
                num_radial_stations=[2],
                chords=ones(2),
                twists=pi / 2.0 * ones(2),
                radial_positions=ones(2),
                omega=ones(2),
                num_blades=1,
            ) for i in 1:2
        ]
        gamma_wake = ones(2, 2)
        vx_rw = [[1 0; 0 1] for i in 1:2, j in 1:2]
        vr_rw = [[1 0; 0 1] for i in 1:2, j in 1:2]

        vx, vr, vt = dt.calculate_induced_velocities_on_rotors(
            blade_elements, gamma_rotor, vx_rw, vr_rw, gamma_wake
        )

        @test all(vt[:, 1] .== 1.0 / (4.0 * pi) * gamma_rotor)
        @test all(vt[:, 2] .== (1.0 / (4.0 * pi) + 1 / (2 * pi)) * gamma_rotor)
        @test all(vx .== 2.0 * ones(2))
        @test all(vr .== 2.0 * ones(2))
    end
end

