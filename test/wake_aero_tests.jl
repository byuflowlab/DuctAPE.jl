@testset "Wake Aerodynamic Tests" begin
    @testset "Enthalpy Jump Calculations" begin

        # - Single Rotor - #
        Gammas = ones(2, 1)
        omega = 1.0
        num_blades = 1
        @test all(
            ones(2) / (2 * pi) .== dt.calculate_enthalpy_jumps(Gammas, omega, num_blades)
        )

        # - Multi Rotor - #
        Gammas = ones(2, 2)
        omega = ones(2)
        num_blades = ones(2)
        @test cumsum(0.5 * ones(2, 2) / pi; dims=2) ==
            dt.calculate_enthalpy_jumps(Gammas, omega, num_blades)
    end

    @testset "Net Circulation Calculations" begin
        # - Single Rotor - #
        Gammas = ones(2, 1)
        num_blades = 2
        @test all([2.0, 2.0] .== dt.calculate_net_circulation(Gammas, num_blades))

        # - Multi Rotor - #
        Gammas = ones(2, 2)
        num_blades = 2.0 * ones(2)
        @test [2.0 4.0; 2.0 4.0] == dt.calculate_net_circulation(Gammas, num_blades)
    end

    @testset "Wake Vorticity Calculations" begin
        # - Single Rotor - #
        Vm = ones(3, 1) .* [1.0; 2.0; 3.0]
        Gamma_tilde = [1.0; 2.0; 3.0]
        H_tilde = [1.0; 2.0; 3.0]
        Rr_wake = ones(4, 1)

        wake_vortex_strengths = similar(Rr_wake) .= 0
        dt.calculate_wake_vortex_strengths!(
            wake_vortex_strengths, Rr_wake, Vm, Gamma_tilde, H_tilde
        )

        @test isapprox(
            wake_vortex_strengths[1], 1.0 / 2.0 * (-(1.0 / (2.0 * pi))^2 * 1.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[2], 1.0 / 3.0 * (-(1.0 / (2.0 * pi))^2 * 3.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[3], 1.0 / 5.0 * (-(1.0 / (2.0 * pi))^2 * 5.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[4], 1.0 / 6.0 * ((1.0 / (2.0 * pi))^2 * 9.0 - 6.0)
        )

        # - Multi Rotor - #
        Vm = ones(3, 2) .* [1.0; 2.0; 3.0]
        Gamma_tilde = ones(3, 2) .* [1.0; 2.0; 3.0]
        H_tilde = ones(3, 2) .* [1.0; 2.0; 3.0]
        Rr_wake = ones(4, 2)

        wake_vortex_strengths = similar(Rr_wake) .= 0
        dt.calculate_wake_vortex_strengths!(
            wake_vortex_strengths, Rr_wake, Vm, Gamma_tilde, H_tilde
        )

        @test isapprox(
            wake_vortex_strengths[1, 2], 1.0 / 2.0 * (-(1.0 / (2.0 * pi))^2 * 1.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[2, 2], 1.0 / 3.0 * (-(1.0 / (2.0 * pi))^2 * 3.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[3, 2], 1.0 / 5.0 * (-(1.0 / (2.0 * pi))^2 * 5.0 + 2.0)
        )
        @test isapprox(
            wake_vortex_strengths[4, 2], 1.0 / 6.0 * ((1.0 / (2.0 * pi))^2 * 9.0 - 6.0)
        )
    end
end
