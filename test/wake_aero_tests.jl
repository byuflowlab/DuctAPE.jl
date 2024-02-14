@testset "Wake Aerodynamic Tests" begin
    @testset "Enthalpy Jump Calculations" begin

        # - Single Rotor - #
        Gammas = ones(2, 1)
        Omega = 1.0
        num_blades = 1
        @test all(
            ones(2) / (2 * pi) .== dt.calculate_enthalpy_jumps(Gammas, Omega, num_blades)
        )

        # - Multi Rotor - #
        Gammas = ones(2, 2)
        Omega = ones(2)
        num_blades = ones(2)
        @test cumsum(0.5 * ones(2, 2) / pi; dims=2) ==
            dt.calculate_enthalpy_jumps(Gammas, Omega, num_blades)
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
end
