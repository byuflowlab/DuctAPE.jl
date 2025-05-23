println("\nAIRFOIL CORRECTION TESTS")
@testset "Stall Cutoffs" begin
    # - test that lift is monotonically increasing - #
    # load nominal data
    include("data/naca_4412_raw.jl")

    # apply stall cutoffs
    clext, cdext = dt.c4b.stall_limiters(alpha, cl, cdrag; cl_cutoff_slope=0.1, N=20, blend_hardness=50)

    @test ismonotonic(clext, 1)
end

@testset "Solidity and Stagger Corrections" begin

    # load wallis data
    include("data/wallis_fig6-29_data.jl")

    # check 1/σ=0.5
    solidity = 1.0 / 0.5
    stagger = oneoversig05[:, 1] * pi / 180.0
    for (i, s) in enumerate(stagger)
        @test isapprox(
            oneoversig05[i, 2], dt.c4b.solidity_and_stagger_factor_smooth(solidity, s), atol=1e-2
        )
    end

    # check 1/σ=1.0
    solidity = 1.0
    stagger = oneoversig10[:, 1] * pi / 180.0
    for (i, s) in enumerate(stagger)
        @test isapprox(
            oneoversig10[i, 2], dt.c4b.solidity_and_stagger_factor_smooth(solidity, s), atol=1e-2
        )
    end

    # check 1/σ=1.5
    solidity = 1.0 / 1.5
    stagger = oneoversig15[:, 1] * pi / 180.0
    for (i, s) in enumerate(stagger)
        @test isapprox(
            oneoversig15[i, 2], dt.c4b.solidity_and_stagger_factor_smooth(solidity, s), atol=1e-2
        )
    end

    # check that factor is constant below 20 degrees
    solidity = 2.0
    @test isapprox(
        dt.c4b.solidity_and_stagger_factor_smooth(solidity, 10.0 * pi / 180),
        dt.c4b.solidity_and_stagger_factor_smooth(solidity, 11.0 * pi / 180),
        atol=1e-8,
    )

    # check that factor maxes out at 1.0
    @test 1.0 == dt.c4b.solidity_and_stagger_factor_smooth(0.5, 80.0 * pi / 180)
end

@testset "Prandtl-Glauert Compressibility Corrections" begin
    # - Check that things aren't weird - #
    @test dt.c4b.prandtl_glauert_factor(0.1) == 1.0 / sqrt(1.0 - 0.1^2)

    # - Check that M=1 doesn't break - #
    @test !isinf(dt.c4b.prandtl_glauert_factor(1.0))

    # - Check that M>1 doesn't break - #
    @test isreal(dt.c4b.prandtl_glauert_factor(2.0))

    # - Check that M>1 behavior is as expected - #
    @test dt.c4b.prandtl_glauert_factor(1.1) == 1.0 / sqrt(1.0 - 0.99^2)
end

@testset "Reynolds Number Corrections" begin
    # - check that nothing is weird - #
    @test dt.c4b.re_drag(1.0, 2e6, 5e6; re_exp=0.5) == (5e6 / 2e6)^0.5
end

@testset "Transonic Effects Corrections" begin

    # - Check that lift is being limited after critical mach - #
    @test all(
        dt.c4b.transonic_lift_limiter_smooth!([-1.0, 1.0], 0.9, 0.1, 1.0, -1.0, 2 * pi) .< [1.0]
    )
    @test all(dt.c4b.transonic_lift_limiter_smooth!([-1.0], 0.9, 0.1, 1.0, -1.0, 2 * pi) .> [-1.0])
end
