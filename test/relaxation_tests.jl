#=
Tests for controlled successive under relaxation functions
=#

@testset "Circulation Relaxation" begin
    Gamr = ones(3, 2)
    deltaG_prev = zeros(size(Gamr))
    deltaG = ones(size(Gamr))
    maxBGamr = zeros(2)
    maxdeltaBGamr = zeros(2)
    B = [3, 2]

    dt.relax_Gamr!(
        Gamr,
        deltaG_prev,
        deltaG,
        maxBGamr,
        maxdeltaBGamr,
        B;
        nrf=0.5,
        bt1=0.2,
        bt2=0.6,
        pf1=0.4,
        pf2=0.5,
    )

    @test isapprox(Gamr, [1.72 1.48; 1.72 1.48; 1.72 1.48])
    @test maxBGamr == [3.0, 2.0]
    @test isapprox(maxdeltaBGamr, [0.72, 0.48])
end

@testset "Wake Strength Relaxation" begin
    gamw = ones(3, 2)
    deltag_prev = zeros(size(Gamr))
    deltag = ones(size(Gamr))
    maxdeltagamw = Array{Float64,0}(undef)
    maxdeltagamw[] = 0.0

    dt.relax_gamw!(gamw, deltag_prev, deltag, maxdeltagamw; nrf=0.5, btw=0.6, pfw=1.2)

    @test gamw == 1.6 * ones(3, 2)
    @test maxdeltagamw != 1.0
end

@testset "Convergence Criteria" begin
    conv = MVector{1,Bool}(false)
    maxBGamr = ones(2)
    maxdeltaBGamr = ones(2)
    maxdeltagamw = Array{Float64,0}(undef)
    maxdeltagamw[] = 0.0
    Vref = 1.0
    dt.check_convergence!(
        conv, maxBGamr, maxdeltaBGamr, maxdeltagamw, Vref; f_circ=1e-3, f_dgamw=2e-4
    )
    @test conv[] == false

    dt.check_convergence!(
        conv, 1e4 * ones(2), maxdeltaBGamr, maxdeltagamw, Vref; f_circ=1e-3, f_dgamw=2e-4
    )
    @test conv[] == true
end
