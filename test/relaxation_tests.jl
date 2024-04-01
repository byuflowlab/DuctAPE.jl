#=
Tests for controlled successive under relaxation functions
=#
println("\nCSUR TESTS")

@testset "Circulation Relaxation" begin
    Gamr = ones(3, 2)
    deltaG_prev = zeros(size(Gamr))
    deltaG = ones(size(Gamr))
    maxBGamr = zeros(2)
    maxdeltaBGamr = zeros(2)
    B = [3, 2]

    Gamr, bladeomega, omega = dt.relax_Gamr!(
        Gamr,
        deltaG_prev,
        deltaG,
        maxBGamr,
        maxdeltaBGamr,
        B;
        nrf=0.4,
        bt1=0.2,
        bt2=0.6,
        pf1=0.4,
        pf2=0.5,
        test=true,
    )

    @test bladeomega == [0.4, 0.4] # note this is for each blade
    @test omega == [0.24, 0.24, 0.24] # note this is for the last blade
    @test isapprox(Gamr, 1.24 * ones(size(Gamr)))
    @test maxBGamr == [3.0, 2.0]
    @test isapprox(maxdeltaBGamr, [3.0, 2.0]) #ΔΓ is 1.0, and this is ΔBΓ so it's just the number of blades
end

@testset "Wake Strength Relaxation" begin
    gamw = ones(3)
    deltag_prev = zeros(size(gamw))
    deltag = ones(size(gamw))
    maxdeltagamw = Array{Float64,0}(undef)
    maxdeltagamw[] = 0.0

    dt.relax_gamw!(gamw, deltag_prev, deltag, maxdeltagamw; nrf=0.4, btw=0.6, pfw=1.2)

    @test gamw == 1.48 * ones(3)
    @test maxdeltagamw[] == 1.0
end

@testset "Convergence Criteria" begin
    conv = [false]
    maxBGamr = ones(2)
    maxdeltaBGamr = ones(2)
    maxdeltagamw = Array{Float64,0}(undef)
    maxdeltagamw[] = 0.0
    Vref = 1.0
    dt.check_CSOR_convergence!(
        conv,
        [maximum(abs.(maxdeltaBGamr ./ maxBGamr)); maxdeltagamw];
        f_circ=1e-3,
        f_dgamw=2e-4,
        convergence_type=dt.Relative(),
        verbose=false,
    )
    @test conv[] == false

    dt.check_CSOR_convergence!(
        conv,
        1e-4 * ones(2);
        f_circ=1e-3,
        f_dgamw=2e-4,
        convergence_type=dt.Relative(),
        verbose=false,
    )
    @test conv[] == true
end

@testset "DFDC Relaxation Comparison" begin
    include("data/dfdc_csor/bin/csor_input_values.jl")
    include("data/dfdc_csor/bin/csor_output_values.jl")

    # - Circulation Relaxation - #
    B = [2]
    Gamr_est = BGX' ./ B[1]
    Gamr = BGAM' ./ B[1]
    deltaG_prev = DBGOLD' .* 1.0
    deltaG = (Gamr_est .- Gamr)
    maxBGamr = [0.0]
    maxdeltaBGamr = [0.0]

    Gamr, bladeomega, omega = dt.relax_Gamr!(
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
        test=true,
    )

    @test isapprox(maxdeltaBGamr[1], DBGMAX)
    @test isapprox(maxBGamr[1], BGMAG)
    @test isapprox(bladeomega[1], RLXB)
    for (o, odfdc) in zip(omega, RLXBG)
        @test isapprox(o, odfdc, atol=1e-6)
    end
    for (g, gdfdc) in zip(Gamr, BGAMupdate' ./ B[1])
        @test isapprox(g, gdfdc, atol=1e-6)
    end

    # - Wake Relaxation - #
    gamw_est = GAMTH' .* 1.0
    gamw = GTH' .* 1.0
    deltag_prev = DGOLD' .* 1.0
    deltag = gamw_est .- gamw
    maxdeltagamw = [0.0]

    gamw, omega = dt.relax_gamw!(
        gamw, deltag_prev, deltag, maxdeltagamw; nrf=0.5, btw=0.6, pfw=1.2, test=true
    )

    @test isapprox(maxdeltagamw[1], DGTHMAX)
    @test isapprox(omega[], RLXG[end], atol=1e-6)

    for (g, gdfdc) in zip(gamw, GTHupdate' .* 1.0)
        @test isapprox(g, gdfdc, atol=1e-6)
    end
end
