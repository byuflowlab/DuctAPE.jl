println("\nROTOR AERODYNAMICS TESTS")

@testset "Blade Element Polar Lookups" begin
    @testset "DFDC Parameterization" begin
        dfdcparam = dt.c4b.DFDCairfoil(
            0.0, 1.5, -1.0, 6.28, 0.5, 0.2, 0.012, 0.1, 0.005, 0.0, 200000.0, 0.35, 0.7
        )
        cl, cdrag = dt.lookup_clcd(
            dfdcparam, #inner airfoil
            dfdcparam, #outer airfoil
            0.5, #inner fraction
            10.0, #Wmag_rotor
            0.5, # solidity
            0.5, # stagger
            10 * pi / 180, #alpha (unused for cascade)
            0.0, #inflow
            100000, # reynolds
            0.1, #mach
            341.0, #speed of sound (not used for cascade)
        )

        @test cl == 1.0735559774586954
        @test cdrag == 0.013159959878426141
    end

    @testset "Cascade Lookup" begin
        filenames = Array{String}(undef, 3, 3, 3, 3)
        for i in 1:3
            for j in 1:3
                for k in 1:3
                    for ell in 1:3
                        filenames[i, j, k, ell] = "data/cas$i$j$k$ell.dat"
                    end
                end
            end
        end

        cas = dt.c4b.InReStSoMaCAS(filenames)

        cl, cdrag = dt.lookup_clcd(
            cas, #inner airfoil
            cas, #outer airfoil
            0.5, #inner fraction
            10.0, #Wmag_rotor
            0.5, # solidity
            0.5, # stagger
            10 * pi / 180, #alpha (unused for cascade)
            0.0, #inflow
            100000, # reynolds
            0.1, #mach
            341.0, #speed of sound (not used for cascade)
        )

        @test cl == 0.5
        @test cdrag == 0.05
    end

    @testset "CCBlade Airfoil" begin
        alpha = [-pi; 0.0; pi]
        cl = 2.0 * ones(length(alpha))
        cd = ones(length(alpha))
        af = dt.c4b.AlphaAF(alpha, cl, cd, "TEST")
        @test dt.search_polars(af, pi / 2.0) == (2.0, 1.0)
    end
end

@testset "Blade Element Circulation and Source Strength Calculation Set 1" begin
    alpha = [-pi; 0.0; pi]
    cl = 2.0 * ones(length(alpha))
    cd = ones(length(alpha))
    af = dt.c4b.AlphaAF(alpha, cl, cd, "TEST")

    Vinf = 1.0
    vm = ones(2, 2)
    vtheta = -ones(2, 2)
    blade_elements = [
        (
            num_radial_stations=[2],
            B=1,
            chords=ones(2),
            twists=pi / 2.0 * ones(2),
            stagger=pi .- (pi / 2.0 * ones(2)),
            solidity=ones(2),
            rbe=ones(2),
            Omega=ones(2),
            inner_airfoil=fill(af, 2),
            outer_airfoil=fill(af, 2),
            inner_fraction=[1.0; 1.0],
            is_stator=false,
        ) for i in 1:2
    ]

    freestream = (; Vinf=1.0, asound=343.0, rhoinf=1.225, muinf=1.81e-5)

    Gamma = similar(vm) .= 0
    Sigma = similar(vm, 3, 2) .= 0

    W = sqrt.(vm .^ 2 .+ vtheta .^ 2)

    dt.calculate_gamma_sigma!(Gamma, Sigma, blade_elements, vm, vtheta, W, freestream)

    @test all(Gamma .== W)
    @test all(Sigma .== 1.0 / (4.0 * pi) * W[1] * ones(3, 2))

    Gamma, Sigma = dt.calculate_gamma_sigma(blade_elements, vm, vtheta, W, freestream)
    @test all(Gamma .== W)
    @test all(Sigma .== 1.0 / (4.0 * pi) * W[1] * ones(3, 2))

end
