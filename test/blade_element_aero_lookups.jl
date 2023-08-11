@testset "Blade Element Aero Lookups" begin
    @testset "DFDC Parameterization" begin
        dfdcparam = dt.DFDCairfoil(
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
                        filenames[i, j, k, ell] = "test/data/cas$i$j$k$ell.dat"
                    end
                end
            end
        end

        cas = dt.StaggerInflowReMachSolidityCAS(filenames)

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

    # TODO: ccblade airfoil capabilities should be working, but want to unit test at some point.
    @testset "CCBlade Airfoil" begin end
end
