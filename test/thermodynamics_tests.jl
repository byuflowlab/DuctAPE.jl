@testset "Thermodynamic Relations" begin
    @testset "Standard Atmosphere" begin

        # - SI Units - #
        T, P, rho, mu = dt.standard_atmosphere(0.0)
        @test isapprox(T, 14.99, atol=1e-12)
        @test isapprox(P, 101400.93090454886, atol=1e-12)
        @test isapprox(rho, 1.22661378741097, atol=1e-12)
        @test isapprox(mu, 1.7889517587366847e-5, atol=1e-12)

        # - Imperial Units - #
        T, P, rho, mu = dt.standard_atmosphere(dt.Imperial(), 0.0)
        @test isapprox(T, 58.982, atol=1e-12)
        @test isapprox(P, 2117.8024776319244, atol=1e-12)
        @test isapprox(rho, 0.002380023263989253, atol=1e-12)
        @test isapprox(mu, 3.7363034244069305e-7, atol=1e-12)
    end

    @testset "OperatingPoint" begin

        # - SI Units - #
        si = dt.OperatingPoint(1.0, 1.0)
        @test isapprox(si.Ttot[], 14.990025904158804, atol=1e-12)
        @test isapprox(si.Ptot[], 101401.54421276736, atol=1e-12)
        @test isapprox(si.asound[], 340.1974608958744, atol=1e-12)
        @test isapprox(si.Minf[], si.Vinf[] / si.asound[], atol=1e-12)

        # - Imperial Units - #
        imp = dt.OperatingPoint(dt.Imperial(), 1.0, 1.0)
        @test isapprox(imp.Ttot[], 58.98200946928549, atol=1e-12)
        @test isapprox(imp.Ptot[], 2117.8036676437946, atol=1e-12)
        @test isapprox(imp.asound[], 1116.133498437807, atol=1e-12)
        @test isapprox(imp.Minf[], imp.Vinf[] / imp.asound[], atol=1e-12)
    end
end
