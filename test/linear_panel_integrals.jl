@testset "Nominal Single Panel Integration" begin

    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/nominal_velocities.jl")
    node1 = p1
    node2 = p2
    influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    controlpoint = pf

    # Calculate Integral
    V = dt.nominal_vortex_panel_integration(
        node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
    )

    # Compare with DFDC integration values
    @test isapprox(V[1, 1], Vgammai[1], atol=1e-6)
    @test isapprox(V[1, 2], Vgammai[2], atol=1e-6)
    @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-6)
    @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-6)
end

@testset "Single Panel Self-Induction Integration" begin

    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/self_velocities.jl")
    node1 = p1
    node2 = p2
    influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    controlpoint = ps

    # Calculate Integral
    V = dt.self_vortex_panel_integration(
        node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
    )

    # Compare with DFDC integration values
    @test isapprox(V[1, 1], Vgammai[1], atol=1e-6)
    @test isapprox(V[1, 2], Vgammai[2], atol=1e-6)
    @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-6)
    @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-6)
end
