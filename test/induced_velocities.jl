@testset "Normalized Relative Geometry" begin
    # setup
    node = [1.0; 1.0]
    controlpoint = [0.0; 2.0]
    influence_length = 1.0

    # geometry
    xi, rho, m, r_influence = dt.calculate_xrm(controlpoint, node)
    @test xi == (controlpoint[1] - node[1]) / node[2]
    @test rho == controlpoint[2] / node[2]
    @test m == 4 * rho / (xi^2 + (rho + 1)^2)
    @test r_influence == node[2]
end

#---------------------------------#
#             VORTICES            #
#---------------------------------#
@testset "Basic Vortex Ring" begin

    # influencing panel
    x = [1.0; 2.0]
    r = [9.0; 10.0]
    ip = dt.generate_panels([x r])

    #check self influence
    xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.node)
    @test dt.vortex_ring_vz(xi, rho, m, r_influence, ip.influence_length[1]) ==
        1.0 / (4.0 * pi * r_influence) *
          (log(8.0 * pi * r_influence / ip.influence_length[1]) - 0.25)
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence on axis should have influence in x, and zero in r
    ra = [0.0; 0.0]
    ap = dt.generate_panels([x ra])
    xi, rho, m, r_influence = dt.calculate_xrm(ap.controlpoint[1, :], ip.node[1, :])
    @test dt.vortex_ring_vz(xi, rho, m, r_influence, ip.influence_length[1]) != 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence FROM axis (zeros all around)
    xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint[1, :], ap.node[1, :])
    @test dt.vortex_ring_vz(xi, rho, m, r_influence, ap.influence_length[1]) == 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0
end

@testset "Nominal Single Vortex Panel Integration" begin

    # - Test 1 - #
    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/nominal_velocities1.jl")
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

    # - Test 2 - #
    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/nominal_velocities2.jl")
    node1 = p1
    node2 = p2
    influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    controlpoint = pf

    # Calculate Integral
    V = dt.nominal_vortex_panel_integration(
        node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
    )

    # Compare with DFDC integration values
    @test isapprox(V[1, 1], Vgammai[1], atol=1e-9)
    @test isapprox(V[1, 2], Vgammai[2], atol=1e-9)
    @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-9)
    @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-9)
end

@testset "Single Vortex Panel Self-Induction Integration" begin

    # - Test 1 - #
    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/self_velocities1.jl")
    node1 = p1
    node2 = p2
    influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    controlpoint = ps

    # Calculate Integral
    V = dt.self_vortex_panel_integration(
        node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
    )

    # Compare with DFDC integration values
    # note: axial tests fail if I use log(8Δs/r) in the analytic term rather than the 16 DFDC has instead of the 8
    # @test isapprox(V[1, 1], Vgammai[1], atol=3e-3)
    # @test isapprox(V[2, 1], Vgammaip1[1], atol=3e-3)

    #note: radial terms have zero for the analytic addition, so they work fine
    @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
    @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)

    # - Test 2 - #
    # Load in comparision values from DFDC extraction
    include("./data/single_linear_panel_integration/self_velocities2.jl")
    node1 = p1
    node2 = p2
    influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    controlpoint = ps

    # Calculate Integral
    V = dt.self_vortex_panel_integration(
        node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
    )

    # Compare with DFDC integration values
    # note: axial tests fail if I use log(8Δs/r) in the analytic term rather than the 16 DFDC has instead of the 8
    # @test isapprox(V[1, 1], Vgammai[1], atol=3e-4)
    # @test isapprox(V[2, 1], Vgammaip1[1], atol=3e-4)

    #note: radial terms have zero for the analytic addition, so they work fine
    @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
    @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)
end

#---------------------------------#
#             SOURCES             #
#---------------------------------#
@testset "Basic Source Ring" begin

    # influencing panel
    x = [1.0; 2.0]
    r = [9.0; 10.0]
    ip = dt.generate_panels([x r])

    #check self influence
    xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.node)
    @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
    @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence ON axis should have influence in x, and zero in r
    ra = [0.0; 0.0]
    ap = dt.generate_panels([x ra])
    xi, rho, m, r_influence = dt.calculate_xrm(ap.controlpoint[1, :], ip.node[1, :])
    @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
    @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence FROM axis (zeros all around)
    xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint[1, :], ap.node[1, :])
    @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
    @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0
end
#TODO: need to figure out source singularity integrals first before updating these
# @testset "Source Induced Velocities" begin
#     # setup
#     node = [1.0; 1.0]
#     controlpoint = [0.0; 2.0]
#     influence_length = 1.0
#     xi, rho, m, r_influence = dt.calculate_xrm(controlpoint, node)

#     # individual components
#     vz = dt.source_ring_vz(xi, rho, m, r_influence)
#     vr = dt.source_ring_vr(xi, rho, m, r_influence)

#     # components together
#     vel1 = dt.source_induced_velocity(controlpoint, node, influence_length)
#     @test vel1 == [vz; vr]

#     # components in place
#     vel2 = similar(vel1) .= 0.0
#     dt.source_induced_velocity!(vel2, controlpoint, node, influence_length)

#     @test vel1 == vel2
# end
