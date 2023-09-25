@testset "Vortex Panel" begin

    # influencing panel
    x = [1.0; 2.0]
    r = [9.0; 10.0]
    ip = dt.generate_panels([x r])

    # no self influence needed anymore?
    ##check self influence
    #xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.controlpoint)
    #@test dt.vortex_ring_vx(xi, rho, m, r_influence, ip.influence_length[1]) ==
    #    -1.0 / (4.0 * pi * rho) * (log(8.0 * pi * rho / ip.influence_length[1]) - 0.25)
    #@test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence on axis should have influence in x, and zero in r
    ra = [0.0; 0.0]
    ap = dt.generate_panels([x ra])
    xi, rho, m, r_influence = dt.calculate_xrm(ip.node[1,:], ap.controlpoint[1,:])
    @test dt.vortex_ring_vx(xi, rho, m, r_influence, ip.influence_length[1]) != 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence FROM axis (zeros all around)
    xi, rho, m, r_influence = dt.calculate_xrm(ap.node[1,:], ip.controlpoint[1,:])
    @test dt.vortex_ring_vx(xi, rho, m, r_influence, ap.influence_length[1]) == 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

end

