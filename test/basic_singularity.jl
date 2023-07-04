@testset "Vortex Panel" begin

    # influencing panel
    x = [1.0; 2.0]
    r = [9.0; 10.0]
    ip = dt.generate_panels([x r])

    #check self influence
    xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint, ip.controlpoint)
    @test dt.vortex_ring_vx(xi, rho, m, r_influence, ip.len[1]) ==
        -1.0 / (4.0 * pi * rho) * (log(8.0 * pi * rho / ip.len[1]) - 0.25)
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence on axis should have influence in x, and zero in r
    ra = [0.0; 0.0]
    ap = dt.generate_panels([x ra])
    xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint, ap.controlpoint)
    @test dt.vortex_ring_vx(xi, rho, m, r_influence, ip.len[1]) != 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

    #check influence FROM axis (zeros all around)
    xi, rho, m, r_influence = dt.calculate_xrm(ap.controlpoint, ip.controlpoint)
    @test dt.vortex_ring_vx(xi, rho, m, r_influence, ip.len[1]) == 0.0
    @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

end

