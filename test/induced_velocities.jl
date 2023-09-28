@testset "Basic Vortex Ring" begin

    # influencing panel
    x = [1.0; 2.0]
    r = [9.0; 10.0]
    ip = dt.generate_panels([x r])

    #check self influence
    xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.node)
    @test dt.vortex_ring_vz(xi, rho, m, r_influence, ip.influence_length[1]) ==
        -1.0 / (4.0 * pi * r_influence) * (log(8.0 * pi * r_influence / ip.influence_length[1]) - 0.25)
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

@testset "Vortex Induced Velocities" begin
    # setup
    node = [1.0; 1.0]
    controlpoint = [0.0; 2.0]
    influence_length = 1.0
    xi, rho, m, r_influence = dt.calculate_xrm(controlpoint, node)

    # individual components
    vz = dt.vortex_ring_vz(xi, rho, m, r_influence, influence_length)
    vr = dt.vortex_ring_vr(xi, rho, m, r_influence)

    # components together
    vel1 = dt.vortex_induced_velocity(controlpoint, node, influence_length)
    @test vel1 == [vz; vr]

    # components in place
    vel2 = similar(vel1) .= 0.0
    dt.vortex_induced_velocity!(vel2, controlpoint, node, influence_length)

    @test vel1 == vel2
end

@testset "Source Induced Velocities" begin
    # setup
    node = [1.0; 1.0]
    controlpoint = [0.0; 2.0]
    influence_length = 1.0
    xi, rho, m, r_influence = dt.calculate_xrm(controlpoint, node)

    # individual components
    vz = dt.source_ring_vz(xi, rho, m, r_influence)
    vr = dt.source_ring_vr(xi, rho, m, r_influence)

    # components together
    vel1 = dt.source_induced_velocity(controlpoint, node, influence_length)
    @test vel1 == [vz; vr]

    # components in place
    vel2 = similar(vel1) .= 0.0
    dt.source_induced_velocity!(vel2, controlpoint, node, influence_length)

    @test vel1 == vel2
end

@testset "Vortex Coefficient Functions" begin
    controlpoints = [1.0 1.0; 2.0 1.0]
    nodes = [1.0 1.0; 1.0 2.0]
    influence_lengths = 2.0 * ones(2)
    strengths = 2.0 * ones(2)

    #TODO: set up control points(s) and node(s), and get individual values to test this first one, then the others are all relative to that.
    v = zeros(2, 2, 2)

    for (ic, cp) in enumerate(eachrow(controlpoints))
        for (in, np) in enumerate(eachrow(nodes))
            dt.vortex_induced_velocity!(
                view(v, ic, in, :), cp, np, influence_lengths[in], strengths[in]
            )
        end
    end

    AICcomp1 = dt.influencefromvortices(controlpoints, nodes, influence_lengths, strengths)
    @test AICcomp1[:, :, 1] == v[:,:,1]
    @test AICcomp1[:, :, 2] == v[:,:,2]

    AICcomp2 = similar(AICcomp1) .= 0.0
    dt.influencefromvortices!(AICcomp2, controlpoints, nodes, influence_lengths, strengths)
    @test AICcomp1 == AICcomp2

    V1 = dt.vfromvortices(controlpoints, nodes, influence_lengths, strengths)
    @test reduce(vcat, V1) == reduce(vcat, sum(AICcomp2; dims=2))

    V2 = similar(V1) .= 0.0
    dt.vfromvortices!(V2, controlpoints, nodes, influence_lengths, strengths)
    @test V2 == V1
end
