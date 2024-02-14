println("\nINDUCED VELOCITY TESTS")

######################################################################
#                                                                    #
#                        RING INDUCED VELOCITIES                     #
#                                                                    #
######################################################################
@testset "Basic Ring Induced Velocities" begin
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

    @testset "Basic Vortex Ring" begin

        # influencing panel
        x = [1.0; 2.0]
        r = [9.0; 10.0]
        ip = dt.generate_panels([x'; r'])

        #check self influence
        xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.node)
        @test dt.vortex_ring_vz(xi, rho, m, r_influence, ip.influence_length[1]) ==
            1.0 / (4.0 * pi * r_influence) *
              (log(8.0 * pi * r_influence / ip.influence_length[1]) - 0.25)
        # @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

        #check influence on axis should have influence in x, and zero in r
        ra = [0.0; 0.0]
        ap = dt.generate_panels([x'; ra'])
        xi, rho, m, r_influence = dt.calculate_xrm(ap.controlpoint[:, 1], ip.node[:, 1])
        @test dt.vortex_ring_vz(xi, rho, m, r_influence, ip.influence_length[1]) != 0.0
        @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0

        #check influence FROM axis (zeros all around)
        xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint[:, 1], ap.node[:, 1])
        @test dt.vortex_ring_vz(xi, rho, m, r_influence, ap.influence_length[1]) == 0.0
        @test dt.vortex_ring_vr(xi, rho, m, r_influence) == 0.0
    end

    @testset "Basic Source Ring" begin

        # influencing panel
        x = [1.0; 2.0]
        r = [9.0; 10.0]
        ip = dt.generate_panels([x'; r'])

        #check self influence
        xi, rho, m, r_influence = dt.calculate_xrm(ip.node, ip.node)
        @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
        @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0

        #check influence ON axis should have influence in x, and zero in r
        ra = [0.0; 0.0]
        ap = dt.generate_panels([x'; ra'])
        xi, rho, m, r_influence = dt.calculate_xrm(ap.controlpoint[:, 1], ip.node[:, 1])
        @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
        @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0

        #check influence FROM axis (zeros all around)
        xi, rho, m, r_influence = dt.calculate_xrm(ip.controlpoint[:, 1], ap.node[:, 1])
        @test dt.source_ring_vz(xi, rho, m, r_influence) == 0.0
        @test dt.source_ring_vr(xi, rho, m, r_influence) == 0.0
    end
end

######################################################################
#                                                                    #
#                       PANEL INDUCED VELOCITIES                     #
#                                                                    #
######################################################################
@testset "Single Linear Panel Induced Velocities" begin

    #---------------------------------#
    #             VORTICES            #
    #---------------------------------#
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

        # - Test 3 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities3.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_vortex_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)
        @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)

        # - Test 4 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities4.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_vortex_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)
        @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)
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
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)

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
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)

        #note: radial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)

        # - Test 3 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities3.jl")
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
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)

        #note: radial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)

        # - Test 4 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities4.jl")
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
        @test isapprox(V[1, 1], Vgammai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vgammaip1[1], atol=1e-5)

        #note: radial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 2], Vgammai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vgammaip1[2], atol=1e-5)
    end

    #---------------------------------#
    #             SOURCES             #
    #---------------------------------#

    @testset "Nominal Single Source Panel Integration" begin

        # - Test 1 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities1.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-6)
        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-6)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-6)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-6)

        # - Test 2 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities2.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-9)
        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-9)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-9)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-9)

        # - Test 3 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities3.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)

        # - Test 4 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/nominal_velocities4.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = pf

        # Calculate Integral
        V = dt.nominal_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)
    end

    @testset "Single Source Panel Self-Induction Integration" begin

        # - Test 1 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities1.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = ps

        # Calculate Integral
        V = dt.self_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        #note: axial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)

        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)

        # - Test 2 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities2.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = ps

        # Calculate Integral
        V = dt.self_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        #note: axial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)

        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)

        # - Test 3 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities3.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = ps

        # Calculate Integral
        V = dt.self_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        #note: axial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)

        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)

        # - Test 4 - #
        # Load in comparision values from DFDC extraction
        include("./data/single_linear_panel_integration/self_velocities4.jl")
        node1 = p1
        node2 = p2
        influence_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        controlpoint = ps

        # Calculate Integral
        V = dt.self_source_panel_integration(
            node1, node2, influence_length, controlpoint; nondimrange=[0.0; 1.0]
        )

        # Compare with DFDC integration values
        #note: axial terms have zero for the analytic addition, so they work fine
        @test isapprox(V[1, 1], Vsigmai[1], atol=1e-5)
        @test isapprox(V[2, 1], Vsigmaip1[1], atol=1e-5)

        @test isapprox(V[1, 2], Vsigmai[2], atol=1e-5)
        @test isapprox(V[2, 2], Vsigmaip1[2], atol=1e-5)
    end
end

######################################################################
#                                                                    #
#                   Multi-PANEL INDUCED VELOCITIES                   #
#                                                                    #
######################################################################
@testset "Multi-Panel, Multi-Target Induced Velocity Matrices" begin
    @testset "Multiple Vortex Panel Induced Velocities on Multiple Points" begin

        # define control points
        controlpoints = [0.0 1.0; 1.0 1.0]'

        # define nodes
        nodes = [-0.5 1.0; 0.5 1.0; 1.0 2.0]'

        # define node map
        nodemap = [1 2; 2 3]'

        # define influence lengths
        influence_lengths = [1.0; 1.0]

        # use unit strengths
        strengths = [1.0 1.0; 1.0 1.0]

        # [cp, n, x/r]
        VEL = dt.induced_velocities_from_vortex_panels_on_points(
            controlpoints, nodes, nodemap, influence_lengths, strengths
        )

        # [vz1 vr1; vz2 vr2]
        vn12cp1 = dt.self_vortex_panel_integration(
            nodes[:, 1], nodes[:, 2], influence_lengths[1], controlpoints[:, 1]
        )

        vn12cp2 = dt.nominal_vortex_panel_integration(
            nodes[:, 1], nodes[:, 2], influence_lengths[1], controlpoints[:, 2]
        )

        vn23cp1 = dt.nominal_vortex_panel_integration(
            nodes[:, 2], nodes[:, 3], influence_lengths[1], controlpoints[:, 1]
        )

        vn23cp2 = dt.nominal_vortex_panel_integration(
            nodes[:, 2], nodes[:, 3], influence_lengths[1], controlpoints[:, 2]
        )

        @test all(VEL[1, 1, :] .== vn12cp1[1, :])
        #note that node connected to 2 panels has influence contributions from both panels
        @test all(VEL[1, 2, :] .== vn12cp1[2, :] .+ vn23cp1[1, :])
        @test all(VEL[1, 3, :] .== vn23cp1[2, :])
        @test all(VEL[2, 1, :] .== vn12cp2[1, :])
        @test all(VEL[2, 2, :] .== vn12cp2[2, :] .+ vn23cp2[1, :])
        @test all(VEL[2, 3, :] .== vn23cp2[2, :])
    end

    @testset "Multiple Source Panel Induced Velocities on Multiple Points" begin

        # define control points
        controlpoints = [0.0 1.0; 1.0 1.0]'

        # define nodes
        nodes = [-0.5 1.0; 0.5 1.0; 1.0 2.0]'

        # define node map
        nodemap = [1 2; 2 3]'

        # define influence lengths
        influence_lengths = [1.0; 1.0]

        # use unit strengths
        strengths = [1.0 1.0; 1.0 1.0]

        # [cp, n, x/r]
        VEL = dt.induced_velocities_from_source_panels_on_points(
            controlpoints, nodes, nodemap, influence_lengths, strengths
        )

        # [vz1 vr1; vz2 vr2]
        vn12cp1 = dt.self_source_panel_integration(
            nodes[:, 1], nodes[:, 2], influence_lengths[1], controlpoints[:, 1]
        )

        vn12cp2 = dt.nominal_source_panel_integration(
            nodes[:, 1], nodes[:, 2], influence_lengths[1], controlpoints[:, 2]
        )

        vn23cp1 = dt.nominal_source_panel_integration(
            nodes[:, 2], nodes[:, 3], influence_lengths[1], controlpoints[:, 1]
        )

        vn23cp2 = dt.nominal_source_panel_integration(
            nodes[:, 2], nodes[:, 3], influence_lengths[1], controlpoints[:, 2]
        )

        @test all(VEL[1, 1, :] .== vn12cp1[1, :])
        @test all(VEL[1, 2, :] .== vn12cp1[2, :] .+ vn23cp1[1, :])
        @test all(VEL[1, 3, :] .== vn23cp1[2, :])
        @test all(VEL[2, 1, :] .== vn12cp2[1, :])
        @test all(VEL[2, 2, :] .== vn12cp2[2, :] .+ vn23cp2[1, :])
        @test all(VEL[2, 3, :] .== vn23cp2[2, :])
    end
end
