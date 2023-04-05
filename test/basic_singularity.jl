@testset "Vortex Panel" begin

    # influencing panel
    influence_panel = (
        panel_center = [10.0 10.0],
        panel_length = [1.0],
        npanels = 1
    )

    # affected panels (placed in a square around the influencing panel)
    affect_panels = (
        panel_center = [
            12.0 10.0; # right
            10.0  8.0; # bottom
                8.0 10.0; # left
            10.0 12.0; # top
            12.0  8.0; # bottom right
                8.0  8.0; # bottom left
                8.0 12.0; # top left
            12.0 12.0; # top right
            10.0 10.0  # influencing panel
            ],
        panel_length = ones(9),
        npanels = 9,
    )

    # mesh from influencing panel to affected panels
    mesh = DuctTAPE.generate_one_way_mesh([influence_panel], [affect_panels])

    # calculate induced velocities (for unit strength)
    vxs = zeros(9)
    vrs = zeros(9)
    for i in 1:9
        vxs[i], vrs[i] = DuctTAPE.calculate_ring_vortex_influence_off_body(
            affect_panels, influence_panel, mesh, i, 1
        )
    end

    # check that directions are correct (positive clockwise)
    @test vxs[1] < 0  # -0.021258123409045343
    @test vxs[2] < 0  # -0.11285411253919467
    @test vxs[3] < 0  # -0.021258123409045343
    @test vxs[4] > 0  #  0.05324236599283371
    @test vxs[5] < 0  # -0.06495742006793827
    @test vxs[6] < 0  # -0.06495742006793827
    @test vxs[7] > 0  #  0.019406110202735118
    @test vxs[8] > 0  #  0.019406110202735118
    @test vxs[9] == 0 #  0.0

    @test vrs[1] < 0  # -0.07618670901637903
    @test vrs[2] == 0 # -0.0
    @test vrs[3] > 0  #  0.07618670901637903
    @test vrs[4] == 0 # -0.0
    @test vrs[5] < 0  # -0.04053553866884062
    @test vrs[6] > 0  #  0.04053553866884062
    @test vrs[7] > 0  #  0.033980984230275554
    @test vrs[8] < 0  # -0.033980984230275554
    @test vrs[9] == 0 #  0.0

    # # vector field plot
    # using Plots
    # pyplot()
    # strength = 10.0
    # x = affect_panels.panel_center[:,1]
    # y = affect_panels.panel_center[:,2]
    # u = strength*vxs
    # v = strength*vrs
    # quiver(x, y, quiver=(u, v), aspect_ratio=1.0)

end

@testset "Source Panel" begin

    # influencing panel
    influence_panel = (
        panel_center = [10.0 10.0],
        panel_length = [1.0],
        npanels = 1
    )

    # affected panels (placed in a grid around the influencing panel)
    affect_panels = (
        panel_center = [
            12.0 10.0; # right
            10.0  8.0; # bottom
                8.0 10.0; # left
            10.0 12.0; # top
            12.0  8.0; # bottom right
                8.0  8.0; # bottom left
                8.0 12.0; # top left
            12.0 12.0; # top right
            10.0 10.0  # influencing panel
            ],
        panel_length = ones(9),
        npanels = 9,
    )

    # mesh from influencing panel to affected panels
    mesh = DuctTAPE.generate_one_way_mesh([influence_panel], [affect_panels])

    # calculate induced velocities (for unit strength)
    vxs = zeros(9)
    vrs = zeros(9)
    for i in 1:9
        vxs[i], vrs[i] = DuctTAPE.calculate_ring_source_influence_off_body(
            affect_panels, influence_panel, mesh, i, 1
        )
    end

    # check directions (positive outward)
    @test vxs[1] > 0  #  0.08043833369818809
    @test vxs[2] == 0 #  0.0
    @test vxs[3] < 0  # -0.08043833369818809
    @test vxs[4] == 0 #  0.0
    @test vxs[5] > 0  #  0.045419914948660156
    @test vxs[6] < 0  # -0.045419914948660156
    @test vxs[7] < 0  # -0.03689595903578364
    @test vxs[8] > 0  #  0.03689595903578364
    @test vxs[9] == 0 #  0.0

    @test vrs[1] > 0  #  0.021258123409045343
    @test vrs[2] < 0  # -0.06167706566641068
    @test vrs[3] > 0  #  0.021258123409045343
    @test vrs[4] > 0  #  0.09006500512891033
    @test vrs[5] < 0  # -0.020998033549562502
    @test vrs[6] < 0  # -0.020998033549562502
    @test vrs[7] > 0  #  0.05147083306332407
    @test vrs[8] > 0  #  0.05147083306332407
    @test vrs[9] == 0 #  0.0

    # check magnitudes
    @test isapprox(vxs[1], -vxs[3])
    @test isapprox(vrs[1],  vrs[3])
    @test isapprox(vxs[5], -vxs[6])
    @test isapprox(vrs[5],  vrs[6])
    @test isapprox(vxs[7], -vxs[8])
    @test isapprox(vrs[7],  vrs[8])

    # # vector field plot
    # using Plots
    # pyplot()
    # strength = 10.0
    # x = affect_panels.panel_center[:,1]
    # y = affect_panels.panel_center[:,2]
    # u = strength*vxs
    # v = strength*vrs
    # quiver(x, y, quiver=(u, v), aspect_ratio=1.0)

end