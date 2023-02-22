@testset "Basic Induced Velocity Tests" begin
    @testset "Single Vortex Tests" begin
        # set up an influence "panel" (think of this as a vortex ring, which it is...)
        influence_panel = (panel_center=[10 10], panel_length=[1.0], npanels=1)

        # set up 2 circles of affect "panels" around influence panel.
        affect_panels = (
            panel_center=[12 10; 10 8; 8 10; 10 12; 12 8; 8 8; 8 12; 12 12; 10 10] * 1.0,
            panel_length=ones(9),
            npanels=9,
        )

        # get mesh
        mesh = dt.generate_one_way_mesh([influence_panel], [affect_panels])

        # Get Induced Velocities for unit strength
        vxs = zeros(9)
        vrs = zeros(9)
        for i in 1:9
            vxs[i], vrs[i] = dt.calculate_ring_vortex_influence_off_body(
                affect_panels, influence_panel, mesh, i, 1
            )
       end

        # Test that all the directions are correct according to your chosen convention (positive clockwise)
        @test sign(vxs[1]) == -1
        @test sign(vxs[2]) == -1
        @test sign(vxs[3]) == -1
        @test sign(vxs[4]) == 1
        @test sign(vxs[5]) == -1
        @test sign(vxs[6]) == -1
        @test sign(vxs[7]) == 1
        @test sign(vxs[8]) == 1
        @test sign(vrs[1]) == -1
        @test vrs[2] == 0.0
        @test sign(vrs[3]) == 1
        @test vrs[4] == 0.0
        @test sign(vrs[5]) == -1
        @test sign(vrs[6]) == 1
        @test sign(vrs[7]) == 1
        @test sign(vrs[8]) == -1
    end

    @testset "Single Source Tests" begin
        # set up an influence "panel" (think of this as a source ring, which it is...)

        # set up 2 circles of affect "panels" around influence panel.

        # Test that the magnitudes are the same at the same radii

        # Test that all the directions are correct according to your chosen convention (positive outward)
    end
end
