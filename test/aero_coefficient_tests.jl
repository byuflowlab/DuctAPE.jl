@testset "One-way Mesh Tests" begin
    @testset "Vortex Meshes" begin

        # - Single Bodies - #
        influence_panels = [(
            panel_center=[0.25 1.0; 0.75 1.0], panel_length=[0.5; 0.5], npanels=2
        )]
        affect_panels = [(
            panel_center=[0.25 2.0; 0.75 2.0], panel_length=[0.5; 0.5], npanels=2
        )]

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="vortex"
        )

        @test mesh.x == [0.0 -0.5; 0.5 0.0]
        @test mesh.r == 2.0 * ones(2, 2)
        @test mesh.mesh2panel_influence == [1; 2]
        @test mesh.mesh2panel_affect == [1; 2]

        # - Multi Body Affect - #
        affect_panels = [
            (panel_center=[0.25 2.0; 0.75 2.0], panel_length=[0.5; 0.5], npanels=2) for
            i in 1:2
        ]

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="vortex"
        )

        @test mesh.x == [0.0 -0.5; 0.5 0.0; 0.0 -0.5; 0.5 0.0]
        @test mesh.r == 2.0 * ones(4, 2)
        @test mesh.mesh2panel_influence == [1; 2]
        @test mesh.mesh2panel_affect == [1; 2; 1; 2]
    end

    @testset "Source Meshes" begin

        # - Single Bodies - #
        influence_panels = [(
            panel_center=[0.25 1.0; 0.75 1.0], panel_length=[0.5; 0.5], npanels=2
        )]
        affect_panels = [(
            panel_center=[0.25 2.0; 0.75 2.0], panel_length=[0.5; 0.5], npanels=2
        )]

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="source"
        )

        @test mesh.x == [0.0 -0.25; 0.25 0.0]
        @test mesh.r == 2.0 * ones(2, 2)
        @test mesh.mesh2panel_influence == [1; 2]
        @test mesh.mesh2panel_affect == [1; 2]

        # - Multi Body Affect - #
        affect_panels = [
            (panel_center=[0.25 2.0; 0.75 2.0], panel_length=[0.5; 0.5], npanels=2) for
            i in 1:2
        ]

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="source"
        )

        @test mesh.x == [0.0 -0.25; 0.25 0.0; 0.0 -0.25; 0.25 0.0]
        @test mesh.r == 2.0 * ones(4, 2)
        @test mesh.mesh2panel_influence == [1; 2]
        @test mesh.mesh2panel_affect == [1; 2; 1; 2]
    end
end

@testset "One-way Coefficient Matrix Tests" begin
    @testset "Self-induction " begin

        # - Self-inducement Tests - #
        method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
        influence_panels = [ff.generate_panels(method, [0.0 1.0; 1.0 1.0])]
        affect_panels = [ff.generate_panels(method, [0.0 1.0; 1.0 1.0])]

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="vortex"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, affect_panels; singularity="vortex"
        )

        cons =
            4.0 * pi * influence_panels[1].panel_center[mesh.mesh2panel_influence[1], 2] /
            influence_panels[1].panel_length[mesh.mesh2panel_influence[1]]
        @test aic[1] ==
            -0.5 - influence_panels[1].panel_curvature[mesh.mesh2panel_influence[1]] -
              (log(2.0 * cons) - 0.25) / cons *
              cos(influence_panels[1].panel_angle[mesh.mesh2panel_influence[1]])
    end

    @testset "Vortex Coefficients, u only" begin

        # - Single Body - #

        method = ff.AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true])
        influence_panels = ff.generate_panels(method, [0.0 1.0; 0.5 1.0; 1.0 1.0])
        affect_panels = ff.generate_panels(method, [0.0 2.0; 0.5 2.0; 1.0 2.0])

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="vortex"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, affect_panels; singularity="vortex"
        )

        @test size(aic) == (2, 2)

        for i in 1:2
            for j in 1:2
                u = ff.get_u_ring_vortex(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    influence_panels.panel_length[mesh.mesh2panel_influence[j]],
                    mesh.m[i, j],
                )

                @test aic[i, j] == u * 0.5
            end
        end

        # - Multiple Body Affect - #

        mesh = dt.generate_one_way_mesh(
            influence_panels, [affect_panels; affect_panels]; singularity="vortex"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, [affect_panels; affect_panels]; singularity="vortex"
        )

        @test size(aic) == (4, 2)

        for i in 1:4
            for j in 1:2
                u = ff.get_u_ring_vortex(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    influence_panels.panel_length[mesh.mesh2panel_influence[j]],
                    mesh.m[i, j],
                )

                @test aic[i, j] == u * 0.5
            end
        end
    end

    @testset "Vortex Coefficients, v only" begin

        # - Single Body - #

        method = ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [true])
        influence_panels = ff.generate_panels(method, [1.0 0.0; 1.0 0.5; 1.0 1.0])
        affect_panels = ff.generate_panels(method, [2.0 0.0; 2.0 0.5; 2.0 1.0])

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="vortex"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, affect_panels; singularity="vortex"
        )

        @test size(aic) == (2, 2)

        for i in 1:2
            for j in 1:2
                v = ff.get_v_ring_vortex(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test isapprox(aic[i, j], v * 0.5)
            end
        end

        # - Multiple Body Affect - #

        mesh = dt.generate_one_way_mesh(
            influence_panels, [affect_panels; affect_panels]; singularity="vortex"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, [affect_panels; affect_panels]; singularity="vortex"
        )

        @test size(aic) == (4, 2)

        for i in 1:4
            for j in 1:2
                v = ff.get_v_ring_vortex(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test isapprox(aic[i, j], v * 0.5)
            end
        end
    end

    @testset "Source Coefficients, v only" begin

        # - Single Body - #

        method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])
        influence_panels = ff.generate_panels(method, [0.0 1.0; 0.5 1.0; 1.0 1.0])
        affect_panels = ff.generate_panels(method, [0.0 2.0; 0.5 2.0; 1.0 2.0])

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="source"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, affect_panels; singularity="source"
        )

        @test size(aic) == (2, 2)

        for i in 1:2
            for j in 1:2
                v = ff.get_v_ring_source(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test aic[i, j] == v * 0.5
            end
        end

        # - Multiple Body Affect - #

        mesh = dt.generate_one_way_mesh(
            influence_panels, [affect_panels; affect_panels]; singularity="source"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, [affect_panels; affect_panels]; singularity="source"
        )

        @test size(aic) == (4, 2)

        for i in 1:4
            for j in 1:2
                v = ff.get_v_ring_source(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test aic[i, j] == v * 0.5
            end
        end
    end

    @testset "Source Coefficients, u only" begin

        # - Single Body - #

        method = ff.AxisymmetricProblem(Source(Constant()), Neumann(), [true])
        influence_panels = ff.generate_panels(method, [1.0 0.0; 1.0 0.5; 1.0 1.0])
        affect_panels = ff.generate_panels(method, [2.0 0.0; 2.0 0.5; 2.0 1.0])

        mesh = dt.generate_one_way_mesh(
            influence_panels, affect_panels; singularity="source"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, affect_panels; singularity="source"
        )

        @test size(aic) == (2, 2)

        for i in 1:2
            for j in 1:2
                u = ff.get_u_ring_source(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test isapprox(aic[i, j], -u * 0.5)
            end
        end

        # - Multiple Body Affect - #

        mesh = dt.generate_one_way_mesh(
            influence_panels, [affect_panels; affect_panels]; singularity="source"
        )

        aic = dt.assemble_one_way_coefficient_matrix(
            mesh, influence_panels, [affect_panels; affect_panels]; singularity="source"
        )

        @test size(aic) == (4, 2)

        for i in 1:4
            for j in 1:2
                u = ff.get_u_ring_source(
                    mesh.x[i, j],
                    mesh.r[i, j],
                    influence_panels.panel_center[mesh.mesh2panel_influence[j], 2],
                    mesh.m[i, j],
                )

                @test isapprox(aic[i, j], -u * 0.5)
            end
        end
    end
end
