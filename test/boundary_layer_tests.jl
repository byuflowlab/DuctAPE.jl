using Test

@testset "Stagnation Point Identification" begin
    cp_duct = range(-10, 10; step=1)
    u, l = dt.split_at_stagnation_point(cp_duct)
    @test u == 12:21
    @test l == 12:-1:1
end

@testset "Radius of Curvature Calculation" begin
    s = range(1, 10)
    controlpoint = [s'; s' .^ 2]
    ss = 4.0
    R = dt.calc_radius_of_curvature(s, controlpoint, ss)
    @test isapprox(R, 1.0 / 0.003816453371976)
end

@testset "Step Setter" begin
    N = 10
    first_step_size = 1
    total_length = 100
    steps = dt.set_bl_steps(N, first_step_size, total_length)

    @test steps[1] == first_step_size
    @test steps[end] == total_length
    @test length(steps) == N
end

@testset "Raw Arc Length Calcuations" begin
    panel_lengths = ones(10)
    bl_ids = 1:6
    s = dt.arc_lengths_from_panel_lengths(panel_lengths, bl_ids)

    @test s == 0:5
end

@testset "Boundary Layer Initialization" begin
    x = 1.0
    Rex = 1.0
    @test 0.036 == dt.d2_init(x, Rex)

    @test dt.calc_H12bar0(2.0, 0.0) == 1.0 / (1.0 - 6.55)

    @test dt.calc_CEeq(1.0, 0.0, 1.0, 1.0) == sqrt((1.0 - 0.32) / 1.2 + 0.0001) - 0.01

    s_init = 0.1
    r_init = 1.0
    Ue = 1.0
    M = 0.1
    rhoe = 1.0
    mue = 1e-7

    states, Cf_init, H12_init = dt.initialize_turbulent_boundary_layer_states(
        s_init, r_init, Ue, M, rhoe, mue
    )

    @test isapprox(
        states, [0.00022714464401286951, 1.3836007395668006, 0.047177584216703636]
    )
    @test isapprox(Cf_init, 0.0035818798451037327)
    @test isapprox(H12_init, 1.3883679410459342)
end

# @testset "Cf0" begin
#     Red2 = 1.0
#     M = 1.0
#     Cf0 = dt.calc_Cf0(Red2, M; hardness=50)
#     @test isapprox(Cf0, (0.01013 / (1.05 - 1.02) - 0.00075) / dt.Fc(M))
#     Red2 = 100.0
#     M = 100.0
#     Cf0 = dt.calc_Cf0(Red2, M; hardness=50)
#     @test isapprox(Cf0, (0.01013 / (log(10, dt.FR(M) * Red2) - 1.02) - 0.00075) / dt.Fc(M))
# end
