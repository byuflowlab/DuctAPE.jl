@testset "Boundary Layer Functions" begin
    # - Stagnation Point - #
    cp_duct = abs.(range(-10, 10; step=1)) .+ 1.1
    panel_lengths = ones(length(cp_duct))
    s, s_stagnation, lower_length, upper_length = dt.split_at_stagnation_point(
        panel_lengths, cp_duct
    )
    @test s == cumsum([0.0; panel_lengths[1:(end - 1)]])
    @test s_stagnation == 10.0
    @test lower_length == 10.0
    @test upper_length == 10.0

    # @testset "Radius of Curvature Calculation" begin
    s = range(1, 10)
    controlpoint = [s'; s' .^ 2]
    ss = 4.0
    R = dt.calculate_radius_of_curvature(s, controlpoint, ss)
    @test isapprox(R, 1.0 / 0.003816453371976)

    # @testset "Step Setter" begin
    N = 10
    first_step_size = 1
    total_length = 100
    steps = dt.set_boundary_layer_steps(N, first_step_size, total_length)

    @test steps[1] == first_step_size
    @test steps[end] == total_length
    @test length(steps) == N

    # @testset "Raw Arc Length Calcuations" begin
    panel_lengths = ones(10)
    s = dt.arc_lengths_from_panel_lengths(panel_lengths)

    @test s == 0:9

    # @testset "Boundary Layer Initialization" begin
    x = 1.0
    Rex = 1.0
    @test 0.036 == dt.d2_init(x, Rex)

    @test dt.calculate_H12bar0(2.0, 0.0) == 1.0 / (1.0 - 6.55)

    @test dt.calculate_CEeq(1.0, 0.0, 1.0, 1.0) == sqrt((1.0 - 0.32) / 1.2 + 0.0001) - 0.01

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
        states, [0.00022714464401286951, 1.3836007395668004, 0.10980454099439806]
    )
    @test isapprox(Cf_init, 0.003581879845103732)
    @test isapprox(H12_init, 1.3883679410459338)

    # @testset "Schlichting" begin
    @test dt.d2_init(1.0, 1.0) == 0.036

    # - H12bar_init - #
    @test dt.H12bar_init(1.0, 1.0) == dt.calculate_H12bar0(1.0, 1.0)

    # - CE_init - #
    @test dt.CE_init(1.0, 0.1, 0.1) == dt.calculate_CEeq(1.0, 0.1, 1, 0.1)

    # - Fc - #
    @test dt.Fc(0) == 1.0

    # - FR - #
    @test dt.FR(0) == 1.0

    # - Reynolds - #
    @test dt.calculate_Re(rhoe, Ue, 1.0, mue) == rhoe * Ue / mue

    # - H12bar0 - #
    @test dt.calculate_H12bar0(0, 1.0) == 1.0

    # - Cf - #
    @test isapprox(dt.calculate_Cf(2.2, 1.0, 1.0), 0.0, atol=eps())

    # - H12 - #
    @test dt.calculate_H12(0.0, sqrt(5)) == 1.0

    # - Ctau - #
    @test dt.calculate_Ctau(1.0, 0.0, sqrt(10)) == 2.0 * (0.024 + 1.2)

    # - F - #
    @test dt.calculate_F(1.0, 1.0) == (0.02 + 1.0 + 0.8 * 1.0 / 3) / 1.01

    # - H1 - #
    @test dt.calculate_H1(2.0) == 3.15 + 1.72 - 0.01

    # - dH12bardH1 - #
    @test dt.calculate_dH12bardH1(2.0) == -1.0 / (1.72 + 0.02)

    # - Ri - #
    @test dt.calculate_richardson_number(1.0, 1.0, 0.0, 1.0, 1.0) == 2.0 / 3.0 * (1.0 + 0.3)

    # - Secondary Influences - #
    Ri = 1.0
    @test dt.longitudinal_curvature_influence(M, Ri) == 8.014
    @test dt.lateral_strain_influence(ones(6)...) == -8.0 - 1 / 3
    @test dt.dilation_influence(ones(7)...) == 10 + 1 / 3

    # - d2dUedsUeeq0 - #
    @test dt.calculate_d2dUedsUeeq0(ones(4)...) == 0.625

    # - CEeq0 - #
    @test dt.calculate_CEeq0(ones(4)...) == 0.0

    # - Ctaueq0 - #
    @test isapprox(dt.calculate_Ctaueq0(ones(3)...), 1.936)

    # - CEeq - #
    @test isapprox(dt.calculate_CEeq(ones(4)...), 0.6907204085147591)

    # - d2dUedsUeeq - #
    @test dt.calculate_d2dUedsUeeq(ones(4)...) == -0.25
end
