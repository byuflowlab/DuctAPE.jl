@testset "Body System" begin

    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    # x2 = [0.0; 0.5; 1.0]
    # r2 = [0.0; 0.5; 0.0]

    c1 = [x1 r1]
    # c2 = [x2 r2]

    # coordinates = [c1,c2]
    coordinates = [c1]

    # generate panels
    panels = dt.generate_panels(coordinates)

    # initialize LHS
    LHS = dt.vortex_influence_matrix(
        panels.controlpoint, panels.normal, panels.node, panels.influence_length
    )
    LHSraw = copy(LHS)

    # kutta LHS
    dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

    #only the first and last columns should be changed
    @test LHS[:, 2:(end - 1)] == LHSraw[:, 2:(end - 1)]

    Vs = [ones(length(x1)) zeros(length(r1))]

    RHS = dt.freestream_influence_vector(panels.normal, Vs)
    @test all(isapprox.(RHS, [-1.0; 1.0; 1.0; -1.0] * 0.5 * sqrt(2.0)))
end

