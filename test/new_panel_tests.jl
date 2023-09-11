@testset "Panel Geometry" begin

    # single panel
    coordinates = [0.0 1.0; 1.0 1.0]
    panels = dt.generate_panels(coordinates)

    @test panels.controlpoint[1, :] == [0.5; 1.0]
    @test panels.normal[1, :] == [0.0; 1.0]
    @test panels.nodes[1, :, :] == coordinates
    #note: for single panel, second node remains at initial zeros
    @test panels.endpoints[1, :, :] == [0.0 1.0; 0.0 0.0]
    #note: for single panel, second node remains at initial ones
    @test panels.endpointidxs == [1 1]

    ## more panels and bodies
    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.0]

    c1 = [x1 r1]
    c2 = [x2 r2]

    coordinates = [c1, c2]

    # generate panels
    panels = dt.generate_panels(coordinates)

    # tests
    @test panels.endpoints[1, :, :] == [1.0 2.0; 1.0 2.0]
    @test panels.endpoints[2, :, :] == [0.0 0.0; 1.0 0.0]
    @test panels.endpointidxs[1, :] == [1; 4]
    @test panels.endpointidxs[2, :] == [5; 6]
end

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
    LHS = dt.doublet_panel_influence_matrix(panels.nodes, panels)
    LHSraw = copy(LHS)

    # kutta LHS
    dt.body_lhs_kutta!(LHS, panels; tol=1e1 * eps(), verbose=true)

    #only the first and last columns should be changed
    @test LHS[:, 2:(end - 1)] == LHSraw[:, 2:(end - 1)]

    Vs = [ones(length(x1)) zeros(length(r1))]

    RHS = dt.freestream_influence_vector(panels.normal, Vs)
    @test all(isapprox.(RHS, [-1.0; 1.0; 1.0; -1.0] * 0.5 * sqrt(2.0)))
end

@testset "Book-keeping Checks" begin

    # load basic 2-rotor geometry that can be hand counted
    include("data/basic_two_rotor_for_test.jl")

    # get inputs
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        hub_coordinates,
        paneling_constants,
        rotor_parameters,
        freestream,
        reference_parameters,
    )

    # check interface indexing and such
    @test inputs.hubwakeinterfaceid == 1:3
    @test inputs.ductwakeinterfaceid == 15:17
    @test inputs.num_wake_x_panels == 7

    # put in wake discritization function test instead of here. they are used inside the precomputed inputs function, but should not be passed out.
    # @test rotor_indices_in_wake == [1; 3]
    # @test hubTE_index = 4
    # @test ductTE_index = 4
end
