@testset "Linear System Assembly" begin

    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    # - Test 1 - #
    # closed hub will have 2 prescribed nodes
    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.0]

    c1 = [x1 r1]
    c2 = [x2 r2]

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates)

    # Generate LHS matrix and RHS vector
    AICn, _ = dt.vortex_aic_boundary_on_boundary(
        panels.controlpoint,
        panels.normal,
        panels.tangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    dt.add_te_gap_aic!(
        AICn,
        AICt,
        panels.controlpoint,
        panels.normal,
        panels.tangent,
        panels.tenode,
        panels.teinfluence_length,
        panels.tendotn,
        panels.tencrossn,
        panels.teadjnodeidxs,
    )

    AICpcp, unused = dt.vortex_aic_boundary_on_field(
        panels.itcontrolpoint,
        panels.itnormal,
        panels.ittangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    dt.add_te_gap_aic!(
        AICpcp,
        unused,
        panels.itcontrolpoint,
        panels.itnormal,
        panels.ittangent,
        panels.tenode,
        panels.teinfluence_length,
        panels.tendotn,
        panels.tencrossn,
        panels.teadjnodeidxs,
    )

    Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    vdnb = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel))
    vdnpcp = dt.freestream_influence_vector(
        panels.itnormal, repeat(Vs, size(panels.itcontrolpoint, 1))
    )

    LHS = dt.assemble_lhs_matrix(AICn, AICpcp, panels; dummyval=1.0)
    RHS = dt.assemble_rhs_matrix(vdnb, vdnpcp, panels)

    #check that added columns are correct
    @test LHS[:, 9:10] == [
        1.0 0.0
        1.0 0.0
        1.0 0.0
        1.0 0.0
        0.0 1.0
        0.0 1.0
        0.0 0.0
        0.0 0.0
        0.0 0.0
        0.0 0.0
    ]

    #check that prescribed rows are correct
    @test LHS[8:10, :] == [
        1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0# kutta
        0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0# prescribed node 1
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0# prescribed node 2
    ]

    # check that zeros are in correct places on RHS
    @test all(RHS[1:7] .!= 0.0)
    @test all(RHS[8:10] .== 0.0)

    # - Test 2 - #
    # open hub will have 1 prescribed node and 1 internal point
    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.25]

    c1 = [x1 r1]
    c2 = [x2 r2]

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates)

    # Generate LHS matrix and RHS vector
    AICn, _ = dt.vortex_aic_boundary_on_boundary(
        panels.controlpoint,
        panels.normal,
        panels.tangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    AICpcp, _ = dt.vortex_aic_boundary_on_field(
        panels.itcontrolpoint,
        panels.itnormal,
        panels.ittangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    vdnb = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel))
    vdnpcp = dt.freestream_influence_vector(
        panels.itnormal, repeat(Vs, size(panels.itcontrolpoint, 1))
    )

    LHS = dt.assemble_lhs_matrix(AICn, AICpcp, panels; dummyval=1.0)
    RHS = dt.assemble_rhs_matrix(vdnb, vdnpcp, panels)

    #check that added columns are correct
    @test LHS[:, 9:10] == [
        1.0 0.0
        1.0 0.0
        1.0 0.0
        1.0 0.0
        0.0 1.0
        0.0 1.0
        0.0 0.0
        0.0 0.0
        0.0 0.0
        0.0 0.0
    ]

    #check that prescribed rows are correct
    @test LHS[8:9, :] == [
        1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0# kutta
        0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0# prescribed node 1
    ]

    #check that the internal panel influences go used instead of the prescribed 2nd panel
    @test !all(LHS[10, :] .!= [0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0])

    # check that zeros are in correct places on RHS
    @test all(RHS[1:7] .!= 0.0)
    @test all(RHS[8:10] .== 0.0)
end

