println("\nLINEAR SYSTEM ASSEMBLY TESTS")
@testset "Linear System Assembly" begin
    integration_options = (; nominal=dt.GaussLegendre(), singular=dt.GaussLegendre())

    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    # - Test 1 - #
    # closed hub will have 2 prescribed nodes
    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.0]

    c1 = [x1'; r1']
    c2 = [x2'; r2']

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates; isbody=true)

    # Generate LHS matrix and RHS vector
    AICn = dt.vortex_aic_boundary_on_boundary(
        panels.controlpoint,
        panels.normal,
        panels.node,
        panels.nodemap,
        panels.influence_length,
        integration_options,
    )

    dt.add_te_gap_aic!(
        AICn,
        panels.controlpoint,
        panels.normal,
        panels.tenode,
        panels.teinfluence_length,
        panels.tendotn,
        panels.tencrossn,
        panels.teadjnodeidxs,
        integration_options
    )

    AICpcp = dt.vortex_aic_boundary_on_field(
        panels.itcontrolpoint,
        panels.itnormal,
        panels.node,
        panels.nodemap,
        panels.influence_length,
        integration_options
    )

    dt.add_te_gap_aic!(
        AICpcp,
        panels.itcontrolpoint,
        panels.itnormal,
        panels.tenode,
        panels.teinfluence_length,
        panels.tendotn,
        panels.tencrossn,
        panels.teadjnodeidxs,
        integration_options
    )

    Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    vdnb = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel[])')
    vdnpcp = dt.freestream_influence_vector(
        panels.itnormal, repeat(Vs, size(panels.itcontrolpoint, 2))
    )

    LHS = dt.assemble_lhs_matrix(
        AICn,
        AICpcp,
        panels.npanel,
        panels.nnode,
        panels.totpanel[],
        panels.totnode[],
        panels.prescribednodeidxs;
        dummyval=1.0,
    )
    RHS = dt.assemble_rhs_matrix(
        vdnb,
        vdnpcp,
        panels.npanel,
        panels.nnode,
        panels.totpanel[],
        panels.totnode[],
        panels.prescribednodeidxs,
    )

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

    c1 = [x1'; r1']
    c2 = [x2'; r2']

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates)

    # Generate LHS matrix and RHS vector
    AICn = dt.vortex_aic_boundary_on_boundary(
        panels.controlpoint,
        panels.normal,
        panels.node,
        panels.nodemap,
        panels.influence_length,
        integration_options
    )

    AICpcp = dt.vortex_aic_boundary_on_field(
        panels.itcontrolpoint,
        panels.itnormal,
        panels.node,
        panels.nodemap,
        panels.influence_length,
        integration_options
    )

    Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    vdnb = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel[])')
    vdnpcp = dt.freestream_influence_vector(
        panels.itnormal, repeat(Vs, size(panels.itcontrolpoint, 1))
    )

    LHS = dt.assemble_lhs_matrix(
        AICn,
        AICpcp,
        panels.npanel,
        panels.nnode,
        panels.totpanel[],
        panels.totnode[],
        panels.prescribednodeidxs;
        dummyval=1.0,
    )
    RHS = dt.assemble_rhs_matrix(
        vdnb,
        vdnpcp,
        panels.npanel,
        panels.nnode,
        panels.totpanel[],
        panels.totnode[],
        panels.prescribednodeidxs,
    )

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

