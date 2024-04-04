println("\nPANEL INITIALIZATION TESTS")
@testset "Panel Geometry" begin

    #=
    Note:
    Currently we assume the following about the inputs:
    for bodies, duct coordinates are given first
    for bodies, coordinates are given from TE to TE clockwise
    =#

    # single panel
    coordinates = [[0.0 1.0; 1.0 1.0]]
    panels = dt.generate_panels(coordinates)

    @test panels.controlpoint[:, 1] == [0.5; 1.0]
    @test panels.normal[:, 1] == [0.0; 1.0]
    @test panels.tangent[:, 1] == [1.0; 0.0]
    @test panels.node == coordinates[1]
    @test panels.nodemap == [1; 2;;]
    #note: for single panel, second node remains at initial zeros
    @test panels.endnodes[1, :, :] == [0.0 1.0; 0.0 0.0]
    #note: for single panel, second node remains at initial ones
    @test panels.endnodeidxs == [1; 1;;]
    @test panels.endpanelidxs == [1; 1;;]
    @test panels.influence_length == [1.0]
    @test panels.totpanel == 1
    @test panels.totnode == 2
    @test panels.npanel == [1]
    @test panels.nnode == [2]
    @test panels.tenode == [0.0 0.0;;; 0.0 0.0] # this test doesn't make sense in this case, but is present for completeness and to let me know if this behavior of the second index staying zeros changes.
    @test panels.tenormal == [1.0; 0.0;;]
    @test panels.teadjnodeidxs == [1; 1;;] #for single panel, second node remains at inital ones
    @test panels.tendotn == [0.0; 0.0;;]
    @test panels.tencrossn == -[1.0; 1.0;;]

    ## more panels and bodies
    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.0]

    c1 = [x1'; r1']
    c2 = [x2'; r2']

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates)

    # tests
    @test panels.influence_length == [
        sqrt(2) / 2
        sqrt(2) / 2
        sqrt(2) / 2
        sqrt(2) / 2
        sqrt(2) / 2
        sqrt(2) / 2
    ]
    @test panels.controlpoint == [
        0.75 1.75
        0.25 1.75
        0.25 2.25
        0.75 2.25
        0.25 0.25
        0.75 0.25
    ]'
    testnodes = reduce(hcat, coordinates)
    @test panels.node == testnodes
    @test panels.nodemap == [1 2; 2 3; 3 4; 4 5; 6 7; 7 8]'
    testnorm = sqrt(2) / 2 .* [1 -1; -1 -1; -1 1; 1 1; -1 1; 1 1]'
    for (pn, tn) in zip(eachcol(panels.normal), eachcol(testnorm))
        @test isapprox(pn, tn)
    end
    testtan = sqrt(2) / 2 .* [-1 -1; -1 1; 1 1; 1 -1; 1 1; 1 -1]'
    for (pn, tt) in zip(eachcol(panels.tangent), eachcol(testtan))
        @test isapprox(pn, tt)
    end
    @test panels.endnodes[1, :, :] == [1.0 2.0; 1.0 2.0]
    @test panels.endnodes[2, :, :] == [0.0 0.0; 1.0 0.0]
    @test panels.endnodeidxs[:, 1] == [1; 5]
    @test panels.endnodeidxs[:, 2] == [6; 8]
    @test panels.endpanelidxs[:, 1] == [1; 4]
    @test panels.endpanelidxs[:, 2] == [5; 6]
    @test panels.npanel == [4; 2]
    @test panels.nnode == [5; 3]
    @test panels.itcontrolpoint == [0.9646446609406726 2.0; 0.5 0.03535533905932738]'
    @test panels.itnormal == [1.0 0.0; 0.0 -1.0]'
    @test panels.tenode[1, :, :] == [1.0 2.0; 1.0 2.0]
    @test panels.tenode[2, :, :] == [1.0 0.0; 1.0 0.0]
    @test panels.tenormal[:, 1] == [0.0, 0.0] # doesn't make sense for sharp trailing edge, but as long as it's not something wacky, we're good.
    @test panels.tenormal[:, 2] == [1.0, 0.0] # note this is hardcoded in this case.
    @test panels.teadjnodeidxs[:, 1] == [5, 1]
    @test panels.teadjnodeidxs[:, 2] == [8, 8]
    @test panels.tendotn[:, 1] == [0.0, 0.0] # again, this and the cross product don't make sense for a sharp trailing edge, but need to check it's not behaving unstably
    @test all(isapprox.(panels.tendotn[:, 2], [sqrt(2) / 2, 0.0]))
    @test panels.tencrossn[:, 1] == [0.0, 0.0]
    @test all(isapprox.(panels.tencrossn[:, 2], -[sqrt(2) / 2, sqrt(2) / 2]))
    @test all(isapprox.(panels.teinfluence_length, [0.0, 0.0]))
    @test panels.totnode == 8
    @test panels.totpanel == 6

    # - TE panel specific tests - #
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [1.9; 1.5; 2.0; 2.5; 2.1]

    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.1]

    c1 = [x1'; r1']
    c2 = [x2'; r2']

    coordinates = [c1, c2]

    # generate nominal (wake-like) panels
    panels = dt.generate_panels(coordinates)
    @test panels.tenode[1, :, :] == [1.0 2.1; 1.0 1.9]
    @test panels.tenode[2, :, :] == [1.0 0.1; 1.0 0.0]
    @test panels.tenormal[:, 1] == [1.0, 0.0]
    @test panels.tenormal[:, 2] == [1.0, 0.0] # note this is hardcoded in this case.
    @test panels.teadjnodeidxs[:, 1] == [5, 1]
    @test panels.teadjnodeidxs[:, 2] == [8, 8]
    @test panels.tendotn[1, 1] ==
        dt.dot(panels.tenormal[:, 1], panels.normal[:, panels.endpanelidxs[2, 1]])
    @test isapprox(
        panels.tendotn[1, 2],
        dt.dot(panels.tenormal[:, 1], panels.normal[:, panels.endpanelidxs[1, 1]]),
    )
    @test isapprox(
        panels.tendotn[2, 1],
        dt.dot(panels.tenormal[:, 2], panels.normal[:, panels.endpanelidxs[2, 2]]),
    )
    @test panels.tendotn[2, 2] == 0.0
    @test panels.tencrossn[1, 1] ==
        -dt.cross2mag(panels.tenormal[:, 1], panels.normal[:, panels.endpanelidxs[2, 1]])
    @test panels.tencrossn[2, 1] ==
        -dt.cross2mag(panels.tenormal[:, 1], panels.normal[:, panels.endpanelidxs[1, 1]])
    @test panels.tencrossn[1, 2] ==
        -dt.cross2mag(panels.tenormal[:, 2], panels.normal[:, panels.endpanelidxs[2, 2]])
    @test panels.tencrossn[2, 2] ==
        -dt.cross2mag(panels.tenormal[:, 2], panels.normal[:, panels.endpanelidxs[2, 2]])
    @test all(isapprox.(panels.teinfluence_length, [0.2, 0.1]))
end
