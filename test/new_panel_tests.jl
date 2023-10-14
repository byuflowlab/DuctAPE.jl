@testset "Panel Geometry" begin

    #=
    Note:
    Currently we assume the following about the inputs:
    for bodies, duct coordinates are given first
    for bodies, coordinates are given from TE to TE clockwise
    =#

    # single panel
    coordinates = [0.0 1.0; 1.0 1.0]
    panels = dt.generate_panels(coordinates)

    @test panels.controlpoint[1, :] == [0.5; 1.0]
    @test panels.normal[1, :] == [0.0; 1.0]
    @test panels.tangent[1, :] == [1.0; 0.0]
    @test panels.node == coordinates
    @test panels.nodemap == [1 2]
    #note: for single panel, second node remains at initial zeros
    @test panels.endpoints[1, :, :] == [0.0 1.0; 0.0 0.0]
    #note: for single panel, second node remains at initial ones
    @test panels.endpointidxs == [1 1]
    @test panels.influence_length == [1.0]

    ## more panels and bodies
    # define coordinates
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r1 = [2.0; 1.5; 2.0; 2.5; 2.0]

    x2 = [0.0; 0.5; 1.0]
    r2 = [0.0; 0.5; 0.0]

    c1 = [x1 r1]
    c2 = [x2 r2]

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
    ]
    testnodes = reduce(vcat, coordinates)
    @test panels.node == testnodes
    @test panels.nodemap == [1 2; 2 3; 3 4; 4 5; 6 7; 7 8]
    testnorm = sqrt(2) / 2 .* [1 -1; -1 -1; -1 1; 1 1; -1 1; 1 1]
    for (pn, tn) in zip(eachrow(panels.normal), eachrow(testnorm))
        @test isapprox(pn, tn)
    end
    testtan = sqrt(2) / 2 .* [-1 -1; -1 1; 1 1; 1 -1; 1 1; 1 -1]
    for (pn, tt) in zip(eachrow(panels.tangent), eachrow(testtan))
        @test isapprox(pn, tt)
    end
    @test panels.endpoints[1, :, :] == [1.0 2.0; 1.0 2.0]
    @test panels.endpoints[2, :, :] == [0.0 0.0; 1.0 0.0]
    @test panels.endpointidxs[1, :] == [1; 5]
    @test panels.endpointidxs[2, :] == [6; 8]

end
