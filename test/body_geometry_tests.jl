@testset "Body Geometry" begin

    # check that your input and output Rtip value for the first rotor is the same
    duct_coordinates = [0.0 0.0; 0.5 0.0; 1.0 0.0]
    hub_coordinates = [0.0 0.0; 0.5 0.0; 1.0 0.0]
    Rtip = 1.0
    tip_gaps = [0.1; 0.0]
    rotorzlocs = [0.5; 1.0]

    dc, rt, rh = dt.place_duct(duct_coordinates, hub_coordinates, Rtip, tip_gaps, rotorzlocs)

    @test isapprox(Rtip, rt[1])
    @test all(dc[:,2] .== 1.1)
end
