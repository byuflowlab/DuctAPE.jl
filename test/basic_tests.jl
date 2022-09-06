
@testset "Basic Geometry" begin

    #test basic duct/hub geometry
    outerwallx = [0; 0.5; 1]
    innerwallx = outerwallx
    outerwallr = [1; 1.5; 1]
    innerwallr = [1.0; 1; 1]
    hubx = outerwallx
    hubr = [0; 0; 0.0]

    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    @test ductgeometry.wallinnerxcoordinates == innerwallx
    @test ductgeometry.wallinnerrcoordinates == innerwallr
    @test ductgeometry.wallouterxcoordinates == outerwallx
    @test ductgeometry.wallouterrcoordinates == outerwallr
    @test ductgeometry.hubxcoordinates == hubx
    @test ductgeometry.hubrcoordinates == hubr
    @test ductgeometry.LEx == 0.0
    @test ductgeometry.TEx == 1.0
    @test ductgeometry.chord == 1.0
    @test ductgeometry.wallbluntTE == false
    @test ductgeometry.hubbluntTE == false

    @test ductsplines.wallinnerspline(0.5) == 1.0
    @test ductsplines.wallouterspline(0.5) == 1.5
    @test ductsplines.hubspline(0.5) == 0.0
end

@testset "Basic Wake Grid Geometry" begin

    #setup
    outerwallx = [0; 0.5; 1]
    innerwallx = outerwallx
    outerwallr = [1; 1.5; 1]
    innerwallr = [1.0; 1; 1]
    hubx = outerwallx
    hubr = [0; 0; 0.0]

    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    #test grid point generation
    grid_options = DuctTAPE.defineGridOptions(
        3; inlet_length=0.0, wake_expansion_factor=1.0
    )

    @test grid_options.num_radial_stations == 3
    @test grid_options.inlet_length == 0.0
    @test grid_options.wake_length == 1.0
    @test grid_options.wake_expansion_factor == 1.0

    x_grid_points, r_grid_points, nx, nr, wallTEidx, hubTEidx, rotoridxs, wall_xstations, hub_xstations = DuctTAPE.generate_grid_points(
        ductgeometry, ductsplines, [], grid_options
    )

    @test x_grid_points == [0.0 0.0 0.0; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test nx == 4
    @test nr == 3
    @test wallTEidx == 2
    @test hubTEidx == 2
    @test wall_xstations == [0.0; 1.0]
    @test hub_xstations == [0.0; 1.0]

    #test grid relaxation
    xr, rr = DuctTAPE.relax_grid(x_grid_points, r_grid_points, nx, nr)

    #test everything together
    grid = DuctTAPE.initialize_wakegrid(ductgeometry, ductsplines, [], grid_options)

    @test xr == x_grid_points
    @test rr == r_grid_points

    @test grid.x_grid_points == [0.0 0.0 0.0; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test grid.r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test grid.nx == 4
    @test grid.nr == 3
    @test grid.wallTEidx == 2
    @test grid.hubTEidx == 2
    @test grid.wall_xstations == [0.0; 1.0]
    @test grid.hub_xstations == [0.0; 1.0]
end
