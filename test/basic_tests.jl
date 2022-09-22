
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

@testset "Basic Wake Grid Geometry: No Rotor" begin

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
    grid = DuctTAPE.generate_wake_grid(ductgeometry, ductsplines, [], grid_options)

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

@testset "Basic Airfoil" begin

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    cl, cd = DuctTAPE.get_clcd(af1, 5.0, 1e6, 0.0, -1.0)
    @test cl == 0.0
    @test cd == 0.0
end

@testset "Basic Rotor: Single" begin

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.25; 0.5; 0.75; 1.0]
    chords = [0.1; 0.1; 0.1; 0.1; 0.1]
    twists = 90.0 * ones(5)                 #test rotor and blade definitions.
    rotor = DuctTAPE.RotorGeometry(
        0.5, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 5), nothing, 1000
    )

    #generate duct geometry (copy from above)
    outerwallx = [0; 0.5; 1]
    innerwallx = outerwallx
    outerwallr = [1; 1.5; 1]
    innerwallr = [1.0; 1; 1]
    hubx = outerwallx
    hubr = [0; 0; 0.0]

    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    #get blade based on grid
    bladedims = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor)

    #test blade geometries.
    @test bladedims.hubr == 0.0
    @test bladedims.tipr == 1.0
    @test bladedims.rdim == rs
    @test bladedims.tdim == twists
    @test isapprox(bladedims.sweptannulus, pi)
    @test isapprox(bladedims.sweptarea, pi)
end

@testset "Basic Rotor: Double" begin

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.25; 0.5; 0.75; 1.0]
    chords = [0.1; 0.1; 0.1; 0.1; 0.1]
    twists = 90.0 * ones(5)                 #test rotor and blade definitions.
    rotor1 = DuctTAPE.RotorGeometry(
        0.25, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 5), nothing, 1000
    )
    rotor2 = DuctTAPE.RotorGeometry(
        0.75, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 5), nothing, 1000
    )

    #generate duct geometry (copy from above)
    outerwallx = [0; 0.5; 1]
    innerwallx = outerwallx
    outerwallr = [1; 1.5; 1]
    innerwallr = [1.0; 1; 1]
    hubx = outerwallx
    hubr = [0; 0; 0.0]

    ductgeometry, ductsplines = DuctTAPE.defineDuctGeometry(
        innerwallx, innerwallr, outerwallx, outerwallr, hubx, hubr
    )

    #get blade based on grid
    bladedims = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor2)

    #test blade geometries.
    @test bladedims.hubr == 0.0
    @test bladedims.tipr == 1.0
    @test bladedims.rdim == rs
    @test bladedims.tdim == twists
    @test isapprox(bladedims.sweptannulus, pi)
    @test isapprox(bladedims.sweptarea, pi)
end

@testset "Basic Grid Geometry: Single Rotor" begin

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

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.5; 1.0]
    chords = [0.1; 0.1; 0.1]
    twists = 90.0 * ones(3)                 #test rotor and blade definitions.
    rotor = DuctTAPE.RotorGeometry(
        0.5, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 3), nothing, 1000
    )

    x_grid_points, r_grid_points, nx, nr, wallTEidx, hubTEidx, rotoridxs, wall_xstations, hub_xstations = DuctTAPE.generate_grid_points(
        ductgeometry, ductsplines, [rotor], grid_options
    )

    @test x_grid_points == [0.5 0.5 0.5; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test nx == 4
    @test nr == 3
    @test wallTEidx == 2
    @test hubTEidx == 2
    @test wall_xstations == [0.0; 0.5; 1.0]
    @test hub_xstations == [0.0; 0.5; 1.0]

    #test grid relaxation
    xr, rr = DuctTAPE.relax_grid(x_grid_points, r_grid_points, nx, nr)

    #test everything together
    grid = DuctTAPE.generate_wake_grid(ductgeometry, ductsplines, [rotor], grid_options)

    @test xr == x_grid_points
    @test rr == r_grid_points

    @test grid.x_grid_points == [0.5 0.5 0.5; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test grid.r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test grid.nx == 4
    @test grid.nr == 3
    @test grid.wallTEidx == 2
    @test grid.hubTEidx == 2
    @test grid.wall_xstations == [0.0; 0.5; 1.0]
    @test grid.hub_xstations == [0.0; 0.5; 1.0]
end

@testset "Basic Grid Geometry: Dual Rotor" begin

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

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.5; 1.0]
    chords = [0.1; 0.1; 0.1]
    twists = 90.0 * ones(3)                 #test rotor and blade definitions.
    rotor1 = DuctTAPE.RotorGeometry(
        0.25, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 3), nothing, 1000
    )
    rotor2 = DuctTAPE.RotorGeometry(
        0.75, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 3), nothing, 1000
    )

    x_grid_points, r_grid_points, nx, nr, wallTEidx, hubTEidx, rotoridxs, wall_xstations, hub_xstations = DuctTAPE.generate_grid_points(
        ductgeometry, ductsplines, [rotor1; rotor2], grid_options
    )

    @test x_grid_points ==
        [0.25 0.25 0.25; 0.75 0.75 0.75; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test nx == 5
    @test nr == 3
    @test wallTEidx == 3
    @test hubTEidx == 3
    @test wall_xstations == [0.0; 0.25; 0.75; 1.0]
    @test hub_xstations == [0.0; 0.25; 0.75; 1.0]

    #test grid relaxation
    xr, rr = DuctTAPE.relax_grid(x_grid_points, r_grid_points, nx, nr)

    #test everything together
    grid = DuctTAPE.generate_wake_grid(
        ductgeometry, ductsplines, [rotor1, rotor2], grid_options
    )

    @test xr == x_grid_points
    @test rr == r_grid_points

    @test x_grid_points ==
        [0.25 0.25 0.25; 0.75 0.75 0.75; 1.0 1.0 1.0; 1.5 1.5 1.5; 2.0 2.0 2.0]
    @test r_grid_points == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
    @test nx == 5
    @test nr == 3
    @test wallTEidx == 3
    @test hubTEidx == 3
    @test wall_xstations == [0.0; 0.25; 0.75; 1.0]
    @test hub_xstations == [0.0; 0.25; 0.75; 1.0]
end

@testset "Basic Paneling Tests" begin

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

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.5; 1.0]
    chords = [0.1; 0.1; 0.1]
    twists = 90.0 * ones(3)                 #test rotor and blade definitions.
    rotor = DuctTAPE.RotorGeometry(
        0.5, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 3), nothing, 1000
    )
    rotors = [rotor]

    #test everything together
    wakegrid = DuctTAPE.generate_wake_grid(
        ductgeometry, ductsplines, [rotor], grid_options
    )

    #generate panels
    panel_geometries = DuctTAPE.generate_panel_geometries(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    wall_panels = panel_geometries.wall_panels
    hub_panels = panel_geometries.hub_panels
    wake_panels = panel_geometries.wake_panels
    rotor_source_panels = panel_geometries.rotor_source_panels

    @test wall_panels.panel_centers ==
        [(0.75, 1.0); (0.25, 1.0); (0.25, 1.25); (0.75, 1.25)]
    @test hub_panels.panel_centers == [(0.25, 0.0); (0.75, 0.0)]
    @test rotor_source_panels[1].panel_centers == [(0.5, 0.25); (0.5, 0.75)]
    #TODO: add wake panel test

    #TODO: need to add test for unit normals once you figure out which direction they are supposed to go.

end

@testset "Basic System Initialization Tests" begin

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

    # generate airfoils for rotor
    datapath = "./data/"
    af1 = DuctTAPE.ccb.AlphaAF(datapath * "basic_af.dat")

    rs = [0.0; 0.5; 1.0]
    chords = [0.1; 0.1; 0.1]
    twists = 90.0 * ones(3)                 #test rotor and blade definitions.
    rotor = DuctTAPE.RotorGeometry(
        0.5, 1, rs, 0.0, chords, twists, nothing, nothing, fill(af1, 3), [-1; -1; -1], 1000
    )
    rotors = [rotor]

    #test everything together
    wakegrid = DuctTAPE.generate_wake_grid(
        ductgeometry, ductsplines, [rotor], grid_options
    )

    #get blade based on grid
    bladedims = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotor)
    blades = [bladedims]

    #generate panels
    panel_geometries = DuctTAPE.generate_panel_geometries(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    wall_panels = panel_geometries.wall_panels
    hub_panels = panel_geometries.hub_panels
    wake_panels = panel_geometries.wake_panels
    rotor_source_panels = panel_geometries.rotor_source_panels

    vinf = 10.0
    vref = vinf
    rho = 1.225
    vso = 341.0
    mu = 1.5e-8

    freestream = DuctTAPE.Freestream(vinf, vref, rho, vso, mu)

    system_aero, rotor_velocities, average_axial_velocity = initialize_system_aerodynamics(
        rotors, blades, wakegrid, rotor_source_panels, freestream; niter=10, rlx=0.5
    )

    #system_aero tests
    @test system_aero.b_gamma_grid == zeros(size(system_aero.b_gamma_grid))
    @test system_aero.delta_enthalpy_grid == zeros(size(system_aero.delta_enthalpy_grid))
    @test system_aero.delta_entropy_grid == zeros(size(system_aero.delta_entropy_grid))

    @test system_aero.b_circ_rotors == [zeros(size(system_aero.b_circ_rotors[1]))]
    @test system_aero.rotor_source_strengths ==
        zeros(size(system_aero.rotor_source_strengths))

    @test system_aero.control_point_velocities ==
        zeros(size(system_aero.control_point_velocities))

    #rotor velocity tests
    @test all(x -> x == 0.0, rotor_velocities[1].induced_axial_velocities)
    @test all(x -> x == 0.0, rotor_velocities[1].induced_radial_velocities)
    @test all(x -> x == 0.0, rotor_velocities[1].induced_tangential_velocities)

    @test all(x -> x == 10.0, rotor_velocities[1].absolute_axial_velocities)
    @test all(x -> x == 0.0, rotor_velocities[1].absolute_radial_velocities)
    @test all(x -> x == 0.0, rotor_velocities[1].absolute_tangential_velocities)

    @test all(x -> x == 10.0, rotor_velocities[1].relative_axial_velocities)
    @test all(x -> x == 0.0, rotor_velocities[1].relative_radial_velocities)
    @test rotor_velocities[1].relative_tangential_velocities ==
        -DuctTAPE.get_omega(1000.0) .* bladedims.rdim

    #average axial velocity test
    @test average_axial_velocity == 10.0 #vinf only
end
