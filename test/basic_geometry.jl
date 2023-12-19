@testset "Basic Geometry" begin

    #---------------------------------#
    #             BODIES              #
    #---------------------------------#
    # duct geometry
    duct_x = [1.0; 0.5; 0.0; 0.5; 1.0]
    duct_r = [0.0; -0.5; 0.0; 0.5; 0.0]
    duct_coordinates = [duct_x duct_r]

    # hub geometry
    hub_x = [0.0; 0.5; 1.0]
    hub_r = [0.0; 0.5; 0.0]
    hub_coordinates = [hub_x hub_r]

    # required rotor information
    rotorzlocs = [0.5; 0.75]
    Rtip = 1.5 # leading rotor tip radius

    # non-dimensional wake length
    wake_length = 1.0

    # number of wakes for each rotor
    nwake_sheets = 3

    # number of panels between discrete points
    # in this case, 2 panels between rotors, 2 panels between last rotor and hub te, and 2 panels from duct TE to end of wake
    npanels = [2; 2; 2]

    # discretize the wake
    xwake, rotor_indices = DuctAPE.discretize_wake(
        duct_coordinates, hub_coordinates, rotorzlocs, wake_length, nwake, npanels
    )

    # test that the correct number of x-stations are defined in the wake
    @test length(xwake) == sum(npanels) + 1
    # test that the correct placements of the x-stations are defined in the wake
    @test xwake == [0.5, 0.625, 0.75, 0.875, 1.0, 1.5, 2.0]
    # test that the correct wake indices are defined for the rotor locations
    @test rotor_indices == [1, 3]

    # number of panels between hub leading edge and first rotor
    nhub = 1

    # number of panels between duct leading edge and first rotor
    nduct_inner = 1

    # number of panels on duct outer surface
    nduct_outer = 2

    # update the body paneling to match the wake discretization
    new_duct_xr, new_hub_xr = DuctAPE.update_body_geometry(
        duct_coordinates,
        hub_coordinates,
        xwake,
        nhub,
        nduct_inner,
        nduct_outer;
        finterp=FLOWMath.linear,
    )

    # test that the inlet hub values are correct
    @test new_hub_xr[1:2, :] == hub_coordinates[1:2, :]
    # test that the hub coordinates use the wake discretization
    @test all(new_hub_xr[:, 1] .== [hub_x[1:2]; xwake[2:5]])
    # test that the radial hub coordinates make sense
    @test all(new_hub_xr[:, 2] .== [0.0; 0.5; 0.375; 0.25; 0.125; 0.0])

    # Move duct geometry into position and get rotor radii
    trans_duct_xr, Rtips, Rhubs = DuctAPE.place_duct(
        new_duct_xr, new_hub_xr, Rtip, rotorzlocs
    )

    # check the duct geometry
    @test all(
        new_duct_xr .≈ [
            1.0 2.0
            0.875 1.875
            0.75 1.75
            0.625 1.625
            0.5 1.5
            0.0 2.0
            0.5 2.5
            1.0 2.0
        ],
    )

    # test the rotor radii
    @test Rtips == [1.5; 1.75]
    @test Rhubs == [0.5; 0.25]

    #---------------------------------#
    #              Wake               #
    #---------------------------------#
    # get discretization of wakes at leading rotor position
    rwake = range(Rhubs[1], Rtips[1], nwake_sheets)

    # initialized wake grid
    xgrid, rgrid = DuctAPE.initialize_wake_grid(trans_duct_xr, new_hub_xr, xwake, rwake)

    #work out conservation of mass stuff by hand
    # full areas
    a1 = pi * (Rtips[1]^2 - Rhubs[1]^2)
    a2 = pi * (Rtips[2]^2 - Rhubs[2]^2)

    # section area on first rotor
    a11 = pi * (1.0^2 - 0.5^2)

    # conserve mass on second rotor
    # non-dimensionalize based on v1
    # assume constant velocities
    v1= 1.0
    v2 = a1/a2
    a21 = a11 / v2

    #solve for required radius
    r2mid = sqrt(a21 / (pi) + Rhubs[2]^2)

    #test that things line up right on the second rotor
    @test rgrid[3,2] == r2mid

    # Relax Grid
    # define simpler geometry for testing
    xs = [0.0; 0.5; 1.0]
    xgridtrivial = repeat(xs; inner=(1, 3))
    rs = [0.0; 0.5; 1.0]
    rgridtrivial = repeat(rs'; inner=(3, 1))
    DuctAPE.relax_grid!(
        xgridtrivial, rgridtrivial; max_iterations=100, tol=1e-9, verbose=false
    )

    #test that nothing changed in the update (chose geometry to do this)
    @test xgridtrivial == repeat(xs; inner=(1, 3))
    @test rgridtrivial == repeat(rs'; inner=(3, 1))

    #---------------------------------#
    #               ROTOR             #
    #---------------------------------#

    # rotor parameters
    rotor1_parameters = (;
        B=2,
        Omega=50,
        rotorzloc=rotorzlocs[1],
        rblade=[0.0, 1.0],
        chords=[0.25, 0.25],
        twists=[0.0, 0.0],
        airfoils=[nothing, nothing],
    )

    # stator parameters
    rotor2_parameters = (; rotor1_parameters..., rotorzloc=rotorzlocs[2])

    # array with rotor and stator parameters
    rotor_parameters = [rotor1_parameters, rotor2_parameters]

    # generate rotor source panel objects
    rotor_source_panels = [
        DuctAPE.generate_rotor_panels(rotorzlocs[i], rgrid[rotor_indices[i], :]) for
        i in 1:length(rotorzlocs)
    ]

    #test that panel centers are where they should be
    @test rotor_source_panels[1].panel_center == [0.5 0.75; 0.5 1.25]
    @test rotor_source_panels[2].panel_center == [
        0.75 (rgrid[rotor_indices[2], 2] + rgrid[rotor_indices[2], 1])/2
        0.75 (rgrid[rotor_indices[2], 2] + rgrid[rotor_indices[2], 3])/2
    ]
    #test that other rotor panel values are correct
    @test rotor_source_panels[1].totpanel == 2
    @test rotor_source_panels[2].totpanel == 2
    @test all(rotor_source_panels[1].panel_angle .== pi / 2)
    @test all(rotor_source_panels[2].panel_angle .== pi / 2)
    @test all(rotor_source_panels[1].panel_curvature .== 0.0)
    @test all(rotor_source_panels[2].panel_curvature .== 0.0)
    @test all(rotor_source_panels[1].panel_length .== 0.5)
    @test all(
        rotor_source_panels[2].panel_length .== [
            (rgrid[rotor_indices[2], 2] - rgrid[rotor_indices[2], 1])
            (rgrid[rotor_indices[2], 3] - rgrid[rotor_indices[2], 2])
        ],
    )
    @test all(rotor_source_panels[1].panel_normal .== [-1.0 0.0; -1.0 0.0])
    @test all(rotor_source_panels[2].panel_normal .== [-1.0 0.0; -1.0 0.0])

    # blade elements for rotors
    blade_elements = [
        DuctAPE.generate_blade_elements(
            rotor_parameters[i].B,
            rotor_parameters[i].Omega,
            rotor_parameters[i].rotorzloc,
            rotor_parameters[i].rblade,
            rotor_parameters[i].chords,
            rotor_parameters[i].twists,
            rotor_parameters[i].airfoils,
            Rtips[i],
            Rhubs[i],
            rotor_source_panels[i].panel_center[:, 2],
        ) for i in 1:length(rotorzlocs)
    ]

    #test that things got put together correctly
    @test blade_elements[1].B == 2
    @test blade_elements[1].rbe == [0.75; 1.25]
    @test blade_elements[1].chords == [0.25; 0.25]
    @test blade_elements[1].twists == [0.0; 0.0]
    @test blade_elements[1].inner_fraction == [0.25; 0.75]
    @test blade_elements[1].solidities ==
        rotor_parameters[1].chords ./
          (2 * pi * rotor_source_panels[1].panel_center[:, 2] / 2)
end
