#=

Tests with Basic Geometry to make sure things are working initially

=#

# @testset "Basic Geometry Tests" begin

    # - Put together basic Geometry - #
    # Duct
    duct_x = [1.0; 0.5; 0.0; 0.5; 1.0]
    duct_r = [1.0; 0.75; 1.0; 1.25; 1.0]

    # Hub
    hub_x = [0.0; 0.5; 1.0]
    hub_r = [0.0; 0.25; 0.0]

    # - Put togethe Basic Rotor - #
    rotor_x_position = 0.5
    nbe = 2
    nbe_fine = 3
    chords = [0.25; 0.25]
    twists = [0.0; 0.0]
    airfoils = [nothing; nothing]
    radial_positions = [0.0; 1.0]
    B = 2

    function fun(duct_x)
        duct_coordinates = [duct_x duct_r]
        hub_coordinates = [hub_x hub_r]

        body_geometry, body_panels = dt.generate_body_geometry(
            duct_coordinates, hub_coordinates
        )

        rotor_blade_elements, rotor_panels = dt.generate_blade_elements(
            rotor_x_position,
            radial_positions,
            chords,
            twists,
            airfoils,
            nbe_fine,
            B,
            body_geometry,
        )

        stator_blade_elements, stator_panels = dt.generate_blade_elements(
            rotor_x_position + 0.25,
            radial_positions,
            chords,
            twists,
            airfoils,
            nbe_fine,
            B,
            body_geometry,
        )

        x_grid_points, r_grid_points, nx, nr, rotoridxs = dt.initialize_grid_points(
            body_geometry,
            [rotor_blade_elements; stator_blade_elements];
            # [rotor_blade_elements];
            wake_length=1.0,
            debug=false,
        )

        dt.relax_grid!(x_grid_points, r_grid_points, nx, nr)

        return r_grid_points[:, 2]
    end

    findiff_j = fnd.finite_difference_jacobian(fun, duct_x)
    fordiff_j = frd.jacobian(fun, duct_x)
# end
