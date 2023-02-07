@testset "Geometry Initialization" begin
    @testset "Rotor Geometry" begin

        # Initialize Functions

        num_rotors = 2

        body_geometry = (
            duct_inner_spline=fm.Akima([0.0; 2.0], [0.75; 0.75]),
            hub_spline=fm.Akima([0.0; 2.0], [0.25; 0.25]),
            duct_range=[0.0; 2.0],
            hub_range=[0.0; 2.0],
        )

        rotor_parameters = (
            rotor_x_position=1.0,
            radial_positions=[0.0; 0.25; 0.5; 0.75; 1.0],
            chords=ones(5),
            twists=ones(5),
            airfoils=ones(5),
            num_radial_stations=6,
            num_blades=2,
            omega=1.0,
        )

        blade_elements = dt.initialize_blade_elements(
            rotor_parameters, body_geometry, num_rotors
        )

        @test length(blade_elements) == 2
        @test blade_elements[1].rotor_x_position == blade_elements[2].rotor_x_position
        @test all(blade_elements[1].twists .== pi / 180.0 * ones(6))
        @test blade_elements[1].radial_positions ==
            0.25 .+ 0.5 * [0.0; 0.2; 0.4; 0.6; 0.8; 1.0]
        @test blade_elements[1].omega == [1.0]
        @test blade_elements[1].num_blades == [2]
        @test blade_elements[1].num_radial_stations == [6]

        dp = dt.initialize_dummy_rotor_panels(blade_elements[1], num_rotors)

        @test all(dp[1].panel_center[:, 1] .== 1.0)
        @test all(dp[1].panel_center[:, 2] .== 0.25 .+ 0.5 * [0.0; 0.2; 0.4; 0.6; 0.8; 1.0])

        # Update Functions

        rotor_parameters = (
            rotor_x_position=0.5,
            radial_positions=[0.0; 0.25; 0.5; 0.75; 1.0],
            chords=2.0 * ones(5),
            twists=2.0 * ones(5),
            airfoils=2.0 * ones(5),
            num_radial_stations=7,
            num_blades=4,
            omega=2.0,
        )

        dt.generate_blade_elements!(
            blade_elements[2],
            rotor_parameters.rotor_x_position,
            rotor_parameters.radial_positions,
            rotor_parameters.chords,
            rotor_parameters.twists,
            rotor_parameters.airfoils,
            rotor_parameters.num_radial_stations,
            rotor_parameters.num_blades,
            rotor_parameters.omega,
            body_geometry;
            updated_radial_positions=blade_elements[1].radial_positions,
        )

        @test blade_elements[1].rotor_x_position != blade_elements[2].rotor_x_position
        @test blade_elements[2].rotor_x_position == [0.5]
        @test all(blade_elements[2].twists .== 2.0 * pi / 180.0 * ones(6))
        @test blade_elements[2].radial_positions ==
            0.25 .+ 0.5 * [0.0; 0.2; 0.4; 0.6; 0.8; 1.0]
        @test blade_elements[2].omega == [2.0]
        @test blade_elements[2].num_blades == [4]
        @test blade_elements[2].num_radial_stations == [6]

        dt.update_dummy_rotor_panels!(
            blade_elements[2], dp[2], blade_elements[2].radial_positions
        )

        @test all(dp[2].panel_center[:, 1] .== 0.5)
        @test all(dp[1].panel_center[:, 2] .== 0.25 .+ 0.5 * [0.0; 0.2; 0.4; 0.6; 0.8; 1.0])
    end
end
