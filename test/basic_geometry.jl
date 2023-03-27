@testset "Body Geometry" begin

    # duct geometry
    duct_x = [1.25; 0.375; 0.0; 0.375; 1.25]
    duct_r = [1.0; 0.875; 1.0; 1.375; 1.0]
    duct_xr = [duct_x duct_r]

    # hub geometry
    hub_x = [0.0; 0.5; 1.0]
    hub_r = [0.0; 0.25; 0.0]
    hub_xr = [hub_x hub_r]

    # rotor parameters
    rotor1_parameters = (;
        B=2,
        omega=50,
        xrotor=0.5,
        rblade=[0.0, 1.0],
        chords=[0.25, 0.25],
        twists=[0.0, 0.0],
        solidity=[0.0],
        airfoils=[nothing, nothing],
    )

    # stator parameters
    rotor2_parameters = (; rotor1_parameters..., xrotor=0.75)

    # array with rotor and stator parameters
    rotor_parameters = [rotor1_parameters, rotor2_parameters]

    # non-dimensional wake length
    wake_length = 1.0

    # number of wakes for each rotor
    nwake = 10

    # discretize the wake
    xwake = DuctTAPE.discretize_wake(duct_xr, hub_xr, rotor_parameters, wake_length, nwake)

    # number of panels between hub leading edge and first rotor
    nhub = 2

    # number of panels between duct leading edge and first rotor
    nduct_inner = 2

    # number of panels on duct outer surface
    nduct_outer = 4

    # update the body paneling to match the wake discretization
    new_duct_xr, new_hub_xr = DuctTAPE.update_body_geometry(
        duct_xr, hub_xr, xwake, nhub, nduct_inner, nduct_outer; finterp=FLOWMath.linear
    )

    # check the duct geometry
    @test all(
        new_duct_xr .≈ [
            1.0 1.0
            0.95 0.975
            0.9 0.95
            0.85 0.925
            0.8 0.9
            0.75 0.875
            0.7 0.85
            0.65 0.825
            0.6 0.8
            0.55 0.775
            0.5 0.75
            0.25 0.875
            0.0 1.0
            0.25 1.125
            0.5 1.25
            0.75 1.125
            1.0 1.0
        ],
    )

    # check the hub geometry
    @test all(
        new_hub_xr .≈ [
            0.0 0.0
            0.25 0.125
            0.5 0.25
            0.55 0.225
            0.6 0.2
            0.65 0.175
            0.7 0.15
            0.75 0.125
            0.8 0.1
            0.85 0.075
            0.9 0.05
            0.95 0.025
            1.0 0.0
        ],
    )

    # uncomment to plot
    using Plots
    pyplot()
    plot(duct_x, duct_r; label="duct input", color=:blue, linestyle=:dot)
    plot!(new_duct_xr[:, 1], new_duct_xr[:, 2]; color=:black, marker=true, label="")
    plot!(new_hub_xr[:, 1], new_hub_xr[:, 2]; color=:black, marker=true, label="")
    plot!([0.5, 0.5], [0.25, 0.75]; color=:black, linestyle=:dash, label="")
    plot!([0.75, 0.75], [0.125, 0.875]; color=:black, linestyle=:dash, label="")
    plot!(; aspect_ratio=1.0)

    # discretize first rotor
    xrotor = rotor_parameters[1].xrotor
    _, leidx = findmin(view(new_duct_xr, :, 1))
    _, ihub = findmin(x -> abs(x - xrotor), view(new_hub_xr, :, 1))
    _, iduct = findmin(x -> abs(x - xrotor), view(new_duct_xr, 1:leidx, 1))
    Rhub = new_hub_xr[ihub, 2]
    Rtip = new_duct_xr[iduct, 2]
    rwake = range(Rhub, Rtip, nwake)

    # generate blade elements for rotor
    first_blade_elements = DuctTAPE.generate_blade_elements(
        rotor_parameters[1].B,
        rotor_parameters[1].omega,
        rotor_parameters[1].xrotor,
        rotor_parameters[1].rblade,
        rotor_parameters[1].chords,
        rotor_parameters[1].twists,
        rotor_parameters[1].solidity,
        rotor_parameters[1].airfoils,
        new_duct_xr,
        new_hub_xr,
        xwake,
        rwake,
    )

    # initialize wake grid
    xgrid, rgrid = DuctTAPE.initialize_wake_grid(new_duct_xr, new_hub_xr, xwake, rwake)

    # uncomment to plot
    using Plots
    pyplot()
    plot(new_duct_xr[:, 1], new_duct_xr[:, 2]; color=:black, marker=true, label="")
    plot!(new_hub_xr[:, 1], new_hub_xr[:, 2]; color=:black, marker=true, label="")
    plot!([0.5, 0.5], [0.25, 0.75]; color=:black, linestyle=:dash, label="")
    plot!([0.75, 0.75], [0.125, 0.875]; color=:black, linestyle=:dash, label="")
    plot!(xgrid, rgrid; color=:blue, label="")
    plot!(; aspect_ratio=1.0)
end

# rotor geometry
xrotor = 0.5

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
end

findiff_j = fnd.finite_difference_jacobian(fun, duct_x)
fordiff_j = frd.jacobian(fun, duct_x)
# end
end

findiff_j = fnd.finite_difference_jacobian(fun, duct_x)
fordiff_j = frd.jacobian(fun, duct_x)
# end
