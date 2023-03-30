@testset "Body Geometry" begin

    # duct geometry
    duct_x = [1.25; 0.375; 0.0; 0.375; 1.25]
    duct_r = [1.0; 0.875; 1.0; 1.375; 1.0]
    duct_coordinates = [duct_x duct_r]

    # hub geometry
    hub_x = [0.0; 0.35; 1.0]
    hub_r = [0.0; 0.25; 0.0]
    hub_coordinates = [hub_x hub_r]

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

    xrotors = [0.5; 0.75]
    Rtip = 1.0 # leading rotor tip radius

    # non-dimensional wake length
    wake_length = 1.0

    # number of wakes for each rotor
    nwake = 10

    # number of panels between discrete points
    # in this case, 5 panels between rotors, 5 panels between last rotor and hub te, 3 panels between hub and duct te's, and 20 panels from duct TE to end of wake
    npanels = [5; 5; 3; 20]

    # discretize the wake
    xwake, rotor_indices = DuctTAPE.discretize_wake(
        duct_coordinates, hub_coordinates, xrotors, wake_length, nwake, npanels
    )

    # number of panels between hub leading edge and first rotor
    nhub = 5

    # number of panels between duct leading edge and first rotor
    nduct_inner = 5

    # number of panels on duct outer surface
    nduct_outer = 10

    # update the body paneling to match the wake discretization
    new_duct_xr, new_hub_xr = DuctTAPE.update_body_geometry(
        duct_coordinates,
        hub_coordinates,
        xwake,
        nhub,
        nduct_inner,
        nduct_outer;
        finterp=FLOWMath.linear,
    )

    trans_duct_xr, Rtips, Rhubs = DuctTAPE.place_duct(
        new_duct_xr, new_hub_xr, Rtip, xrotors
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
    # pyplot()
    plot(
        duct_x,
        duct_r;
        label="input geometry",
        color=mycolors[1],
        marker=true,
        linestyle=:dot,
    )
    plot!(hub_x, hub_r; label="", color=mycolors[1], marker=true, linestyle=:dot)
    plot!(
        xwake,
        0.75 * ones(length(xwake));
        seriestype=:scatter,
        markersize=2,
        label="wake x-coordinates",
    )
    plot!(
        trans_duct_xr[:, 1],
        trans_duct_xr[:, 2];
        color=mycolors[2],
        markersize=2,
        marker=true,
        label="re-paneled, shifted bodies",
    )
    plot!(
        new_hub_xr[:, 1],
        new_hub_xr[:, 2];
        color=mycolors[2],
        marker=true,
        markersize=2,
        label="",
    )
    plot!(
        xrotors[1] * ones(2),
        [Rhubs[1]; Rtips[1]];
        color=:black,
        linestyle=:dash,
        label="rotor locations",
    )
    plot!(
        xrotors[2] * ones(2), [Rhubs[2]; Rtips[2]]; color=:black, linestyle=:dash, label=""
    )

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
    xgrid, rgrid = DuctTAPE.initialize_wake_grid(trans_duct_xr, new_hub_xr, xwake, rwake)

    # uncomment to plot
    plot!(xgrid, rgrid; markershape=:vline, color=mycolors[4], label="")

    # Relax Grid
    # DuctTAPE.relax_grid!(xgrid, rgrid; max_iterations=100, tol=1e-9, verbose=false)
    # plot!(xgrid, rgrid; markershape=:vline, color=mycolors[5], label="")

    #
end
