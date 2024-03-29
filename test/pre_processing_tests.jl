#=
Tests for all the components of precomputed_inputs function
=#
println("\nPRECOMPUTED ROTOR & WAKE INPUTS")

#simple geometry to work with:
# include("data/basic_two_rotor_for_test_NEW.jl")

@testset "Rotor/Wake Geometry Initialization" begin
    # get input data
    include("data/basic_two_rotor_for_test_NEW.jl")

    zwake, rotor_indices_in_wake = dt.discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotorstator_parameters.rotorzloc, # rotor axial locations
        paneling_constants.wake_length,
        paneling_constants.npanels,
        paneling_constants.dte_minus_cbte;
    )

    @test zwake == [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    @test rotor_indices_in_wake == [1, 3]

    # - Get Problem Dimensions - #
    problem_dimensions = dt.get_problem_dimensions(paneling_constants)

    (;
        nrotor, # number of rotors
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nws,    # number of wake sheets (also rotor nodes)
        nwsn,   # number of nodes in each wake sheet
    ) = problem_dimensions

    TF = Float64
    wake_grid = zeros(TF, 2, nwsn, nws)
    rp_duct_coordinates = zeros(TF, 2, ndn)
    rp_centerbody_coordinates = zeros(TF, 2, ncbn)
    rotor_indices_in_wake = ones(Int, nrotor)

    dt.reinterpolate_bodies!(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        duct_coordinates,
        centerbody_coordinates,
        zwake,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    @test size(rp_duct_coordinates, 1) < size(rp_duct_coordinates, 2)
    @test size(rp_centerbody_coordinates, 1) < size(rp_centerbody_coordinates, 2)
    @test size(rp_duct_coordinates, 2) == 2 * size(duct_coordinates, 1) - 1
    @test size(rp_centerbody_coordinates, 2) == 2 * size(centerbody_coordinates, 1) - 1

    @test rp_duct_coordinates == [
        1.0 0.75 0.5 0.25 0.0 0.25 0.5 0.75 1.0
        2.0 1.75 1.5 1.75 2.0 2.25 2.5 2.25 2.0
    ]
    @test rp_centerbody_coordinates == [
        0.0 0.25 0.5 0.75 1.0
        0.0 0.25 0.5 0.25 0.0
    ]

    rpb4 = copy(rp_duct_coordinates)

    dt.place_duct!(
        rp_duct_coordinates,
        rotorstator_parameters.Rtip[1],
        rotorstator_parameters.rotorzloc[1],
        rotorstator_parameters.tip_gap[1],
    )

    @test rp_duct_coordinates[1, :] == rpb4[1, :]
    @test rp_duct_coordinates[2, :] == rpb4[2, :] .- 0.75

    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotorstator_parameters.tip_gap,
        rotorstator_parameters.rotorzloc,
    )

    @test all(Rtips .== 1.0)
    @test all(Rhubs .== 0.25)

    rwake = range(Rhubs[1], Rtips[1]; length=paneling_constants.nwake_sheets)

    # Check grid initialization
    grid = dt.initialize_wake_grid(
        rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )

    @test grid[:, :, 1] == [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        0.25 0.5 0.25 0.0 0.0 0.0 0.0 0.0
    ]
    @test isapprox(
        grid[:, :, 2],
        [
            0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
            0.625 0.599479 0.625 0.73951 0.73951 0.73951 0.73951 0.73951
        ],
        atol=1e-6,
    )
    @test grid[:, :, 3] == [
        0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0
        1.0 0.75 1.0 1.25 1.25 1.25 1.25 1.25
    ]

    # check wake paneling
    wake_vortex_panels = dt.generate_wake_panels(grid[:, :, 1:length(rwake)])

    # # TODO: add any checks for things that get used.
    # @test wake_vortex_panels.node[1,:] == reduce(vcat,grid[1,:,:])
    # @test wake_vortex_panels.node[2,:] == reduce(vcat,grid[2,:,:])
    # @test wake_vortex_panels.endpanelidxs == [1 8 15; 7 14 21]
    # @test wake_vortex_panels.endnodeidxs == [1 9 17; 8 16 24]

    # check wakeK calcualtion
    wakeK = dt.get_wake_k(wake_vortex_panels.node[2, :], problem_dimensions.nwn)

    @test isapprox(
        wakeK,
        [
            -0.20264236728467555
            -0.05066059182116889
            -0.20264236728467555
            0.0
            0.0
            0.0
            0.0
            0.0
            -0.03242277876554809
            -0.03524215083211749
            -0.03242277876554809
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.0231591276896772
            -0.012665147955292222
            -0.02251581858718617
            -0.012665147955292222
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
            -0.008105694691387022
        ],
    )

    # just make sure it converges
    dt.relax_grid!(dt.SLORWakeSolverOptions(), grid; verbose=false, silence_warnings=true)

    # Check grid initialization
    # re-set up initial grid for easier testing
    grid = dt.initialize_wake_grid(
        rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )

    num_rotors = length(rotorstator_parameters)

    # rotor source panel objects
    # TODO: incorrect function

    rotor_source_panels = dt.generate_rotor_panels(
        rotorzloc, grid, [1, 3], paneling_constants.nwake_sheets
    )

    # rotor blade element objects
    blade_elements = dt.interpolate_blade_elements(
        rotorstator_parameters,
        Rtips,
        Rhubs,
        rotor_source_panels.controlpoint[2, :],
        problem_dimensions.nbe,
    )

    @test blade_elements.inner_fraction == [0.75 0.75; 0.25 0.25]
    @test blade_elements.stagger == [70 70; 70 70] * pi / 180
    @test blade_elements.rotor_panel_centers == [0.4375 0.4375; 0.8125 0.8125]
    @test blade_elements.B == [2; 4]
    @test blade_elements.chords == [0.1 0.1; 0.1 0.1]
    @test blade_elements.twists == [20 20; 20 20] * pi / 180
    @test isapprox(
        blade_elements.solidity,
        [
            0.07275654541343787 0.14551309082687575
            0.039176601376466544 0.07835320275293309
        ],
        atol=1e-6,
    )
    @test blade_elements.fliplift == [false, false]
    @test blade_elements.Rhub == rotorstator_parameters.Rhub
    @test blade_elements.Rtip == rotorstator_parameters.Rtip
end

# TODO: put together a hole new test here for each initialization function
# Make sure that the initial conditions make sense
@testset "Rotor/Wake Aero Initialization" begin
    include("data/basic_two_rotor_for_test_NEW.jl")
    options = dt.set_options()

    # Allocate Cache
    solve_parameter_caching = dt.allocate_solve_parameter_cache(
        options.solve_options, propulsor.paneling_constants
    )

    # separate out caching items
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, Float64(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = dt.withdraw_solve_parameter_cache(
        options.solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(propulsor.operating_point))
        solve_parameter_tuple.operating_point[f] .= getfield(propulsor.operating_point, f)
    end

    # - Do precomputations - #
    ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps, panels, problem_dimensions = dt.precompute_parameters_iad!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        propulsor;
        wake_solve_options=options.wake_options,
        autoshiftduct=options.autoshiftduct,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    (; vz_rotor, vtheta_rotor, Cm_wake, operating_point, linsys, ivr, ivw) =
        solve_parameter_tuple

    # initialize velocities
    dt.initialize_velocities!(
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        solve_parameter_tuple.operating_point,
        (; solve_parameter_tuple.blade_elements..., airfoils...),
        (; solve_parameter_tuple.linsys..., A_bb_LU),
        ivr,
        ivw,
        idmaps.body_totnodes,
        idmaps.wake_panel_sheet_be_map,
    )

    # check that the wake velocities are taken straight back from one rotor to the next, then to the end of the wake
    cmcheck = reshape(Cm_wake, 7, 3)
    @test all(cmcheck[1, :] .== cmcheck[2, :])
    for i in 3:6
        @test all(cmcheck[i, :] .== cmcheck[i + 1, :])
    end
    @test all(cmcheck .!= 0.0)

    # TODO: need to add a better test with more realistic geometry that you can draw more conclusions from

    # CSOR Solve initialization
    options = dt.set_options(; solve_options=dt.CSORSolverOptions(), wake_options=dt.SLORWakeSolverOptions())

    # Allocate Cache
    solve_parameter_caching = dt.allocate_solve_parameter_cache(
        options.solve_options, propulsor.paneling_constants
    )

    # separate out caching items
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, Float64(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = dt.withdraw_solve_parameter_cache(
        options.solve_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(propulsor.operating_point))
        solve_parameter_tuple.operating_point[f] .= getfield(propulsor.operating_point, f)
    end

    # - Do precomputations - #
    ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps, panels, problem_dimensions = dt.precompute_parameters_iad!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        propulsor;
        wake_solve_options=options.wake_options,
        autoshiftduct=options.autoshiftduct,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    (; Gamr, sigr, gamw, operating_point, linsys, ivr, ivw, wakeK) = solve_parameter_tuple

    # initialize velocities
    dt.initialize_strengths!(
        Gamr,
        sigr,
        gamw,
        operating_point,
        (; blade_elements..., airfoils...),
        (; linsys..., A_bb_LU),
        ivr,
        ivw,
        wakeK,
        idmaps.body_totnodes,
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
        idmaps.wake_panel_sheet_be_map,
        idmaps.wake_node_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface,
    )

    #centerbody trailing edge is on axis, so all gamw along hub should be zero
    @test all(iszero.(gamw[1:8]))

    # check that the end of the wake behaves as expected
    gamwcheck = reshape(gamw, 8, 3)
    for i in 4:7
        @test isapprox(gamwcheck[i, 2], gamwcheck[i + 1, 2], atol=1e-1)
        @test gamwcheck[i, 3] .== gamwcheck[i + 1, 3]
    end

    # TODO: need to add a better test with more realistic geometry that you can draw more conclusions from

end

