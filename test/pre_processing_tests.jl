#=
Tests for all the components of precomputed_inputs function
=#
println("\nPRECOMPUTED ROTOR & WAKE INPUTS")

#simple geometry to work with:
# include("data/basic_two_rotor_for_test_NEW.jl")

@testset "Rotor/Wake Geometry Initialization" begin
    # get input data
    include("data/basic_two_rotor_for_test_NEW.jl")

    zwake, rotor_indices_in_wake, duct_le_coordinates = dt.discretize_wake(
        duct_coordinates,
        centerbody_coordinates,
        rotor.rotorzloc, # rotor axial locations
        paneling_constants.wake_length,
        paneling_constants.npanels,
        paneling_constants.dte_minus_cbte;
        le_bracket=1,
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
        duct_le_coordinates,
        paneling_constants.ncenterbody_inlet,
        paneling_constants.nduct_inlet;
        finterp=fm.linear,
    )

    @test size(rp_duct_coordinates, 1) < size(rp_duct_coordinates, 2)
    @test size(rp_centerbody_coordinates, 1) < size(rp_centerbody_coordinates, 2)
    @test size(rp_duct_coordinates, 2) == 2 * size(duct_coordinates, 1) - 1
    @test size(rp_centerbody_coordinates, 2) == 2 * size(centerbody_coordinates, 1) - 1

    @test isapprox(
        rp_duct_coordinates,
        [
            1.0 0.75 0.5 0.25 0.0 0.25 0.5 0.75 1.0
            2.0 1.75 1.5 1.75 2.0 2.25 2.5 2.25 2.0
        ],
    )
    @test rp_centerbody_coordinates == [
        0.0 0.25 0.5 0.75 1.0
        0.0 0.25 0.5 0.25 0.0
    ]

    rpb4 = copy(rp_duct_coordinates)

    dt.place_duct!(rp_duct_coordinates, rotor.Rtip[1], rotor.rotorzloc[1], rotor.tip_gap[1])

    @test rp_duct_coordinates[1, :] == rpb4[1, :]
    @test rp_duct_coordinates[2, :] == rpb4[2, :] .- 0.75

    Rtips, Rhubs = dt.get_blade_ends_from_body_geometry(
        rp_duct_coordinates, rp_centerbody_coordinates, rotor.tip_gap, rotor.rotorzloc
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
    wakeK = dt.get_wake_k(wake_vortex_panels.node[2, :])

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
    dt.relax_grid!(dt.SLORGridSolverOptions(), grid; verbose=false, silence_warnings=true)

    # Check grid initialization
    # re-set up initial grid for easier testing
    grid = dt.initialize_wake_grid(
        rp_duct_coordinates, rp_centerbody_coordinates, zwake, rwake
    )

    num_rotors = length(rotor.B)

    # rotor source panel objects
    rotor_source_panels = dt.generate_rotor_panels(
        rotor.rotorzloc, grid, [1, 3], paneling_constants.nwake_sheets
    )

    # rotor blade element objects
    blade_elements, airfoils = dt.interpolate_blade_elements(
        rotor, Rtips, Rhubs, rotor_source_panels.controlpoint[2, :], problem_dimensions.nbe
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
    @test blade_elements.Rhub == rotor.Rhub
    @test blade_elements.Rtip == rotor.Rtip
end

@testset "Bookkeeping Tests" begin

    # hub is longer than duct
    npanels = [2, 1, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = -1 # paneling_constants.dte_minus_cbte
    wake_nodemap = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    ndp = 8 # pd.ndp
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        ndp,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3, 4]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4, 5]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [2]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [11]

    # hub is shorter than duct
    npanels = [2, 1, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = 1 # paneling_constants.dte_minus_cbte
    wake_nodemape = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    ndp = 8 # pd.ndp
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        ndp,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13, 14]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10, 11]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [3, 2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [3]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [11]

    # hub and duct have same TE axial coordinate
    npanels = [3, 1] #npanels
    ncenterbody_inlet = 2 # paneling_constants.ncenterbody_inlet
    nwake_sheets = 3 # paneling_constants.nwake_sheets
    dte_minus_cbte = 0 # paneling_constants.dte_minus_cbte
    wake_nodemape = [
        1 2 3 4 6 7 8 9 11 12 13 14
        2 3 4 5 7 8 9 10 12 13 14 15
    ]# wake_vortex_panels.nodemap
    wake_endnodeidxs = [1 6 11; 5 10 15] # wake_vortex_panels.endnodeidxs
    nwp = 12 # pd.nwp
    nwsp = 4 # pd.nwsp
    nbn = 12 # pd.nbn
    ndp = 8 # pd.ndp
    rotor_indices_in_wake = [1] # rotor_indices_in_wake
    nrotor = 1 # pd.nrotor

    idmaps = dt.set_index_maps(
        npanels,
        ncenterbody_inlet,
        nwake_sheets,
        dte_minus_cbte,
        wake_nodemap,
        wake_endnodeidxs,
        nwp,
        nwsp,
        nbn,
        ndp,
        rotor_indices_in_wake,
        nrotor,
    )

    @test idmaps.wake_node_sheet_be_map[:, 1] ==
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
    @test all(idmaps.wake_node_sheet_be_map[:, 2] .== 1)
    @test idmaps.wake_panel_sheet_be_map[:, 1] == [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    @test all(idmaps.wake_panel_sheet_be_map[:, 2] .== 1)

    @test idmaps.wake_node_ids_along_casing_wake_interface == [11, 12, 13, 14]
    @test idmaps.wake_node_ids_along_centerbody_wake_interface == [1, 2, 3, 4]

    @test idmaps.wake_panel_ids_along_casing_wake_interface == [9, 10, 11]
    @test idmaps.wake_panel_ids_along_centerbody_wake_interface == [1, 2, 3]

    @test idmaps.centerbody_panel_ids_along_centerbody_wake_interface == [3, 4, 5]
    @test idmaps.duct_panel_ids_along_centerbody_wake_interface == [3, 2, 1]

    @test idmaps.id_of_first_casing_panel_aft_of_each_rotor == [3]
    @test idmaps.id_of_first_centerbody_panel_aft_of_each_rotor == [11]
end

@testset "Rotor/Wake Aero Initialization" begin
    r1 = [0.25; 0.5; 0.75; 1.0]
    Rtip = [1.0, 1.0]
    rnondim1 = r1 ./ Rtip[1]
    rnondim = [rnondim1 rnondim1]
    afparams1 = dt.c4b.DFDCairfoil()
    rotorzloc = [0.25, 0.75]
    r = rnondim
    chords = 0.1 * ones(size(rnondim))
    twists = 20.0 * pi / 180.0 * ones(size(rnondim))
    airfoils = fill(afparams1, 4, 2)
    Rhub = [0.25, 0.25]
    Rtip = Rtip
    tip_gap = [0.0, 0.0]
    B = [2, 4]
    fliplift = [0.0, 0.0]

    rotor = dt.Rotor(
        B, rotorzloc, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, fliplift
    )

    ncenterbody_inlet = 1
    nduct_inlet = 1
    nwake_sheets = 3
    wake_length = 1.0
    npanels = [2, 1, 4]
    dte_minus_cbte = 0

    paneling_constants = dt.PanelingConstants(
        nduct_inlet, ncenterbody_inlet, npanels, dte_minus_cbte, nwake_sheets, wake_length
    )

    Vinf = [10.0]
    rhoinf = [1.226]
    muinf = [1.78e-5]
    asound = [340.0]
    Omega = [5000.0, 0.0] * pi / 30  # convert from RPM to rad/s

    operating_point = dt.OperatingPoint(Vinf, Omega, rhoinf, muinf, asound)

    Vref = [10.0]
    Rref = [Rtip]
    reference_parameters = dt.ReferenceParameters(Vref, Rref)

    duct_coordinates = [1.0 2.0; 0.5 1.5; 0.0 2.0; 0.5 2.5; 1.0 2.0]
    centerbody_coordinates = [0.0 0.0; 0.5 0.5; 1.0 0.0]

    ducted_rotor = dt.DuctedRotor(
        duct_coordinates, centerbody_coordinates, rotor, paneling_constants
    )

    options = dt.set_options()

    # - Get Problem Dimensions - #
    problem_dimensions = dt.get_problem_dimensions(ducted_rotor.paneling_constants)

    # - Set up Pre- and Post-process Cache - #
    # Allocate Cache
    prepost_container_caching = dt.allocate_prepost_container_cache(
        ducted_rotor.paneling_constants
    )

    # unpack the caching
    (; prepost_container_cache, prepost_container_cache_dims) = prepost_container_caching

    # Get correct cached types
    prepost_container_cache_vec = @views PreallocationTools.get_tmp(
        prepost_container_cache, Float64(1.0)
    )

    # reset cache
    prepost_container_cache_vec .= 0

    # Reshape Cache
    prepost_containers = dt.withdraw_prepost_container_cache(
        prepost_container_cache_vec, prepost_container_cache_dims
    )

    # - Set up Solver Sensitivity Paramter Cache - #

    # Allocate Cache
    solve_parameter_caching = dt.allocate_solve_parameter_cache(
        options.solver_options, ducted_rotor.paneling_constants
    )

    # unpack caching
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, Float64(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = dt.withdraw_solve_parameter_cache(
        options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(operating_point))
        if f != :units
            solve_parameter_tuple.operating_point[f] .= getfield(operating_point, f)
        end
    end

    # - Do preprocessutations - #
    if options.verbose
        println("Pre-computing Parameters")
    end

    ##### ----- PERFORM PREPROCESSING COMPUTATIONS ----- #####

    # - Preprocess - #
    A_bb_LU, lu_decomp_flag, airfoils, idmaps, _ = dt.precompute_parameters!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        ducted_rotor,
        operating_point,
        prepost_containers,
        problem_dimensions;
        grid_solver_options=options.grid_solver_options,
        integration_options=options.integration_options,
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
        options.solver_options,
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
    options = dt.DFDC_options()

    # - Get Problem Dimensions - #
    problem_dimensions = dt.get_problem_dimensions(ducted_rotor.paneling_constants)

    # - Set up Pre- and Post-process Cache - #
    # Allocate Cache
    prepost_container_caching = dt.allocate_prepost_container_cache(
        ducted_rotor.paneling_constants
    )

    # unpack the caching
    (; prepost_container_cache, prepost_container_cache_dims) = prepost_container_caching

    # Get correct cached types
    prepost_container_cache_vec = @views PreallocationTools.get_tmp(
        prepost_container_cache, Float64(1.0)
    )

    # reset cache
    prepost_container_cache_vec .= 0

    # Reshape Cache
    prepost_containers = dt.withdraw_prepost_container_cache(
        prepost_container_cache_vec, prepost_container_cache_dims
    )

    # - Set up Solver Sensitivity Paramter Cache - #

    # Allocate Cache
    solve_parameter_caching = dt.allocate_solve_parameter_cache(
        options.solver_options, ducted_rotor.paneling_constants
    )

    # unpack caching
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, Float64(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = dt.withdraw_solve_parameter_cache(
        options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(operating_point))
        if f != :units
            solve_parameter_tuple.operating_point[f] .= getfield(operating_point, f)
        end
    end

    # - Do preprocessutations - #
    if options.verbose
        println("Pre-computing Parameters")
    end

    ##### ----- PERFORM PREPROCESSING COMPUTATIONS ----- #####

    # - Preprocess - #
    A_bb_LU, lu_decomp_flag, airfoils, idmaps, _ = dt.precompute_parameters!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        ducted_rotor,
        operating_point,
        prepost_containers,
        problem_dimensions;
        grid_solver_options=options.grid_solver_options,
        integration_options=options.integration_options,
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
        options.solver_options,
        Gamr,
        sigr,
        gamw,
        operating_point,
        (; solve_parameter_tuple.blade_elements..., airfoils...),
        (; solve_parameter_tuple.linsys..., A_bb_LU),
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

@testset "Wake Geometry Solver" begin
    function f(x)
        TF = eltype(x)
        trcz = x[1] #top right corner
        trcr = x[2] #top right corner

        # converging geometry
        wake_grid = zeros(TF, 2, 3, 3)
        wake_grid[1, :, 1] = [0.0; 0.75; 1.75] * trcz
        wake_grid[1, :, 2] = [0.0; 0.75; 1.75] * trcz
        wake_grid[1, :, 3] = [0.0; 0.75; 1.75] * trcz

        wake_grid[2, 1, :] = [0.3; 1.0; 1.7] * trcr
        wake_grid[2, 2, :] = [0.4; 1.0; 1.6] * trcr
        wake_grid[2, 3, :] = [0.5; 1.0; 1.5] * trcr

        zs = wake_grid[1, :, :]
        rs = wake_grid[2, :, :]

        # - Use NLsolve to obtain grid solution - #
        dt.solve_elliptic_grid!(
            wake_grid;
            algorithm=:trust_region,
            atol=1e-15,
            iteration_limit=100,
            converged=[false],
            verbose=false,
        )

        return reshape(wake_grid, :)
    end

    x0 = [1.0; 1.0]

    wg = f(x0)

    steps = 1e-6
    central_diff_jac = FiniteDiff.finite_difference_jacobian(
        f, x0, Val(:central); relstep=steps
    )

    forwardAD_jac = ForwardDiff.jacobian(f, x0)

    maxj, _ = findmax(
        abs.((forwardAD_jac .- central_diff_jac) ./ (1.0 .+ central_diff_jac))
    )

    @test isapprox(
        wg,
        [
            0.0
            0.3
            0.75
            0.4
            1.75
            0.5
            0.0
            1.0
            0.75
            0.9587508251715947
            1.75
            0.9387512252547923
            0.0
            1.7
            0.75
            1.6
            1.75
            1.5
        ],
        atol=1e-12,
    )

    @test maxj < 1e-9
end
