println("\nITERATION STEP THROUGH TESTS")

@testset "Iteration Steps" begin
    ##### ----- Set up ----- #####

    datapath = "data/dfdc_init_iter1/"

    # - read in index file - #
    # in the case here, element 1 is the hub, 2 is the duct, 3 is the rotor, 4 is the hub wake, 14 is the duct wake, and the rest are the interior wake sheets.
    include(datapath * "element_indices.jl")

    # get ranges of element idices
    nidx = dfdc_eid[:, 2:3]
    pidx = dfdc_eid[:, 4:5]
    nidranges = nidx[:, 2] .- nidx[:, 1] .+ 1
    pidranges = pidx[:, 2] .- pidx[:, 1] .+ 1

    # - read in geometry file to help find indices - #
    include(datapath * "dfdc_geometry.jl")
    include(datapath * "ductape_formatted_dfdc_geometry.jl")
    include(datapath * "ductape_parameters.jl")
    #NOTE, for some reason, julia doesn't recognize that this was defined in the above file...
    operating_point = dt.OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)
    propulsor = dt.Propulsor(
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotorstator_parameters,
        operating_point,
        paneling_constants,
        reference_parameters,
    )

    options = dt.DFDC_options(;
        integration_options=dt.IntegrationOptions(;
            nominal=dt.GaussLegendre(30), singular=dt.GaussLegendre(30)
        ),
    )

    # rotor is at first wake index
    rotor_indices_in_wake = [1]

    problem_dimensions = dt.get_problem_dimensions(
        dt.generate_all_panels(
            rp_duct_coordinates,
            rp_centerbody_coordinates,
            paneling_constants.nwake_sheets,
            rotor_indices_in_wake,
            rotorstator_parameters.rotorzloc,
            wake_grid;
            itcpshift=options.itcpshift,
            axistol=options.axistol,
            tegaptol=options.tegaptol,
            silence_warnings=options.silence_warnings,
        )...,
    )

    # - Set up Pre- and Post-process Cache - #
    # Allocate Cache
    prepost_container_caching = dt.allocate_prepost_container_cache(problem_dimensions)

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
        options.solver_options, problem_dimensions
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
    for f in fieldnames(typeof(propulsor.operating_point))
        solve_parameter_tuple.operating_point[f] .= getfield(propulsor.operating_point, f)
    end

    # generate inputs
    ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps, problem_dimensions = dt.precompute_parameters!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        wake_grid,
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        rotor_indices_in_wake,
        rotorstator_parameters,
        paneling_constants,
        operating_point,
        prepost_containers,
        problem_dimensions;
        integration_options=options.integration_options,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    (; blade_elements, linsys, operating_point, ivr, ivw, wakeK) = solve_parameter_tuple

    solve_container_caching = dt.allocate_solve_container_cache(
        options.solver_options, problem_dimensions
    )

    (; solve_container_cache, solve_container_cache_dims) = solve_container_caching

    solve_container_cache_vector = @views PreallocationTools.get_tmp(
        solve_container_cache, Float64(1.0)
    )
    # reset cache
    solve_container_cache_vector .= 0
    solve_containers = dt.withdraw_solve_container_cache(
        options.solver_options, solve_container_cache_vector, solve_container_cache_dims
    )

    # - Parameters - #
    B = [5] #there are 5 blades

    # get wake sheet length
    num_wake_z_panels = num_wake_z_nodes - 1

    # get number of nodes and panels on bodies interface with wake
    # NOTE: these happen to be wrong, but I think it works out later...
    nhub_interface_nodes = num_wake_z_nodes - nidranges[4]
    nduct_interface_nodes = num_wake_z_nodes - nidranges[end]
    nhub_interface_panels = num_wake_z_panels - pidranges[4]
    nduct_interface_panels = num_wake_z_panels - pidranges[end]

    ##### ----- read in functions to reformat ----- #####
    include(datapath * "reformatting_functions.jl")

    # we need iter 2 sigr since that get's updated at the beginning of iteration 2
    include(datapath * "iter2_sigr.jl")
    sigr2mat = [sigr2;;]
    # we need the relaxed Gamr and gamw from iter 1, since that's what's used for the body at the beginning of iter 2
    include(datapath * "iter1_relaxed_bgam.jl")
    include(datapath * "iter1_relaxed_gamw.jl")
    Gamr1 = reformat_circulation(bgamr1, B)
    Gamr1mat = [Gamr1;;]
    gamw1 = reformat_gamw(relaxed_GTH, nidx, nhub_interface_nodes, nduct_interface_nodes)

    # solve body with those inputs and test body strengths
    ##### ----- Test body strengths ----- #####
    dt.calculate_body_vortex_strengths!(
        solve_containers.gamb,
        A_bb_LU,
        linsys.b_bf,
        -gamw1,
        linsys.A_bw,
        linsys.A_pw,
        -sigr2mat,
        linsys.A_br,
        linsys.A_pr,
        linsys.A_bb,
        solve_containers.rhs,
    )

    # - test body strengths against DFDC written values - #
    # get the body strengths from DFDC
    include(datapath * "iter2_gamb.jl")
    gamb2 = reformat_sol(res2, nidx)
    @test maximum(abs.(solve_containers.gamb .+ gamb2)) < 1.5

    ##### ----- Test blade element values----- #####
    include(datapath * "iter2_blade_element_values.jl")
    b2 = reformat_rotor_velocities(rotor_vels2)
    include(datapath * "iter2_extended_blade_element_values.jl")
    be2 = reformat_blade_elements(extended_blade_element_values2)

    dt.calculate_induced_velocities_on_rotors!(
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        Gamr1mat,
        -gamw1,
        -sigr2mat,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        @view(ivr.v_rw[:, :, 1]),
        @view(ivr.v_rr[:, :, 1]),
        @view(ivr.v_rb[:, :, 1]),
        blade_elements.B,
        blade_elements.rotor_panel_centers,
    )
    dt.reframe_rotor_velocities!(
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        operating_point.Vinf[],
        operating_point.Omega,
        blade_elements.rotor_panel_centers,
    )

    @test isapprox(solve_containers.Cz_rotor, be2.Wz, atol=5e-1)
    @test isapprox(solve_containers.Cmag_rotor, be2.Wmag, atol=1e-1)
    @test isapprox(solve_containers.Ctheta_rotor, be2.Wtheta, atol=1e-4)

    ##### ----- Test cl, cd, and angles ----- #####
    # - Calculate Blade Element Values - #
    dt.calculate_blade_element_coefficients!(
        solve_containers.cl,
        solve_containers.cd,
        solve_containers.beta1,
        solve_containers.alpha,
        solve_containers.reynolds,
        solve_containers.mach,
        (; blade_elements..., airfoils...),
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        operating_point;
        verbose=options.verbose,
    )

    @test isapprox(solve_containers.cl, be2.cl, atol=1e-2)
    @test isapprox(solve_containers.cd, be2.cd, atol=1e-4)
    @test isapprox(pi / 2.0 .- solve_containers.beta1, be2.phi, atol=1e-3)
    @test isapprox(solve_containers.alpha, be2.alpha, atol=1e-3)

    ##### ----- Test estimated Gamr ----- #####

    Gamr_est = similar(Gamr1mat) .= 0
    dt.calculate_rotor_circulation_strengths!(
        Gamr_est, solve_containers.Cmag_rotor, blade_elements.chords, solve_containers.cl
    )
    @test isapprox(Gamr_est, be2.Gamr_est ./ B, atol=1e-1)

    ##### ----- Test relaxed Gamr ----- #####
    # with estimated Gamr test relaxed Gamr
    include(datapath * "iter2_relaxed_bgam.jl")
    Gamr2 = reformat_circulation(bgamr2, B)
    Gamr2mat = [Gamr2;;]

    deltaG_prev =
        [
            -1.37291718; -1.12498820; -0.889800072; -0.504527688; -0.191246957; -3.48181650E-02; 3.92462686E-02; 7.07832426E-02; 8.18593130E-02; 8.51513743E-02;;
        ] ./ B
    maxBGamr = MVector{1,Float64}(0.0)
    maxdeltaBGamr = MVector{1,Float64}(0.0)
    deltaG = Gamr_est - Gamr1mat

    Gamr2test = copy(Gamr1mat)

    # relax Gamr values
    Gamr2test, bladeomega, omega = dt.relax_Gamr!(
        Gamr2test, deltaG_prev, deltaG, maxBGamr, maxdeltaBGamr, blade_elements.B; test=true
    )

    @test isapprox(Gamr2test, Gamr2, atol=1e-2)

    ##### ----- Test wake velocities ----- #####
    include(datapath * "iter2_wm_wake_panels.jl")
    Wm_wake2 = reformat_wake_panel_vels(
        VMAV, nhub_interface_panels, nduct_interface_panels, pidranges
    )

    dt.calculate_wake_velocities!(
        solve_containers.Cm_wake,
        solve_containers.vz_wake,
        solve_containers.vr_wake,
        -gamw1,
        sigr2mat,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        ivw,
        operating_point.Vinf[],
    )

    Wm_wake_test = copy(solve_containers.Cm_wake)
    Wm_wake_test[idmaps.wake_panel_ids_along_casing_wake_interface] .= 0.0
    Wm_wake_test[idmaps.wake_panel_ids_along_centerbody_wake_interface] .= 0.0

    findmax(abs.(Wm_wake_test .- Wm_wake2))

    ##### ----- Test estimated gamw ----- #####
    include(datapath * "iter2_gamw_est.jl")
    gamw2_est = reformat_gamw(GAMTH2, nidx, nhub_interface_nodes, nduct_interface_nodes)

    # dt.calculate_wake_vortex_strengths!(gamw_est_test, Gamr2test, Wm_wake, inputs)

    dt.average_wake_velocities!(
        solve_containers.Cm_avg,
        solve_containers.Cm_wake,
        idmaps.wake_nodemap,
        idmaps.wake_endnodeidxs,
    )

    dt.calculate_wake_vortex_strengths!(
        solve_containers.gamw_est,
        solve_containers.Gamma_tilde,
        solve_containers.H_tilde,
        solve_containers.deltaGamma2,
        solve_containers.deltaH,
        Gamr2test,
        solve_containers.Cm_avg,
        blade_elements.B,
        operating_point.Omega,
        wakeK,
        idmaps.wake_node_sheet_be_map,
        idmaps.wake_node_ids_along_casing_wake_interface,
        idmaps.wake_node_ids_along_centerbody_wake_interface;
    )

    @test maximum(abs.(solve_containers.gamw_est .+ gamw2_est)) < 0.75

    ##### ----- Test relaxed gam ----- #####
    include(datapath * "iter2_relaxed_gamw.jl")
    gamw_relax2 = reformat_gamw(
        relaxed_GTH, nidx, nhub_interface_nodes, nduct_interface_nodes
    )
    include(datapath * "iter2_dgold.jl")
    deltag_prev = reformat_gamw(-dgold2, nidx, nhub_interface_nodes, nduct_interface_nodes)

    deltag = -(gamw2_est - gamw1)
    maxdeltagamw = MVector{1,Float64}(0.0)

    # relax gamw values
    gamw2_test, omega = dt.relax_gamw!(
        -copy(gamw1), deltag_prev, deltag, maxdeltagamw; test=true
    )

    @test maximum(abs.(gamw2_test .+ gamw_relax2)) < 1e-5

    ##### ----- Test Convergence Criteria ----- #####
    @test isapprox(maxdeltagamw[], -13.8755674)
    @test isapprox(maxdeltaBGamr[], -12.9095955, atol=1e-1)
    @test isapprox(maxBGamr[], 12.3562546)

    ##### ----- Test ----- #####
    # Update rotor blade element velocities without body influence
    include(datapath * "iter3_sigr.jl")
    sigr3mat = [sigr3;;]

    dt.calculate_induced_velocities_on_rotors!(
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        Gamr2test,
        gamw2_test,
        sigr2mat,
        @view(solve_containers.gamb[1:(idmaps.body_totnodes)]),
        @view(ivr.v_rw[:, :, 1]),
        @view(ivr.v_rr[:, :, 1]),
        @view(ivr.v_rb[:, :, 1]),
        blade_elements.B,
        blade_elements.rotor_panel_centers,
    )

    # - Get Absolute Rotor Velocities - #
    dt.reframe_rotor_velocities!(
        solve_containers.Cz_rotor,
        solve_containers.Ctheta_rotor,
        solve_containers.Cmag_rotor,
        solve_containers.vz_rotor,
        solve_containers.vtheta_rotor,
        operating_point.Vinf[],
        operating_point.Omega,
        blade_elements.rotor_panel_centers,
    )

    # - Calculate Rotor Panel Strengths - #
    dt.calculate_rotor_source_strengths!(
        sigr2mat,
        solve_containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        solve_containers.cd,
    )

    @test isapprox(sigr2mat, sigr3mat, atol=1e-2)
end
