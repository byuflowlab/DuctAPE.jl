println("\nITERATION STEP THROUGH TESTS")

@testset "Iteration Steps" begin
    ##### ----- Set up ----- #####
    datapath = "data/dfdc_init_iter1/"
    include(datapath * "ductape_parameters_OLD.jl")
    include(datapath * "ductape_formatted_dfdc_geometry_OLD.jl")

    # generate inputs
    inputs = dt.precomputed_inputs(
        system_geometry,
        rotorstator_parameters, #vector of named tuples
        freestream,
        reference_parameters;
        debug=false,
    )

    # - Parameters - #
    B = 5 #there are 5 blades

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

    # get wake sheet length
    num_wake_z_panels = num_wake_z_nodes - 1

    # get number of nodes and panels on bodies interface with wake
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
        inputs.gamb,
        inputs.A_bb,
        inputs.b_bf,
        -gamw1,
        inputs.A_bw,
        inputs.A_pw,
        sigr2mat,
        inputs.A_br,
        inputs.A_pr,
        inputs.LHS;
        post=false,
    )

    # - test body strengths against DFDC written values - #
    # get the body strengths from DFDC
    include(datapath * "iter2_gamb.jl")
    gamb2 = reformat_sol(res2, nidx)
    @test maximum(inputs.gamb .+ gamb2) < 0.5

    ##### ----- Test blade element values----- #####
    include(datapath * "iter2_blade_element_values.jl")
    b2 = reformat_rotor_velocities(rotor_vels2)
    include(datapath * "iter2_extended_blade_element_values.jl")
    be2 = reformat_blade_elements(extended_blade_element_values2)

    _, vr_rotor, _, Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = dt.calculate_rotor_velocities(
        Gamr1mat,
        -gamw1,
        sigr2mat,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)],
        inputs,
    )

    @test isapprox(vr_rotor, b2.Wr, atol=1e-1)
    @test isapprox(Wm_rotor, b2.Wm, atol=1e-1)
    @test isapprox(Wz_rotor, be2.Wz, atol=1e-1)
    @test isapprox(Wmag_rotor, be2.Wmag, atol=1e-1)
    @test isapprox(Wtheta_rotor, be2.Wtheta, atol=1e-4)

    ##### ----- Test cl, cd, and angles ----- #####
    cl, cd, phi, alpha = dt.calculate_blade_element_coefficients(
        inputs.blade_elements,
        Wz_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        inputs.freestream;
        post=true,
    )

    @test isapprox(cl, be2.cl, atol=1e-3)
    @test isapprox(cd, be2.cd, atol=1e-5)
    @test isapprox(pi / 2.0 .- phi, be2.phi, atol=1e-3)
    @test isapprox(alpha, be2.alpha, atol=1e-3)

    ##### ----- Test estimated Gamr ----- #####

    Gamr_est = similar(Gamr1mat) .= 0
    dt.calculate_rotor_circulation_strengths!(
        Gamr_est, Wmag_rotor, inputs.blade_elements, cl
    )
    @test isapprox(Gamr_est, be2.Gamr_est ./ B, atol=1e-2)

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
        Gamr2test,
        deltaG_prev,
        deltaG,
        maxBGamr,
        maxdeltaBGamr,
        inputs.blade_elements.B;
        test=true,
    )

    @test isapprox(Gamr2test, Gamr2, atol=1e-3)

    ##### ----- Test wake velocities ----- #####
    include(datapath * "iter2_wm_wake_panels.jl")
    Wm_wake2 = reformat_wake_panel_vels(
        VMAV, nhub_interface_panels, nduct_interface_panels, pidranges
    )

    Wm_wake = dt.calculate_wake_velocities(
        -gamw1, sigr2mat, inputs.gamb[1:(inputs.body_vortex_panels.totnode)], inputs
    )

    Wm_wake_test = copy(Wm_wake)
    Wm_wake_test[inputs.ductwakeinterfaceid] .= 0.0
    Wm_wake_test[inputs.hubwakeinterfaceid] .= 0.0

    @test isapprox(Wm_wake_test, Wm_wake2, atol=5e-2)

    ##### ----- Test estimated gamw ----- #####
    include(datapath * "iter2_gamw_est.jl")
    gamw2_est = reformat_gamw(GAMTH2, nidx, nhub_interface_nodes, nduct_interface_nodes)

    gamw_est_test = copy(gamw2_est) .= 0
    dt.calculate_wake_vortex_strengths!(gamw_est_test, Gamr2test, Wm_wake, inputs)

    @test maximum(gamw_est_test .+ gamw2_est) < 0.11

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

    @test maximum(gamw2_test .+ gamw_relax2) < 1e-5

    ##### ----- Test Convergence Criteria ----- #####
    @test isapprox(maxdeltagamw[], -13.8755674)
    @test isapprox(maxdeltaBGamr[], -12.9095955, atol=1e-2)
    @test isapprox(maxBGamr[], 12.3562546)

    ##### ----- Test ----- #####
    #TODO: with everything done, test sigr to be used for next iteration.
    # Update rotor blade element velocities without body influence
    include(datapath * "iter3_sigr.jl")
    sigr3mat = [sigr3;;]

    _, _, _, _, _, _, Wmag_rotor = dt.calculate_rotor_velocities(
        Gamr2test, gamw2_test, sigr2mat, inputs
    )

    # update sigr in place
    dt.calculate_rotor_source_strengths!(
        sigr2mat, Wmag_rotor, inputs.blade_elements, cd, inputs.freestream.rhoinf
    )

    @test isapprox(sigr2mat, sigr3mat, atol=1e-2)
end
