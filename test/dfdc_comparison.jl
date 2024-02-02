println("\nMISC TESTS")

@testset "Body Surface Pressure" begin
    datapath = "data/dfdc_init_iter1/"
    include(datapath * "ductape_parameters.jl")
    include(datapath * "ductape_formatted_dfdc_geometry.jl")

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

    # - read in functions to reformat - #
    include(datapath * "reformatting_functions.jl")
    include(datapath * "iter2_sigr.jl")
    include(datapath * "iter1_relaxed_gamw.jl")
    gamw1 = reformat_gamw(relaxed_GTH, nidx, nhub_interface_nodes, nduct_interface_nodes)
    include(datapath * "iter1_relaxed_bgam.jl")
    Gamr1 = reformat_circulation(relaxed_BGAM, B)

    # - Body Pressures - #
    include(datapath * "iter2_body_cp_withdeltas.jl")
    duct_cpR2, hub_cpR2 = reformat_body_cp(body_cpR2, pidx)

    Gamr = zeros(length(Gamr1), 1) .= copy(Gamr1)
    sigr = zeros(length(sigr2), 1) .= copy(sigr2)
    gamw = -copy(gamw1)
    outs = dt.post_process(Gamr, sigr, gamw, inputs; write_outputs=false)

    duct_cp = [outs.bodies.cp_casing_out; outs.bodies.cp_nacelle_out]
    centerbody_cp = outs.bodies.cp_centerbody_out
    @test all(isapprox.(duct_cp, reverse(duct_cpR2), atol=1e-2))
    @test all(isapprox.(centerbody_cp, reverse(hub_cpR2), atol=1e-2))
end

@testset "Iteration Steps" begin
    ##### ----- Set up ----- #####
    datapath = "data/dfdc_init_iter1/"
    include(datapath * "ductape_parameters.jl")
    include(datapath * "ductape_formatted_dfdc_geometry.jl")

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
        inputs.RHS;
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
        gamw1,
        sigr2mat,
        inputs.gamb[1:(inputs.body_vortex_panels.totnode)],
        inputs,
    )

    @test isapprox(vr_rotor, b2.Wr, atol=1e-1)
    @test isapprox(Wm_rotor, b2.Wm, atol=1e-1)
    @test isapprox(Wz_rotor, be2.Wz, atol=1e-1)
    @test isapprox(Wmag_rotor, be2.Wmag, atol=1e-1)
    @test isapprox(Wtheta_rotor, be2.Wtheta, atol=1e-4)

    ##### ----- Test ----- #####
    #TODO: with updated rotor velocities, test blade element values (cl, cd, phi, alpha)
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

    ##### ----- Test ----- #####
    #TODO: with updated blade values, test estimated Gamr

    Gamr_est = similar(Gamr1mat) .= 0
    dt.calculate_rotor_circulation_strengths!(
        Gamr_est, Wmag_rotor, inputs.blade_elements, cl
    )
    @test isapprox(Gamr_est, be2.Gamr_est ./ B, atol=1e-2)

    ##### ----- Test ----- #####
    #TODO; with estimated Gamr test relaxed Gamr and convergence criteria
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

    ##### ----- Test ----- #####
    #TODO: with relaxed Gamr, text jump terms

    #
    #TODO: you are here TODO: you are here TODO: you are here
    #

    ##### ----- Test ----- #####
    #TODO; with updated body strengths, test wake velocities

    ##### ----- Test ----- #####
    #TODO; with updated wake velocities, test estimated gamw

    ##### ----- Test ----- #####
    #TODO: with estimated gamw, test relaxed gamw and convergence criteria

end
