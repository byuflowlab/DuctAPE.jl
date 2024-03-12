println("\nSTATE ESTIMATION TESTS")

@testset "State Estimation Step Through" begin
    # some setup
    containers = (;
        Cz_rotor=zeros(3, 2),
        Ctheta_rotor=zeros(3, 2),
        Cmag_rotor=zeros(3, 2),
        cl=zeros(3, 2),
        cd=zeros(3, 2),
    )
    op = (; Vinf=1.0, rhoinf=1.225, muinf=1.81e-5, asound=341.0, Omega=[100.0, 0.0])
    rotor_panel_center = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    vz_rotor = ones(3, 2)
    vtheta_rotor = ones(3, 2)

    # - Test rotor velocity reframing - #

    dt.reframe_rotor_velocities!(
        @view(containers.Cz_rotor[:, :]),
        @view(containers.Ctheta_rotor[:, :]),
        @view(containers.Cmag_rotor[:, :]),
        vz_rotor, # state var
        vtheta_rotor, # state var
        op.Vinf[1],
        op.Omega,
        rotor_panel_center,
    )

    vthetatest1 = vtheta_rotor[:, 1] - op.Omega[1] * rotor_panel_center[:, 1]
    vthetatest2 = vtheta_rotor[:, 2] - op.Omega[2] * rotor_panel_center[:, 2]
    vthetatest = [vthetatest1 vthetatest2]

    @test all(containers.Cz_rotor .== 2.0 * ones(3, 2))
    @test all(containers.Ctheta_rotor .== vthetatest)
    @test all(containers.Cmag_rotor .== sqrt.(vthetatest .^ 2 .+ (2.0 * ones(3, 2)) .^ 2))

    # more setup
    dfdcparam = dt.c4b.DFDCairfoil(
        0.0, 1.5, -1.0, 6.28, 0.5, 0.2, 0.012, 0.1, 0.005, 0.0, 200000.0, 0.35, 0.7
    )

    blade_elements = (;
        B=[5, 5],
        inner_fraction=0.5 * ones(3, 2),
        solidity=0.5 * ones(3, 2),
        stagger=0.5 * ones(3, 2),
        inner_airfoil=fill(dfdcparam, 3, 2),
        outer_airfoil=fill(dfdcparam, 3, 2),
        fliplift=zeros(3, 2),
        chords=ones(3, 2),
        twists=zeros(3, 2),
        rotor_panel_center=rotor_panel_center,
    )

    # - Test cl and cd calcs - #
    dt.calculate_blade_element_coefficients!(
        @view(containers.cl[:, :]),
        @view(containers.cd[:, :]),
        blade_elements,
        ones(3, 2),
        ones(3, 2),
        ones(3, 2),
        op;
        post=false,
        verbose=false,
    )

    @test all(containers.cl .== -1.5630840108818207)
    @test all(containers.cd .== 2.166343708285816)

    # - Test Circulation calcs - #
    Gamr = zeros(3, 2)
    dt.calculate_rotor_circulation_strengths!(
        Gamr, containers.Cmag_rotor, blade_elements.chords, 2.0 * ones(size(Gamr))
    )

    @test all(Gamr .== 0.5 .* containers.Cmag_rotor .* blade_elements.chords * 2.0)

    # - Test source strength calcs - #
    sigr = zeros(4, 2)
    dt.calculate_rotor_source_strengths!(
        sigr,
        containers.Cmag_rotor,
        blade_elements.chords,
        blade_elements.B,
        2.0 * ones(size(sigr)),
    )

    @test all(
        sigr[1, :] .==
        B ./ (4.0 * pi) * containers.Cmag_rotor[1, :] .* blade_elements.chords[1, :] .* 2.0,
    )
    @test all(
        sigr[2, :] .==
        B ./ (4.0 * pi) * (
            containers.Cmag_rotor[1, :] .* blade_elements.chords[1, :] .* 2.0 +
            containers.Cmag_rotor[2, :] .* blade_elements.chords[2, :] .* 2.0
        ) / 2,
    )

    # - test wake velocity averaging - #
    Cm_avg = zeros(3)
    Cm_wake = [1.0; 2.0]
    endnodeidxs = [1.0; 3.0]
    nodemap = [1 2; 2 3]
    dt.average_wake_velocities!(Cm_avg, Cm_wake, nodemap, endnodeidxs)

    @test Cm_avg == [1.0, 1.5, 2.0]

    # - test wake strength calcs - #
    wakeK = ones(16)
    gamw = zeros(16)
    wake_node_ids_along_casing_wake_interface = nothing
    wake_node_ids_along_centerbody_wake_interface = nothing
    rotorwakenodeid = [
        1 1
        1 1
        1 2
        1 2
        2 1
        2 1
        2 2
        2 2
        3 1
        3 1
        3 2
        3 2
        4 1
        4 1
        4 2
        4 2
    ]
    dt.calculate_wake_vortex_strengths!(
        gamw,
        ones(3, 2),
        ones(16),
        blade_elements.B,
        op.Omega,
        wakeK,
        rotorwakenodeid,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_centerbody_wake_interface;
    )
    @test all(gamw[5:12] .== 0.0)
    @test !any(gamw[1:4] .== 0.0)
    @test !any(gamw[13:16] .== 0.0)
    #TODO: may need a better test here at some point.

    # - test Induced velocities on rotors - #
    Vz_est = zeros(3, 2)
    Vtheta_est = zeros(3, 2)
    Gamr = ones(3, 2)
    gamw = ones(16)
    sigr = ones(4, 2)
    gamb = ones(3)
    ivr = (; v_rb=ones(6, 3, 2), v_rw=ones(6, 16, 2), v_rr=ones(6, 8, 2))
    dt.calculate_induced_velocities_on_rotors!(
        @view(Vz_est[:, :]),
        @view(Vtheta_est[:, :]),
        Gamr,
        gamw,
        sigr,
        gamb,
        ivr,
        blade_elements.B,
        blade_elements.rotor_panel_center,
    )

    @test all(Vz_est .== 27.0)
    @test all(
        Vtheta_est .== [
            0.3978873577297384 1.1936620731892151
            0.1989436788648692 0.5968310365946076
            0.1326291192432461 0.3978873577297383
        ],
    )

    vzt, vtt = dt.calculate_induced_velocities_on_rotors(
        blade_elements.B,
        blade_elements.rotor_panel_center,
        Gamr,
        ivr.v_rw[:, :, 1],
        gamw,
        ivr.v_rr[:, :, 1],
        sigr,
        ivr.v_rb[:, :, 1],
        gamb,
    )

    @test all(vzt .== Vz_est)
    @test all(vtt .== Vtheta_est)

    # - test wake velocity calcs - #
    Cm_wake = zeros(12)
    vz_wake = ones(12)
    vr_wake = ones(12)
    Vinf = 1.0
    dt.reframe_wake_velocities!(@view(Cm_wake[:]), vz_wake, vr_wake, Vinf)
    @test all(Cm_wake .== sqrt.((vz_wake .+ Vinf) .^ 2 .+ vr_wake .^ 2))

    ivw = (; v_wb=ones(12, 3, 2), v_wr=ones(12, 8, 2), v_ww=ones(12, 16, 2))
    dt.calculate_wake_velocities!(
        @view(Cm_wake[:]), @view(vz_wake[:]), @view(vr_wake[:]), gamw, sigr, gamb, ivw, Vinf
    )
    @test all(Cm_wake .== 38.897300677553446)
end
