# Thrust Tests
@testset "Post-Processing" begin

    #---------------------------------#
    #           ROTOR THRUST          #
    #---------------------------------#

    # load in the dfdc data for comparison
    include("data/dfdc_dstestr/DFDC_THRUST_INVISCID.jl")
    include("data/dfdc_dstestr/DFDC_THRUST_VISCOUS.jl")

    B = 5
    Vinf = 0.0
    Vref = 50.0
    rho = 1.226
    Omega = 8000 * pi / 30

    # columns: B*Gamma,Wtheta,dr,section thrust
    Wtheta = dfdc_thrust_inviscid[:, 2]
    Gamma_tilde = reshape(dfdc_thrust_inviscid[:, 1], (length(Wtheta), 1))
    rotor_panel_length = dfdc_thrust_inviscid[:, 3]
    dfdc_tinv_dist = dfdc_thrust_inviscid[:, 4]

    tinv, tinv_dist = dt.inviscid_rotor_trust(Wtheta, Gamma_tilde, rotor_panel_length, rho)

    @test all(isapprox.(dfdc_tinv_dist, tinv_dist; atol=1e-6))

    # columns: W,chord,cd,Wx,dr, section thrust
    Wmag_rotor = dfdc_thrust_viscous[:, 1]
    chord = dfdc_thrust_viscous[:, 2]
    cd = dfdc_thrust_viscous[:, 3]
    Wx_rotor = reshape(dfdc_thrust_viscous[:, 4], (length(Wmag_rotor), 1))
    rotor_panel_length = dfdc_thrust_viscous[:, 5]
    dfdc_tvis_dist = dfdc_thrust_viscous[:, 6]

    tvis, tvis_dist = dt.viscous_rotor_thrust(
        Wx_rotor, Wmag_rotor, B, chord, rotor_panel_length, cd, rho
    )

    @test all(isapprox.(dfdc_tvis_dist, tvis_dist; atol=1e-6))

    #---------------------------------#
    #           ROTOR TORQUE          #
    #---------------------------------#

    # load in the dfdc data for comparison
    include("data/dfdc_dstestr/DFDC_TORQUE_INVISCID.jl")
    include("data/dfdc_dstestr/DFDC_TORQUE_VISCOUS.jl")

    # columns: B*Gamma,Wx,rpl,rpc,section torque
    Wx_rotor = dfdc_torque_inviscid[:, 2]
    Gamma_tilde = reshape(dfdc_torque_inviscid[:, 1], (length(Wx_rotor), 1))
    rpl = dfdc_torque_inviscid[:, 3]
    rpc = dfdc_torque_inviscid[:, 4]
    dfdc_qinv_dist = dfdc_torque_inviscid[:, 5]

    rotor_inviscid_torque, rotor_inviscid_torque_dist = dt.inviscid_rotor_torque(
        Wx_rotor, rpc, rpl, Gamma_tilde, rho
    )

    @test all(isapprox.(dfdc_qinv_dist, rotor_inviscid_torque_dist; atol=1e-6))

    # columns: W,chord,cd,Wt,rpl,rpc,section torque
    Wmag_rotor = dfdc_torque_viscous[:, 1]
    chord = dfdc_torque_viscous[:, 2]
    cdrag = dfdc_torque_viscous[:, 3]
    Wtheta_rotor = reshape(dfdc_torque_viscous[:, 4], (length(Wmag_rotor), 1))
    rpl = dfdc_torque_viscous[:, 5]
    rpc = dfdc_torque_viscous[:, 6]
    dfdc_qvis_dist = dfdc_torque_viscous[:, 7]

    rotor_viscous_torque, rotor_viscous_torque_dist = dt.viscous_rotor_torque(
        Wtheta_rotor, Wmag_rotor, B, chord, rpc, rpl, cdrag, rho
    )

    @test all(isapprox.(dfdc_qvis_dist, rotor_viscous_torque_dist; atol=1e-6))

    #---------------------------------#
    #           ROTOR POWER           #
    #---------------------------------#

    pinv = dt.inviscid_rotor_power(rotor_inviscid_torque, Omega)

    pvis = dt.viscous_rotor_power(rotor_viscous_torque, Omega)

    @test pvis == rotor_viscous_torque * Omega
    @test pinv == rotor_inviscid_torque * Omega

    #---------------------------------#
    #           DUCT THRUST           #
    #---------------------------------#

    # load in the dfdc data for comparison
    include("data/dfdc_dstestr/DFDC_DUCT_FORCES.jl")
    include("data/dfdc_dstestr/DFDC_HUB_FORCES.jl")

    # COLUMNS:IC,CPL(IC),CPR(IC),xnromal, DSC(IC),YC(IC)
    cpleft = dfdc_duct_thrust[:, 2]
    cpright = dfdc_duct_thrust[:, 3]
    xnormal = dfdc_duct_thrust[:, 4]
    panel_length = dfdc_duct_thrust[:, 5]
    panel_radial_position = dfdc_duct_thrust[:, 6]
    dfdc_d_thrust = -0.5 * rho * Vref^2 * dfdc_duct_thrust[end, 7]

    ductpanel = (;
        panel_normal=xnormal,
        panel_length,
        panel_center=[zeros(length(panel_radial_position)) panel_radial_position],
    )

    duct_thrust, _ = dt.forces_from_pressure(
        cpleft .- cpright, ductpanel; rho=rho, Vref=Vref
    )

    @test isapprox(-duct_thrust, dfdc_d_thrust; atol=1e-4)

    # COLUMNS:IC,CPL(IC),CPR(IC),xnromal, DSC(IC),YC(IC)
    cpleft = dfdc_hub_thrust[:, 2]
    cpright = dfdc_hub_thrust[:, 3]
    xnormal = dfdc_hub_thrust[:, 4]
    panel_length = dfdc_hub_thrust[:, 5]
    panel_radial_position = dfdc_hub_thrust[:, 6]
    dfdc_h_thrust = -0.5 * rho * Vref^2 * dfdc_hub_thrust[end, 7]

    hubpanel = (;
        panel_normal=xnormal,
        panel_length,
        panel_center=[zeros(length(panel_radial_position)) panel_radial_position],
    )

    hub_thrust, _ = dt.forces_from_pressure(cpleft.-cpright, hubpanel; rho=rho, Vref=Vref)

    @test isapprox(-hub_thrust, dfdc_h_thrust; atol=1e-4)

    display(hub_thrust)
    display(duct_thrust)
    display(hub_thrust + duct_thrust)
end
