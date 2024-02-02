println("\nPOST PROCESSING TESTS")

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
    rhoinf = 1.226
    Omega = 8000 * pi / 30

    # columns: B*Gamma,Wtheta,dr,section thrust
    Wtheta = dfdc_thrust_inviscid[:, 2]
    Gamma_tilde = reshape(dfdc_thrust_inviscid[:, 1], (length(Wtheta), 1))
    rotor_panel_length = dfdc_thrust_inviscid[:, 3]
    dfdc_tinv_dist = dfdc_thrust_inviscid[:, 4]

    tinv, tinv_dist = dt.inviscid_rotor_trust(
        Wtheta, Gamma_tilde, rotor_panel_length, rhoinf
    )

    @test all(isapprox.(dfdc_tinv_dist, tinv_dist; atol=1e-6))

    # columns: W,chord,cd,Wx,dr, section thrust
    Wmag_rotor = dfdc_thrust_viscous[:, 1]
    chord = dfdc_thrust_viscous[:, 2]
    cd = dfdc_thrust_viscous[:, 3]
    Wx_rotor = reshape(dfdc_thrust_viscous[:, 4], (length(Wmag_rotor), 1))
    rotor_panel_length = dfdc_thrust_viscous[:, 5]
    dfdc_tvis_dist = dfdc_thrust_viscous[:, 6]

    tvis, tvis_dist = dt.viscous_rotor_thrust(
        Wx_rotor, Wmag_rotor, B, chord, rotor_panel_length, cd, rhoinf
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
        Wx_rotor, rpc, rpl, Gamma_tilde, rhoinf
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
        Wtheta_rotor, Wmag_rotor, B, chord, rpc, rpl, cdrag, rhoinf
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
    dfdc_d_thrust = -0.5 * rhoinf * Vref^2 * dfdc_duct_thrust[end, 7]
    controlpoint = zeros(2, length(panel_radial_position))
    controlpoint[2, :] = panel_radial_position

    ductpanel = (; nbodies=1, endpanelidxs=[1; length(panel_length);;], normal=xnormal' .* 1.0, influence_length=panel_length, controlpoint)

    duct_thrust, _ = dt.forces_from_pressure(
        cpleft,  cpright, ductpanel; rhoinf=rhoinf, Vref=Vref
    )

    @test isapprox(duct_thrust[1], dfdc_d_thrust; atol=1e-4)

    # COLUMNS:IC,CPL(IC),CPR(IC),xnromal, DSC(IC),YC(IC)
    cpleft = dfdc_hub_thrust[:, 2]
    cpright = dfdc_hub_thrust[:, 3]
    xnormal = dfdc_hub_thrust[:, 4]
    panel_length = dfdc_hub_thrust[:, 5]
    panel_radial_position = dfdc_hub_thrust[:, 6]
    dfdc_h_thrust = -0.5 * rhoinf * Vref^2 * dfdc_hub_thrust[end, 7]
    controlpoint = zeros(2, length(panel_radial_position))
    controlpoint[2, :] = panel_radial_position

    hubpanel = (; nbodies=1, endpanelidxs=[1; length(panel_length);;], normal=xnormal' .* 1.0, influence_length=panel_length, controlpoint)

    hub_thrust, _ = dt.forces_from_pressure(
        cpleft,  cpright, hubpanel; rhoinf=rhoinf, Vref=Vref
    )

    @test isapprox(hub_thrust[1], dfdc_h_thrust; atol=1e-4)

end

#TODO: update this test to use DuctAPE.c4b
# @testset "Rotor Loads" begin

#     # set up ccblade example
#     Rtip = 10 / 2.0 * 0.0254  # inches to meters
#     Rhub = 0.10 * Rtip
#     B = 2  # number of blades
#     rotor = Rotor(Rhub, Rtip, B)

#     propgeom = [
#         0.15 0.130 32.76
#         0.20 0.149 37.19
#         0.25 0.173 33.54
#         0.30 0.189 29.25
#         0.35 0.197 25.64
#         0.40 0.201 22.54
#         0.45 0.200 20.27
#         0.50 0.194 18.46
#         0.55 0.186 17.05
#         0.60 0.174 15.97
#         0.65 0.160 14.87
#         0.70 0.145 14.09
#         0.75 0.128 13.39
#         0.80 0.112 12.84
#         0.85 0.096 12.25
#         0.90 0.081 11.37
#         0.95 0.061 10.19
#         1.00 0.041 8.99
#     ]

#     r = propgeom[:, 1] * Rtip
#     chord = propgeom[:, 2] * Rtip
#     theta = propgeom[:, 3] * pi / 180
#     af = AlphaAF("test/data/naca4412.dat")
#     sections = Section.(r, chord, theta, Ref(af))

#     Vinf = 5.0
#     Omega = 5400 * pi / 30  # convert to rad/s
#     rho = 1.225
#     op = simple_op.(Vinf, Omega, r, rho)

#     # run ccblade
#     out = solve.(Ref(rotor), sections, op)

#     # run ccblade post processing
#     T, Q = thrusttorque(rotor, sections, out)
#     eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")

#     # run DuctAPE post processing
#     blade_elements = (
#         Omega=Omega, B=B, chords=chord, twists=theta, rbe=r, Rtip=Rtip, Rhub=Rhub
#     )
#     aero = dt.get_rotor_loads(out.W, out.phi, out.cl, out.cd, blade_elements, (; Vinf, rho))
#     CTdt = aero.CT
#     CQdt = aero.CQ
#     effdt = aero.eff
#     Npdt = aero.Np
#     Tpdt = aero.Tp

#     # compare ccblade and DuctAPE values (should be very similar if not identical)
#     @test CT == CTdt
#     @test CQ == CQdt
#     @test eff == effdt
# end
