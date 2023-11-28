using CCBlade

@testset "Rotor Loads Etc." begin

    # set up ccblade example
    Rtip = 10 / 2.0 * 0.0254  # inches to meters
    Rhub = 0.10 * Rtip
    B = 2  # number of blades
    rotor = Rotor(Rhub, Rtip, B)

    propgeom = [
        0.15 0.130 32.76
        0.20 0.149 37.19
        0.25 0.173 33.54
        0.30 0.189 29.25
        0.35 0.197 25.64
        0.40 0.201 22.54
        0.45 0.200 20.27
        0.50 0.194 18.46
        0.55 0.186 17.05
        0.60 0.174 15.97
        0.65 0.160 14.87
        0.70 0.145 14.09
        0.75 0.128 13.39
        0.80 0.112 12.84
        0.85 0.096 12.25
        0.90 0.081 11.37
        0.95 0.061 10.19
        1.00 0.041 8.99
    ]

    r = propgeom[:, 1] * Rtip
    chord = propgeom[:, 2] * Rtip
    theta = propgeom[:, 3] * pi / 180
    af = AlphaAF("test/data/naca4412.dat")
    sections = Section.(r, chord, theta, Ref(af))

    Vinf = 5.0
    Omega = 5400 * pi / 30  # convert to rad/s
    rho = 1.225
    op = simple_op.(Vinf, Omega, r, rho)

    # run ccblade
    out = solve.(Ref(rotor), sections, op)

    # run ccblade post processing
    T, Q = thrusttorque(rotor, sections, out)
    eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")

    # run DuctAPE post processing
    blade_elements = (
        Omega=Omega, B=B, chords=chord, twists=theta, rbe=r, Rtip=Rtip, Rhub=Rhub
    )
    aero = dt.get_rotor_loads(out.W, out.phi, out.cl, out.cd, blade_elements, (; Vinf, rho))
    CTdt = aero.CT
    CQdt = aero.CQ
    effdt = aero.eff
    Npdt = aero.Np
    Tpdt = aero.Tp

    # compare ccblade and DuctAPE values (should be very similar if not identical)
    @test CT == CTdt
    @test CQ == CQdt
    @test eff == effdt
end
