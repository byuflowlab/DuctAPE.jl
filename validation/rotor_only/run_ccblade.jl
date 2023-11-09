function run_ccblade(Vinf; airfoil="test/data/naca4412.dat")

    # set up ccblade example
    Rtip = 10 / 2.0 * 0.0254  # inches to meters
    Rhub = 0.10 * Rtip
    B = 2  # number of blades
    rotor = Rotor(Rhub, Rtip, B, tip=nothing)

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
        # 1.00 0.041 8.99
    ]

    r = propgeom[:, 1] * Rtip
    chord = propgeom[:, 2] * Rtip
    theta = propgeom[:, 3] * pi / 180
    # af = AlphaAF("test/data/naca4412.dat")
    # af = AlphaAF("test/data/xrotor_af_test.dat")
    af = AlphaAF(airfoil)
    sections = Section.(r, chord, theta, Ref(af))

    Omega = 5400 * pi / 30  # convert to rad/s
    rho = 1.225
    op = simple_op.(Vinf, Omega, r, rho)

    # run ccblade
    out = ccb.solve.(Ref(rotor), sections, op)

    # run ccblade post processing
    T, Q = thrusttorque(rotor, sections, out)
    eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")

    circ = 0.5.*out.cl.*out.W.*chord

    return (;eff, CT, CQ, out, circ, r)
end
