"""
W : Inflow Magnitude
phi : Inflow angle
"""
function get_rotor_loads(W, phi, cl, cd, blade_element, fs)

    # rename for convenience
    cphi = cos.(phi)
    sphi = sin.(phi)

    # resolve lift and drag into normal and tangential coefficients
    cn = cl .* cphi .- cd .* sphi
    ct = cl .* sphi .+ cd .* cphi

    # get the normal and tangential loads per unit length N' and T'
    Np = cn .* 0.5 .* fs.rho .* W .^ 2 .* blade_element.chords
    Tp = ct .* 0.5 .* fs.rho .* W .^ 2 .* blade_element.chords

    ## -- Integrate Loads to get Thurst and Torque
    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    rfull = [blade_element.Rhub; blade_element.rbe; blade_element.Rtip]
    Npfull = [0.0; Np; 0.0]
    Tpfull = [0.0; Tp; 0.0]

    # thrust and torqe distributions
    thrust = Npfull
    torque = Tpfull .* rfull

    # integrate Thrust and Torque (trapezoidal)
    T = blade_element.B * fm.trapz(rfull, thrust)
    Q = blade_element.B * fm.trapz(rfull, torque)
    P = Q * blade_element.Omega

    ## -- Get Non-dimensional versions -- ##
    n = blade_element.Omega / (2 * pi)
    D = 2 * blade_element.Rtip

    if T < 0
        eff = 0.0  # creating drag not thrust
    else
        eff = T * fs.Vinf / P
    end
    CT = T / (fs.rho * n^2 * D^4)
    CQ = Q / (fs.rho * n^2 * D^5)

    return (; Np, Tp, T, Q, CT, CQ, eff, cn, ct)
end

"""
"""
function states_to_outputs_rotor_only(states, params)
    Gamma, gamma_theta, sigma = extract_rotor_states(states, params)

    wake_vortex_strengths = repeat(gamma_theta; inner=(1, params.nxwake))

    TF = eltype(Gamma)

    # - get the induced velocities at the rotor plane - #
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        params.blade_elements,
        Gamma,
        params.vx_rw,
        params.vr_rw,
        wake_vortex_strengths,
        params.vx_rr,
        params.vr_rr,
        sigma,
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ params.Vinf
    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wθ = vtheta_rotor .- params.blade_elements[1].Omega .* params.blade_elements[1].rbe

    Wm = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # - Get the inflow magnitude at the rotor as the combination of all the components - #
    W = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wθ .^ 2)

    nbe = length(params.blade_elements[1].rbe)
    r = zeros(TF, nbe)
    chord = zeros(TF, nbe)
    twist = zeros(TF, nbe)
    phi = zeros(TF, nbe)
    alpha = zeros(TF, nbe)
    affrac = zeros(TF, nbe)
    cl = zeros(TF, nbe)
    clin = zeros(TF, nbe)
    clout = zeros(TF, nbe)
    cd = zeros(TF, nbe)
    cdin = zeros(TF, nbe)
    cdout = zeros(TF, nbe)

    for ir in 1:(nbe)
        # extract blade element properties
        B = params.blade_elements[1].B # number of blades
        chord[ir] = params.blade_elements[1].chords[ir] # chord length
        twist[ir] = params.blade_elements[1].twists[ir] # twist
        r[ir] = params.blade_elements[1].rbe[ir] # radius


        # calculate angle of attack
        # phi[ir] = atan(Wm[ir], -Wθ[ir])
        phi[ir] = atan(Wm[ir], -Wθ[ir])
        alpha[ir] = twist[ir] - phi[ir]

        affrac[ir] = params.blade_elements[1].inner_fraction[ir]

        # look up lift and drag data for the nearest two input sections
        clin[ir], cdin[ir] = search_polars(
            params.blade_elements[1].inner_airfoil[ir], alpha[ir]
        )
        clout[ir], cdout[ir] = search_polars(
            params.blade_elements[1].outer_airfoil[ir], alpha[ir]
        )
        # linearly interpolate between those two values at your blade element location
        cl[ir] = fm.linear([0.0; 1.0], [clin[ir], clout[ir]], affrac[ir])
        cd[ir] = fm.linear([0.0; 1.0], [cdin[ir], cdout[ir]], affrac[ir])
    end

    return (;
        r,
        chord,
        twist,
        vx_rotor,
        vr_rotor,
        vtheta_rotor,
        # Wm=Wx_rotor,
        Wm,
        Wθ,
        W,
        phi,
        alpha,
        affrac,
        clin,
        clout,
        cl,
        cdin,
        cdout,
        cd,
    )
end