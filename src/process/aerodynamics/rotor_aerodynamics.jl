"""
    calculate_induced_velocities_on_rotors(blade_elements, Gamr, vz_rb, vr_rb, gamb,
        vz_rw, vr_rw, gamw)

Calculate axial, radial, and tangential induced velocity on the rotors.

`Gamr::Matrix{Float}` : rotor circulation values (must be a matrix, even if only 1 rotor is being used). Size is number of blade elements by number of rotors.
`vz_rw::Matrix{Matrix{Float}}` : unit axial induced velocities of wakes on rotors. Must be a matrix of size number of rotors by number of wakes made up of matrices of size number of blade elements by number of wake "panels."
`gamw::Matrix{Float}` : wake circulation values (must be a matrix). Size is number of wake sheets by number of "panels" per sheet.
"""
function calculate_induced_velocities_on_rotors(
    Gamr, gamw, sigr, gamb, vz_rw, vz_rr, vz_rb, B, rotor_panel_centers; post=false
)

    # problem dimensions
    nbe, nrotor = size(Gamr) # number of rotors

    # initialize outputs
    vz_rotor = similar(Gamr) .= 0 # axial induced velocity
    vtheta_rotor = similar(Gamr) .= 0 # tangential induced velocity
    return calculate_induced_velocities_on_rotors!(
        vz_rotor,
        vtheta_rotor,
        Gamr,
        gamw,
        sigr,
        gamb,
        vz_rw,
        vz_rr,
        vz_rb,
        B,
        rotor_panel_centers;
        post=post,
    )
end

"""
- `ivr:NamedTuple` : induced velocities on rotor
"""
function calculate_induced_velocities_on_rotors!(
    vz_rotor,
    vtheta_rotor,
    Gamr,
    gamw,
    sigr,
    gamb,
    vz_rw,
    vz_rr,
    vz_rb,
    B,
    rotor_panel_centers;
    post=false,
)

    # - Clear out inputs - #
    vz_rotor .= 0
    vtheta_rotor .= 0

    # problem dimensions
    nbe, nrotor = size(Gamr) # number of rotors

    if post

        # initialize outputs
        vzb_rotor = similar(Gamr) .= 0 # axial induced velocity

        vzw_rotor = similar(Gamr) .= 0 # axial induced velocity

        vzr_rotor = similar(Gamr) .= 0 # axial induced velocity
    end

    # loop through each affected rotor
    for irotor in 1:nrotor
        irotorrange = (nbe * (irotor - 1) + 1):(nbe * irotor)
        # add body induced velocities
        if gamb != nothing
            @views vz_rotor[:, irotor] .+= vz_rb[irotorrange, :] * gamb
            if post
                @views vzb_rotor[:, irotor] .+= vz_rb[irotorrange, :] * gamb
            end
        end

        # add wake induced velocities
        @views vz_rotor[:, irotor] .+= vz_rw[irotorrange, :] * gamw

        if post
            @views vzw_rotor[:, irotor] .+= vz_rw[irotorrange, :] * gamw[:]
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vz_rotor[:, irotor] .+=
                vz_rr[irotorrange, ((nbe + 1) * (jrotor - 1) + 1):((nbe + 1) * jrotor)] *
                sigr[:, jrotor]

            if post
                @views vzr_rotor[:, irotor] .+=
                    vz_rr[irotorrange, (nbe * (jrotor - 1) + 1):(nbe * jrotor)] *
                    sigr[:, jrotor]
            end
        end

        # add self-induced tangential velocity
        r = @view(rotor_panel_centers[:, irotor])
        @views vtheta_rotor[:, irotor] .+= B[irotor] .* Gamr[:, irotor] ./ (4.0 * pi * r)

        # add induced tangential velocity from upstream rotors
        for jrotor in 1:(irotor - 1)
            r = @view(rotor_panel_centers[:, jrotor])
            @views vtheta_rotor[:, irotor] .+=
                B[irotor] .* Gamr[:, jrotor] ./ (2.0 * pi * r)
        end
    end

    # return raw induced velocities
    if post
        return vz_rotor, vtheta_rotor, vzb_rotor, vzw_rotor, vzr_rotor

    else
        return vz_rotor, vtheta_rotor
    end
end

"""
"""
function reframe_rotor_velocities!(
    Cz_rotor,
    Ctheta_rotor,
    Cmag_rotor,
    vz_rotor,
    vtheta_rotor,
    Vinf,
    Omega,
    rotor_panel_centers,
)

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    @. Cz_rotor = vz_rotor + Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Ctheta_rotor .= vtheta_rotor

    for (io, O) in enumerate(Omega)
        @view(Ctheta_rotor[:, io]) .-= O .* @view(rotor_panel_centers[:, io])
    end

    # Get the inflow magnitude at the rotor as the combination of axial and tangential components
    @. Cmag_rotor = sqrt(Cz_rotor^2 + Ctheta_rotor^2)

    return Cz_rotor, Ctheta_rotor, Cmag_rotor
end

"""
"""
function calculate_rotor_circulation_strengths!(
    Gamr, Wmag_rotor, chords, cls, blade_elements
)
    for (G, c, cl, W, is_stator) in zip(
        eachcol(Gamr),
        eachcol(chords),
        eachcol(cls),
        eachcol(Wmag_rotor),
        blade_elements.is_stator,
    )
        # calculate vortex strength
        @. G = 0.5 * W * c * cl
        # flip if stator
        if !iszero(is_stator)
            G .*= -1.0
        end
    end

    return Gamr
end

# function calculate_rotor_source_strengths!(sigr, Wmag_rotor, blade_elements, rho)
function calculate_rotor_source_strengths!(sigr, Wmag_rotor, chords, B, cd)

    # Loop through rotors
    for irotor in axes(sigr, 2)

        # Loop through blade panel nodes
        for inode in axes(sigr, 1)

            # Get average of blade element values at nodes
            if inode == 1
                W = Wmag_rotor[inode, irotor]
                c = chords[inode, irotor]
                d = cd[inode, irotor]

            elseif inode == size(sigr, 1)
                W = Wmag_rotor[inode - 1, irotor]
                c = chords[inode - 1, irotor]
                d = cd[inode - 1, irotor]

            else
                W = Wmag_rotor[(inode - 1):inode, irotor]
                c = chords[(inode - 1):inode, irotor]
                d = cd[(inode - 1):inode, irotor]
            end

            sigr[inode, irotor] = B[irotor] / (4.0 * pi) * sum(W .* c .* d) / length(d)
        end
    end

    return sigr
end

function calc_reynolds(chord, Wmag_rotor, rho, mu)
    return chord .* abs.(Wmag_rotor) * rho / mu
end

"""
    calculate_gamma_sigma(blade_elements, Vm, vtheta_rotor)

Calculate rotor circulation and source strengths using blade element data and inflow velocities.

# Arguments:
`Vm::Matrix{Float}` : Meridional velocity (including freestream and induced velocity) at rotor planes. Must be a matrix of size number of blade elements by number of rotors.

# Returns:
`Gamr::Matrix{Float}` : Rotor circulations [num blade_elements x num rotors]
`sigr::Matrix{Float}` : Rotor panel source strengths [num blade_elements x num rotors]
"""
function calculate_gamma_sigma(
    blade_elements,
    Wz_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    post=false,
    verbose=false,
)

    # get floating point type
    TF = promote_type(
        eltype(blade_elements[1].chords),
        eltype(blade_elements[1].twists),
        eltype(Wm_rotor),
        eltype(Wtheta_rotor),
    )

    # get problem dimensions (number of radial stations x number of rotors)
    nr, nrotor = size(Wm_rotor)

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr + 1, nrotor)

    # call in-place function
    return calculate_gamma_sigma!(
        Gamr,
        sigr,
        blade_elements,
        Wz_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        post=post,
        verbose=verbose,
    )
end

"""
"""
function calculate_blade_element_coefficients!(
    cl,
    cd,
    beta1,
    alpha,
    reynolds,
    mach,
    blade_elements,
    Wz_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    operating_point;
    post=false,
    verbose=false,
)

    # problem dimensions
    nr, nrotor = size(Wmag_rotor) #num radial stations, num rotors
    # rename for convenience
    be = blade_elements

    if post
        # initialize extra outputs
        phidb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        alphadb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
    end

    # loop through rotors
    for irotor in 1:nrotor

        # loop through radial stations
        for ir in 1:nr

            # extract blade element properties
            # TODO: consider not renaming these to avoid the small allocation
            B = be.B[irotor] # number of blades
            is_stator = be.is_stator[irotor]
            c = be.chords[ir, irotor] # chord length
            stagger = be.stagger[ir, irotor]
            solidity = be.solidity[ir, irotor]

            # calculate angle of attack
            beta1[ir, irotor], alpha[ir, irotor] = calculate_inflow_angles(
                # Wz rather Wm used here in DFDC, presumable because the blade elements are definied in the z-r plane rather than the meridional direction
                # Wm_rotor[ir, irotor], Wtheta_rotor[ir, irotor], stagger
                Wz_rotor[ir, irotor],
                Wtheta_rotor[ir, irotor],
                stagger,
            )

            # get local Reynolds number
            reynolds[ir, irotor] = calc_reynolds(
                c, Wmag_rotor[ir, irotor], operating_point.rhoinf[], operating_point.muinf[]
            )
            # c * abs(Wmag_rotor[ir, irotor]) * operating_point.rhoinf / operating_point.muinf

            # get local Mach number
            mach[ir, irotor] = Wmag_rotor[ir, irotor] / operating_point.asound[]

            # - Update cl and cd values - #
            # be.cl[ir], be.cd[ir] = lookup_clcd(
            cl[ir, irotor], cd[ir, irotor] = lookup_clcd(
                be.inner_airfoil[ir, irotor],
                be.outer_airfoil[ir, irotor],
                be.inner_fraction[ir, irotor],
                c,
                B,
                Wmag_rotor[ir, irotor],
                solidity,
                stagger,
                alpha[ir, irotor],
                beta1[ir, irotor],
                reynolds[ir, irotor],
                mach[ir, irotor],
                operating_point.asound[];
                is_stator=is_stator,
                verbose=verbose,
            )

            if post
                phidb[ir, irotor] = beta1
                alphadb[ir, irotor] = alpha
            end
        end
    end

    if post
        return cl, cd, phidb, alphadb
    else
        return cl, cd
    end
end

function lookup_clcd(
    inner_airfoil,
    outer_airfoil,
    inner_fraction,
    chord,
    B,
    Wmag,
    solidity,
    stagger,
    alpha,
    inflow,
    reynolds,
    mach,
    asound;
    verbose=false,
    is_stator=0,
)

    # look up lift and drag data for the nearest two input sections
    if typeof(inner_airfoil) <: c4b.DFDCairfoil
        # - DFDC Airfoil Parameter - #

        # get inner values
        clin, cdin, _ = c4b.dfdceval(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            inner_airfoil,
            asound;
            verbose=verbose,
            is_stator=is_stator,
        )
        # get outer values
        clout, cdout, _ = c4b.dfdceval(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            outer_airfoil,
            asound;
            verbose=verbose,
            is_stator=is_stator,
        )

    elseif typeof(inner_airfoil) <: c4b.DTCascade
        # - Cascade Lookups - #
        # get inner values
        clin, cdin = c4b.caseval(inner_airfoil, inflow, reynolds, stagger, solidity, mach)
        # get outer values
        clout, cdout = c4b.caseval(outer_airfoil, inflow, reynolds, stagger, solidity, mach)
    elseif typeof(inner_airfoil) <: c4b.AFType
        # - Airfoil Lookups - #
        # get inner values
        clin, cdin = c4b.afeval(inner_airfoil, alpha, reynolds, mach)
        # get outer values
        clout, cdout = c4b.afeval(outer_airfoil, alpha, reynolds, mach)
    elseif typeof(inner_airfoil) <: c4b.ADM
        clin = clfromGamr(inner_airfoil.prescribed_circulation, Wmag, chord)
        clout = clfromGamr(outer_airfoil.prescribed_circulation, Wmag, chord)
        cdin = cdfromsigr(inner_airfoil.prescribed_source_strength, Wmag, chord, B)
        cdout = cdfromsigr(outer_airfoil.prescribed_source_strength, Wmag, chord, B)
    elseif typeof(inner_airfoil) <: c4b.ExternalAirfoil
        # get inner values
        clin, cdin = c4b.external_eval(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            inner_airfoil,
            asound;
            verbose=verbose,
            is_stator=is_stator,
        )

        # get outer values
        clout, cdout = c4b.external_eval(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            outer_airfoil,
            asound;
            verbose=verbose,
            is_stator=is_stator,
        )
    else
        @error "No blade element datatype: $(typeof(inner_airfoil)) defined."
    end

    # interpolate inner and outer values
    cl = FLOWMath.linear([0.0; 1.0], [clin, clout], inner_fraction)
    cd = FLOWMath.linear([0.0; 1.0], [cdin, cdout], inner_fraction)

    return cl, cd
end

function clfromGamr(Gamr, Wmag, c)
    return 2.0 * Gamr / (Wmag * c)
end

function cdfromsigr(sigr, Wmag, c, B)
    return sigr * 4.0 * pi / (B * Wmag * c)
end

"""
"""
function calculate_inflow_angles(Cz_rotor, Ctheta_rotor, stagger)
    #inflow angle
    beta1 = atan.(-Ctheta_rotor, Cz_rotor)

    #angle of attack
    alpha = beta1 - stagger

    return beta1, alpha
end

######################################################################
#                                                                    #
#                           STUFF FOR WAKE                           #
#                                                                    #
######################################################################

"""
    calculate_enthalpy_jumps(Gamr, Omega, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Gamr, Omega, num_blades)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Gamr) .= 0

    calculate_enthalpy_jumps!(H_tilde, Gamr, Omega, num_blades)

    return H_tilde
end

function calculate_enthalpy_jumps!(H_tilde, Gamr, Omega, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)
    # reset H_tilde
    H_tilde .= 0.0

    for ir in 1:nr
        for irotor in 1:nrotor
            omega = Omega[irotor]
            B = num_blades[irotor]

            G = Gamr[ir, irotor]
            H_tilde[ir, irotor] += omega * B * G / (2.0 * pi)

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                H_tilde[ir, irotor] += H_tilde[ir, jrotor]
            end
        end
    end

    return H_tilde
end

"""
    calculate_net_circulation(Gamr, num_blades)

Calculate net circulation from upstream rotors.
"""
function calculate_net_circulation(Gamr, num_blades)

    # calculate net circulations
    Gamma_tilde = similar(Gamr) .= 0

    calculate_net_circulation!(Gamma_tilde, Gamr, num_blades)

    return Gamma_tilde
end

function calculate_net_circulation!(Gamma_tilde, Gamr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)
    # reset Gamma_tilde
    Gamma_tilde .= 0.0

    for ir in 1:nr
        for irotor in 1:nrotor
            B = num_blades[irotor]

            G = Gamr[ir, irotor]
            Gamma_tilde[ir, irotor] += B * G

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                Gamma_tilde[ir, irotor] += Gamma_tilde[ir, jrotor]
            end
        end
    end

    return Gamma_tilde
end
