"""
    calculate_induced_velocities_on_rotors(blade_elements, Gamr, vz_rb, vr_rb, gamb,
        vz_rw, vr_rw, gamw)

Calculate axial, radial, and tangential induced velocity on the rotors.

`Gamr::Matrix{Float}` : rotor circulation values (must be a matrix, even if only 1 rotor is being used). Size is number of blade elements by number of rotors.
`vz_rw::Matrix{Matrix{Float}}` : unit axial induced velocities of wakes on rotors. Must be a matrix of size number of rotors by number of wakes made up of matrices of size number of blade elements by number of wake "panels."
`gamw::Matrix{Float}` : wake circulation values (must be a matrix). Size is number of wake sheets by number of "panels" per sheet.
"""
function calculate_induced_velocities_on_rotors(
    blade_elements,
    Gamr,
    vz_rw,
    vr_rw,
    gamw,
    vz_rr,
    vr_rr,
    sigr,
    vz_rb=nothing,
    vr_rb=nothing,
    gamb=nothing,
    vz_rbte=nothing,
    vr_rbte=nothing,
    TEidxs=nothing;
    debug=false,
)

    # problem dimensions
    _, nrotor = size(Gamr) # number of rotors

    # initialize outputs
    vz_rotor = similar(Gamr) .= 0 # axial induced velocity
    vr_rotor = similar(Gamr) .= 0 # radial induced velocity
    vtheta_rotor = similar(Gamr) .= 0 # tangential induced velocity

    if debug

        # initialize outputs
        vzb_rotor = similar(Gamr) .= 0 # axial induced velocity
        vrb_rotor = similar(Gamr) .= 0 # radial induced velocity

        vzw_rotor = similar(Gamr) .= 0 # axial induced velocity
        vrw_rotor = similar(Gamr) .= 0 # radial induced velocity

        vzr_rotor = similar(Gamr) .= 0 # axial induced velocity
        vrr_rotor = similar(Gamr) .= 0 # radial induced velocity
    end

    # loop through each affected rotor
    for irotor in 1:nrotor

        # add body induced velocities
        if gamb != nothing
            @views vz_rotor[:, irotor] .+= vz_rb[irotor] * gamb
            @views vr_rotor[:, irotor] .+= vr_rb[irotor] * gamb

            if debug
                @views vzb_rotor[:, irotor] .+= vz_rb[irotor] * gamb
                @views vrb_rotor[:, irotor] .+= vr_rb[irotor] * gamb
            end

            #note: v?_rbte[irotor] should have been defined as zeros for endpoints that are not coincident.
            @views vz_rotor .+= vz_rbte[irotor] * gamb[TEidxs]
            @views vr_rotor .+= vr_rbte[irotor] * gamb[TEidxs]

            if debug
                @views vzb_rotor .+= vz_rbte[irotor] * gamb[TEidxs]
                @views vrb_rotor .+= vr_rbte[irotor] * gamb[TEidxs]
            end
        end

        # add wake induced velocities
        @views vz_rotor[:, irotor] .+= vz_rw[irotor] * gamw
        @views vr_rotor[:, irotor] .+= vr_rw[irotor] * gamw

        if debug
            @views vzw_rotor[:, irotor] .+= vz_rw[irotor] * gamw[:]
            @views vrw_rotor[:, irotor] .+= vr_rw[irotor] * gamw[:]
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vz_rotor[:, irotor] .+= vz_rr[irotor, jrotor] * sigr[:, jrotor]
            @views vr_rotor[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]

            if debug
                @views vzr_rotor[:, irotor] .+= vz_rr[irotor, jrotor] * sigr[:, jrotor]
                @views vrr_rotor[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]
            end
        end

        # add self-induced tangential velocity
        B = blade_elements[irotor].B
        r = blade_elements[irotor].rbe
        @views vtheta_rotor[:, irotor] .+= B .* Gamr[:, irotor] ./ (4 * pi * r)

        # add induced tangential velocity from upstream rotors
        for jrotor in 1:(irotor - 1)
            B = blade_elements[jrotor].B
            r = blade_elements[jrotor].rbe
            @views vtheta_rotor[:, irotor] .+= B .* Gamr[:, jrotor] ./ (2 * pi * r)
        end
    end

    # return raw induced velocities
    if debug
        return vz_rotor,
        vr_rotor, vtheta_rotor, vzb_rotor, vrb_rotor, vzw_rotor, vrw_rotor, vzr_rotor,
        vrr_rotor
    else
        return vz_rotor, vr_rotor, vtheta_rotor
    end
end

"""
"""
function reframe_rotor_velocities(
    vz_rotor, vr_rotor, vtheta_rotor, Vinf, Omega, rotor_panel_centers
)

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wz_rotor = vz_rotor .+ Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = similar(vtheta_rotor) .= vtheta_rotor

    for i in 1:length(Omega)
        Wtheta_rotor[:, i] .-= Omega[i] .* rotor_panel_centers[:, i]
    end

    # meridional component
    Wm_rotor = sqrt.(Wz_rotor .^ 2 .+ vr_rotor .^ 2)

    # Get the inflow magnitude at the rotor as the combination of all the components
    Wmag_rotor = sqrt.(Wz_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    return Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor
end

"""
"""
function calculate_rotor_velocities(Gamr, gamw, sigr, gamb, inputs)

    # - Get induced velocities on rotor planes - #
    vz_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vz_rw,
        inputs.vr_rw,
        gamw,
        inputs.vz_rr,
        inputs.vr_rr,
        sigr,
        inputs.vz_rb,
        inputs.vr_rb,
        gamb,
        inputs.vz_rbte,
        inputs.vr_rbte,
        (p -> p.idx).(inputs.body_vortex_panels.TEnodes),
    )

    # - Reframe rotor velocities into blade element frames
    Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vz_rotor,
        vr_rotor,
        vtheta_rotor,
        inputs.freestream.Vinf,
        inputs.blade_elements.Omega,
        inputs.rotor_panel_centers,
    )

    return vz_rotor, vr_rotor, vtheta_rotor, Wz_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor
end


"""
"""
function gamma_from_coeffs!(Gamr, Wmag_rotor, blade_elements, cl)

    for (irotor, be) in enumerate(blade_elements)
    # calculate vortex strength
    Gamr[:,irotor] = @. 1.0 / 2.0 * Wmag_rotor[:,irotor] * be.chords * cl[:,irotor]
end

    return Gamr
end

# function sigma_from_coeffs!(sigr, Wmag_rotor, blade_elements, rho)
function sigma_from_coeffs!(sigr, Wmag_rotor, blade_elements,cd, rho)

    # calculate source strength
    # TODO: need to figure out if there should be a division by r here or not.
    #NOTE: theory sigr doesn't make any sense. consider using actual drag equation until you have good reason not to
    #= it doesn't look like DFDC has the division, even though it's in the theory doc
    #it's possible (but I haven't been able to find it, that they divide by r later for some reason
    =#
    # this one is in the theory doc
    # things converge twice as fast with this one...
    # sigr[:] .= B / (4.0 * pi * r) * Wmag_rotor * c * cd
    # this one is in the DFDC source code
    # things converge slower, likely because the values are so small...
    # sigr[:] .= B / (4.0 * pi) * Wmag_rotor * c * cd

    for (irotor, be) in enumerate(blade_elements)
        B = be.B

        for inode in 1:size(sigr, 1)
            if inode == 1
                W = Wmag_rotor[inode, irotor]
                c = be.chords[inode]
                # cd = be.cd[inode]
                d = cd[inode]

            elseif inode == size(sigr, 1)
                W = Wmag_rotor[inode - 1, irotor]
                c = be.chords[inode - 1]
                # cd = be.cd[inode-1]
                d = cd[inode - 1]

            else
                W = Wmag_rotor[(inode - 1):inode, irotor]
                c = be.chords[(inode - 1):inode]
                # cd = be.cd[inode-1:inode]
                d = cd[(inode - 1):inode]
            end

            sigr[inode, irotor] = B / (4 * pi) * sum(W .* c .* d) / length(d)
            # sigr[:] .= B / (4 * pi) * sum(W .* c .* cd ./ r) / length(cd)
            # sigr[:] .= B / (4 * pi) * rho * sum(W^2 .* c .* cd ./ r) / length(cd)
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

**Arguments:**
`Vm::Matrix{Float}` : Meridional velocity (including freestream and induced velocity) at rotor planes. Must be a matrix of size number of blade elements by number of rotors.

**Returns:**
`Gamr::Matrix{Float}` : Rotor circulations [num blade_elements x num rotors]
`sigr::Matrix{Float}` : Rotor panel source strengths [num blade_elements x num rotors]
"""
function calculate_gamma_sigma(
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
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
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        debug=debug,
        verbose=verbose,
    )
end

"""
    calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, vtheta_rotor)

In-place version of [`calculate_gamma_sigma`](@ref)

Note that circulations and source strengths must be matrices of size number of blade elements by number of rotors. (same dimensions as Vm and vtheta_rotor)
"""
function calculate_gamma_sigma!(
    Gamr,
    sigr,
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
    verbose=false,
)

    #update coeffs
    if debug
        cl, cd, phidb, alphadb, cldb, cddb = update_coeffs!(
            blade_elements,
            Wm_rotor,
            Wtheta_rotor,
            Wmag_rotor,
            freestream;
            debug=debug,
            verbose=verbose,
        )
    else
        cl, cd = update_coeffs!(
            blade_elements,
            Wm_rotor,
            Wtheta_rotor,
            Wmag_rotor,
            freestream;
            debug=debug,
            verbose=verbose,
        )
    end

    # - get circulation strength - #
    gamma_from_coeffs!(Gamr, Wmag_rotor, blade_elements, cl)

    # - get source strength - #
    sigma_from_coeffs!(sigr, Wmag_rotor, blade_elements, cd, freestream.rhoinf)

    return Gamr, sigr
end

function update_coeffs!(
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
    verbose=false,
)

    # problem dimensions
    nr, nrotor = size(Wmag_rotor) #num radial stations, num rotors

    cl = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
    cd = zeros(eltype(Wmag_rotor), size(Wmag_rotor))

    if debug
        # initialize extra outputs
        phidb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        alphadb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        cldb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        cddb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
    end

    # loop through rotors
    for (irotor, be) in enumerate(blade_elements)

        # loop through radial stations
        for ir in 1:nr

            # extract blade element properties
            B = be.B # number of blades
            omega = be.Omega # rotation rate
            fliplift = be.fliplift

            c = be.chords[ir] # chord length
            twist = be.twists[ir] # twist
            stagger = be.stagger[ir]
            r = be.rbe[ir] # radius
            solidity = be.solidity[ir]

            # calculate angle of attack
            phi, alpha = calculate_inflow_angles(
                Wm_rotor[ir, irotor], Wtheta_rotor[ir, irotor], stagger
            )

            # get local Reynolds number
            reynolds = calc_reynolds(
                c, Wmag_rotor[ir, irotor], freestream.rhoinf, freestream.muinf
            )
            # c * abs(Wmag_rotor[ir, irotor]) * freestream.rhoinf / freestream.muinf

            # get local Mach number
            mach = Wmag_rotor[ir, irotor] / freestream.asound

            # - Update cl and cd values - #
            # be.cl[ir], be.cd[ir] = lookup_clcd(
            cl[ir, irotor], cd[ir, irotor] = lookup_clcd(
                blade_elements[irotor].inner_airfoil[ir],
                blade_elements[irotor].outer_airfoil[ir],
                blade_elements[irotor].inner_fraction[ir],
                Wmag_rotor[ir, irotor],
                solidity,
                stagger,
                alpha,
                phi,
                reynolds,
                mach,
                freestream.asound;
                verbose=verbose,
                fliplift=fliplift,
            )

            if debug
                phidb[ir, irotor] = phi
                alphadb[ir, irotor] = alpha
                cldb[ir, irotor] = cl
                cddb[ir, irotor] = cd
            end
        end
    end

    if debug
        return cl, cd, phidb, alphadb, cldb, cddb
        # return phidb, alphadb, cldb, cddb
    else
        # return nothing
        return cl, cd
    end
end

function lookup_clcd(
    inner_airfoil,
    outer_airfoil,
    inner_fraction,
    Wmag,
    solidity,
    stagger,
    alpha,
    inflow,
    reynolds,
    mach,
    asound;
    verbose=false,
    fliplift=false,
)

    # look up lift and drag data for the nearest two input sections
    # TODO: this breaks rotor aero tests... need to update those.
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
            fliplift=fliplift,
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
            fliplift=fliplift,
        )

    elseif typeof(inner_airfoil) <: c4b.DTCascade
        # - Cascade Lookups - #
        # get inner values
        clin, cdin = c4b.caseval(inner_airfoil, stagger, inflow, reynolds, mach, solidity)
        # get outer values
        clout, cdout = c4b.caseval(outer_airfoil, stagger, inflow, reynolds, mach, solidity)
    elseif typeof(inner_airfoil) <: c4b.AFType
        # - Airfoil Lookups - #
        # get inner values
        clin, cdin = search_polars(inner_airfoil, alpha)
        # get outer values
        clout, cdout = search_polars(outer_airfoil, alpha)
    else
        @error "No blade element datatype: $(typeof(inner_airfoil)) defined."
    end

    # interpolate inner and outer values
    cl = fm.linear([0.0; 1.0], [clin, clout], inner_fraction)
    cd = fm.linear([0.0; 1.0], [cdin, cdout], inner_fraction)

    return cl, cd
end

"""
"""
function calculate_inflow_angles(Wm_rotor, Wtheta_rotor, stagger)
    #inflow angle
    # phi = atan.(Wm_rotor, -Wtheta_rotor)
    phi = atan.(-Wtheta_rotor, Wm_rotor)

    #angle of attack
    alpha = phi - stagger
    # alpha = twist .- phi

    return phi, alpha
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
TODO: add in cascade database search at some point.
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = c4b.afeval(airfoil, alpha, re, ma)

######################################################################
#                                                                    #
#                           STUFF FOR WAKE                           #
#                                                                    #
######################################################################

"""
only used in post-process for cp.
expression not in dfdc theory, comes from source code.
todo: derive in theory doc
"""
function calculate_entropy_jumps(sigr, Vm)
    return sigr .* Vm
end

"""
    calculate_enthalpy_jumps(Gamr, Omega, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Gamr, Omega, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Gamr) .= 0
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

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # calculate net circulations
    Gam_tilde = similar(Gamr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            B = num_blades[irotor]

            G = Gamr[ir, irotor]
            Gam_tilde[ir, irotor] += B * G

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                Gam_tilde[ir, irotor] += Gam_tilde[ir, jrotor]
            end
        end
    end

    return Gam_tilde
end
