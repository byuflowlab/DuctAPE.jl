"""
    calculate_induced_velocities_on_rotors(blade_elements, Gamr, vx_rb, vr_rb, gamb,
        vx_rw, vr_rw, gamw)

Calculate axial, radial, and tangential induced velocity on the rotors.

`Gamr::Matrix{Float}` : rotor circulation values (must be a matrix, even if only 1 rotor is being used). Size is number of blade elements by number of rotors.
`vx_rw::Matrix{Matrix{Float}}` : unit axial induced velocities of wakes on rotors. Must be a matrix of size number of rotors by number of wakes made up of matrices of size number of blade elements by number of wake "panels."
`gamw::Matrix{Float}` : wake circulation values (must be a matrix). Size is number of wake sheets by number of "panels" per sheet.
"""
function calculate_induced_velocities_on_rotors(
    blade_elements,
    Gamr,
    vx_rw,
    vr_rw,
    gamw,
    vx_rr,
    vr_rr,
    sigr,
    vx_rb=nothing,
    vr_rb=nothing,
    gamb=nothing;
    debug=false,
)

    # problem dimensions
    _, nrotor = size(Gamr) # number of rotors
    nwake, _ = size(gamw) # number of wake sheets

    # initialize outputs
    vx = similar(Gamr) .= 0 # axial induced velocity
    vr = similar(Gamr) .= 0 # radial induced velocity
    vθ = similar(Gamr) .= 0 # tangential induced velocity

    if debug

        # initialize outputs
        vxb = similar(Gamr) .= 0 # axial induced velocity
        vrb = similar(Gamr) .= 0 # radial induced velocity

        vxw = similar(Gamr) .= 0 # axial induced velocity
        vrw = similar(Gamr) .= 0 # radial induced velocity

        vxr = similar(Gamr) .= 0 # axial induced velocity
        vrr = similar(Gamr) .= 0 # radial induced velocity
    end

    # loop through each affected rotor
    for irotor in 1:nrotor

        # add body induced velocities
        if gamb != nothing
            @views vx[:, irotor] .+= vx_rb[irotor] * gamb
            @views vr[:, irotor] .+= vr_rb[irotor] * gamb

            if debug
                @views vxb[:, irotor] .+= vx_rb[irotor] * gamb
                @views vrb[:, irotor] .+= vr_rb[irotor] * gamb
            end
        end

        # add wake induced velocities
        for jwake in 1:nwake
            @views vx[:, irotor] .+= vx_rw[irotor, jwake] * gamw[jwake, :]
            @views vr[:, irotor] .+= vr_rw[irotor, jwake] * gamw[jwake, :]

            if debug
                @views vxw[:, irotor] .+= vx_rw[irotor, jwake] * gamw[jwake, :]
                @views vrw[:, irotor] .+= vr_rw[irotor, jwake] * gamw[jwake, :]
            end
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vx[:, irotor] .+= vx_rr[irotor, jrotor] * sigr[:, jrotor]
            @views vr[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]

            if debug
                @views vxr[:, irotor] .+= vx_rr[irotor, jrotor] * sigr[:, jrotor]
                @views vrr[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]
            end
        end

        # add self-induced tangential velocity
        B = blade_elements[irotor].B
        r = blade_elements[irotor].rbe
        @views vθ[:, irotor] .+= B .* Gamr[:, irotor] ./ (4 * pi * r)

        # add induced tangential velocity from upstream rotors
        for jrotor in 1:(irotor - 1)
            B = blade_elements[jrotor].B
            r = blade_elements[jrotor].rbe
            @views vθ[:, irotor] .+= B .* Gamr[:, jrotor] ./ (2 * pi * r)
        end
    end

    # return raw induced velocities
    if debug
        return vx, vr, vθ, vxb, vrb, vxw, vrw, vxr, vrr
    else
        return vx, vr, vθ
    end
end

"""
"""
function calculate_rotor_velocities(Gamr, wake_vortex_strengths, sigr, gamb, inputs)

    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        wake_vortex_strengths,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
    )

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ inputs.Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers

    # meridional component
    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # Get the inflow magnitude at the rotor as the combination of all the components
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    return vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor
end

"""
    calculate_gamma_sigma(blade_elements, Vm, Vθ)

Calculate rotor circulation and source strengths using blade element data and inflow velocities.

**Arguments:**
`Vm::Matrix{Float}` : Meridional velocity (including freestream and induced velocity) at rotor planes. Must be a matrix of size number of blade elements by number of rotors.

**Returns:**
`Gamr::Matrix{Float}` : Rotor circulations [num blade_elements x num rotors]
`sigr::Matrix{Float}` : Rotor panel source strengths [num blade_elements x num rotors]
"""
function calculate_gamma_sigma(blade_elements, Wm, Wθ, W; debug=false)

    # get floating point type
    TF = promote_type(
        eltype(blade_elements[1].chords),
        eltype(blade_elements[1].twists),
        eltype(Wm),
        eltype(Wθ),
    )

    # get problem dimensions (number of radial stations x number of rotors)
    nr, nrotor = size(Wm)

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr, nrotor)

    # call in-place function
    return calculate_gamma_sigma!(Gamr, sigr, blade_elements, Wm, Wθ, W; debug=debug)
end

"""
    calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, Vθ)

In-place version of [`calculate_gamma_sigma`](@ref)

Note that circulations and source strengths must be matrices of size number of blade elements by number of rotors. (same dimensions as Vm and Vθ)
"""
function calculate_gamma_sigma!(Gamr, sigr, blade_elements, Wm, Wθ, W; debug=false)

    # problem dimensions
    nr, nrotor = size(Gamr) #num radial stations, num rotors

    if debug
        # initialize extra outputs
        phidb = zeros(eltype(Gamr), size(Gamr))
        alphadb = zeros(eltype(Gamr), size(Gamr))
        cldb = zeros(eltype(Gamr), size(Gamr))
        cddb = zeros(eltype(Gamr), size(Gamr))
    end

    # loop through rotors
    for irotor in 1:nrotor

        # loop through radial stations
        for ir in 1:nr

            # extract blade element properties
            B = blade_elements[irotor].B # number of blades
            c = blade_elements[irotor].chords[ir] # chord length
            twist = blade_elements[irotor].twists[ir] # twist
            r = blade_elements[irotor].rbe[ir] # radius
            Ω = blade_elements[irotor].Omega # rotation rate

            # calculate angle of attack
            phi = atan(Wm[ir, irotor], -Wθ[ir, irotor])
            alpha = twist - phi

            # look up lift and drag data for the nearest two input sections
            # TODO: this breaks rotor aero tests... need to update those.
            clin, cdin = search_polars(blade_elements[irotor].inner_airfoil[ir], alpha)
            clout, cdout = search_polars(blade_elements[irotor].outer_airfoil[ir], alpha)
            # linearly interpolate between those two values at your blade element location
            cl = fm.linear(
                [0.0; 1.0], [clin, clout], blade_elements[irotor].inner_fraction[ir]
            )
            cd = fm.linear(
                [0.0; 1.0], [cdin, cdout], blade_elements[irotor].inner_fraction[ir]
            )

            # calculate vortex strength
            Gamr[ir, irotor] = 1 / 2 * W[ir, irotor] * c * cl

            # calculate source strength
            sigr[ir, irotor] = B / (4 * pi * r) * W[ir, irotor] * c * cd

            if debug
                phidb[ir, irotor] = phi
                alphadb[ir, irotor] = alpha
                cldb[ir, irotor] = cl
                cddb[ir, irotor] = cd
            end
        end
    end

    if debug
        return Gamr, sigr, phidb, alphadb, cldb, cddb
    else
        return Gamr, sigr
    end
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
TODO: add in cascade database search at some point.
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = ccb.afeval(airfoil, alpha, re, ma)
