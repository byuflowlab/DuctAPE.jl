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
    gamb=nothing,
)

    # problem dimensions
    _, nrotor = size(Gamr) # number of rotors
    _, nwake = size(vx_rw) # number of wake sheets

    # initialize outputs
    vx = similar(Gamr) .= 0 # axial induced velocity
    vr = similar(Gamr) .= 0 # radial induced velocity
    vθ = similar(Gamr) .= 0 # tangential induced velocity

    # loop through each affected rotor
    for irotor in 1:nrotor

        # add body induced velocities
        if gamb != nothing
            @views vx[:, irotor] .+= vx_rb[irotor] * gamb
            @views vr[:, irotor] .+= vr_rb[irotor] * gamb
        end

        # add wake induced velocities
        for jwake in 1:nwake
            @views vx[:, irotor] .+= vx_rw[irotor, jwake] * gamw[jwake, :]
            @views vr[:, irotor] .+= vr_rw[irotor, jwake] * gamw[jwake, :]
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vx[:, irotor] .+= vx_rr[irotor, jrotor] * sigr[:, jrotor]
            @views vr[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]
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
    return vx, vr, vθ
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
function calculate_gamma_sigma(blade_elements, Vm, Vθ)

    # get floating point type
    TF = promote_type(
        eltype(blade_elements[1].chords),
        eltype(blade_elements[1].twists),
        eltype(Vm),
        eltype(Vθ),
    )

    # get problem dimensions (number of radial stations x number of rotors)
    nr, nrotor = size(Vm)

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr, nrotor)

    # call in-place function
    return calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, Vθ)
end

"""
    calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, Vθ)

In-place version of [`calculate_gamma_sigma`](@ref)

Note that circulations and source strengths must be matrices of size number of blade elements by number of rotors. (same dimensions as Vm and Vθ)
"""
function calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, Vθ)

    # problem dimensions
    nr, nrotor = size(Gamr) #num radial stations, num rotors

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

            # meridional and tangential inflow velocities
            # TODO: make sure the velocities getting passed in here are what you expect them to be. track them back.
            Wm = Vm[ir, irotor]
            Wθ = Vθ[ir, irotor] #TODO:!!! Vθ ≂̸ Wθ !!! (could just be a naming issue)
            # TODO: this W isn't quite right. probably want to pass that in as well since sqrt(vx^2 + vr^2 + vt^2) ≂̸ sqrt(sqrt(vx^2+vr^2) + vt^2)
            W = sqrt(Wm^2 + Wθ^2)

            # calculate angle of attack
            alpha = twist - atan(Wm, -Wθ)

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
            Gamr[ir, irotor] = 1 / 2 * W * c * cl

            # calculate source strength
            sigr[ir, irotor] = B / (4 * pi * r) * W * c * cd
        end
    end

    return Gamr, sigr
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
TODO: add in cascade database search at some point.
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = ccb.afeval(airfoil, alpha, re, ma)
