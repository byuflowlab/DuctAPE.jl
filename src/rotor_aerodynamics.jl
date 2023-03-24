"""
    calculate_induced_velocities_on_rotors(blade_elements, Γr, vx_rb, vr_rb, Γb,
        vx_rw, vr_rw, Γw)

Calculate meridional and tangential induced velocity on the rotors.

The meridional velocity is defined tangent to a streamline in the x-r plane.  Its magnitude
can therefore be computed as `Vm = Vx + Vr`  **need to check this**
"""
function calculate_induced_velocities_on_rotors(
    blade_elements, Γr, vx_rb, vr_rb, Γb, vx_rw, vr_rw, Γw, vx_rr, vr_rr, Σr
)

    # problem dimensions
    nr, nrotor = size(Γr)

    # initialize outputs
    Vm = similar(Γr) .= 0
    Vθ = similar(Γr) .= 0

    # loop through each rotor
    for jrotor in 1:nrotor

        # add body induced meridional velocities
        @views Vm[:, jrotor] .+= vx_rb[jrotor] * Γb
        @views Vm[:, jrotor] .+= vr_rb[jrotor] * Γb

        # add wake induced meridional velocities
        for iwake in 1:nwake
            @views Vm[:, jrotor] .+= vx_rw[iwake, jrotor] * view(Γw, :, iwake)
            @views Vm[:, jrotor] .+= vr_rw[iwake, jrotor] * view(Γw, :, iwake)
        end

        # add rotor induced meridional velocities
        for irotor in 1:nrotor
            @views Vm[:, jrotor] .+= vx_rr[irotor, jrotor] * view(Σr, :, irotor)
            @views Vm[:, jrotor] .+= vr_rr[irotor, jrotor] * view(Σr, :, irotor)
        end

        # add induced tangential velocity from self
        B = blade_elements[jrotor].num_blades
        r = blade_elements[jrotor].radial_positions
        @views Vθ[:, jrotor] .+= B .* Γr[:, jrotor] ./ (4 * pi * r)

        # add induced tangential velocity from previous rotors
        for krotor in 1:(jrotor - 1)
            B = blade_elements[krotor].num_blades
            r = blade_elements[krotor].radial_positions
            @views Vθ[:, jrotor] .+= B .* Γr[:, krotor] ./ (2 * pi * r)
        end
    end

    return Vm, Vθ
end

"""
    calculate_gamma_sigma(blade_elements, Vinf [, Vm, Vθ])

Calculate rotor circulation and source strengths using airfoil data.
"""
function calculate_gamma_sigma(blade_elements, Vinf, Vm=nothing, Vθ=nothing)

    # get floating point type
    if isnothing(Vm) && isnothing(Vθ)
        TF = promote_type(eltype(blade_elements[1].chords), typeof(Vinf))
    else
        TF = promote_type(
            eltype(blade_elements[1].chords), typeof(Vinf), eltype(Vm), eltype(Vθ)
        )
    end

    # get problem dimensions
    nr = blade_elements[1].num_radial_stations[1]
    nrotor = length(blade_elements)

    # initialize outputs
    Γr = zeros(TF, nr, nrotor)
    Σr = zeros(TF, nr, nrotor)

    # call in-place function
    return calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf, Vm, Vθ)
end

"""
    calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf [, Vm, Vθ])

In-place version of [`calculate_gamma_sigma`](@ref)
"""
function calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf, Vm=nothing, Vθ=nothing)

    # problem dimensions
    nr, nrotor = size(Γr)

    # loop through rotors
    for ir in 1:nr
        for irotor in 1:nrotor

            # extract blade element properties
            B = blade_elements[ir].num_blades # number of blades
            c = blade_elements[ir].chords[irotor] # chord length
            twist = blade_elements[ir].twists[irotor] # twist
            r = blade_elements[ir].radial_positions[irotor] # radius
            Ω = blade_elements[ir].omega[irotor] # rotation rate

            # meridional and tangential inflow velocities
            Wm = isnothing(Vm) ? Vinf : Vm[ir, irotor] + Vinf
            Wθ = isnothing(Vθ) ? -Ω * r : Vθ[ir, irotor] - Ω * r
            W = sqrt(Wm^2 + Wθ^2)

            # calculate angle of attack
            alpha = twist - atan(Wm, -Wθ)

            # look up lift and drag data
            cl, cd = search_polars(blade_elements[ir].airfoils[irotor], alpha)

            # calculate vortex strength
            Γr[ir, irotor] = 1 / 2 * W * c * cl

            # calculate source strength
            Σr[ir, irotor] = B / (4 * pi * r) * W * c * cd
        end
    end

    return Γr, Σr
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = ccb.afeval(airfoil, alpha, re, ma)
