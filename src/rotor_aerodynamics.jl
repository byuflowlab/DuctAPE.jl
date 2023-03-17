"""
    calculate_induced_velocities_on_rotors(blade_elements, Γr, Ax_br, Ar_br, Γb,
        Ax_wr, Ar_wr, Γw)

Calculate meridional and tangential induced velocity on the rotors.

The meridional velocity is defined tangent to a streamline in the x-r plane.  Its magnitude
can therefore be computed as `Vm = Vx + Vr`  **need to check this**
"""
function calculate_induced_velocities_on_rotors(blade_elements, Γr, Ax_br, Ar_br, Γb,
    Ax_wr, Ar_wr, Γw)

    # problem dimensions
    nr = length(blade_elements)

    # initialize outputs
    Vm = similar(BΓr) .= 0
    Vθ = similar(BΓr) .= 0

    # loop through each rotor
    for i in 1:nr

        # add body induced meridional velocities
        Vm[:, i] .+= Ax_br[i] * Γb
        Vm[:, i] .+= Ar_br[i] * Γb

        # add wake induced meridional velocities
        for iw in 1:size(Γw, 1)
            Vm[:, i] .+= Ax_wr[iw, i] * view(Γw, iw, :)
            Vm[:, i] .+= Ar_wr[iw, i] * view(Γw, iw, :)
        end

        # add rotor induced meridional velocities
        for iw in 1:size(Σr, 1)
            Vm[:, i] .+= Ax_rr[iw, i] * view(Σr, iw, :)
            Vm[:, i] .+= Ar_rr[iw, i] * view(Σr, iw, :)
        end

        # add induced tangential velocity from self
        @views Vθ[:, i] .+= B .* Γr[:, i] ./ (4*pi*blade_elements[i].radial_positions)

        # add induced tangential velocity from previous rotors
        for j = 1:i-1
            @views Vθ[:, i] .+= B .* Γr[:, j] ./ (2*pi*blade_elements[j].radial_positions)
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
        TF = promote_type(eltype(blade_elements[1].chords), typeof(Vinf), eltype(Vm), eltype(Vθ))
    end

    # get problem dimensions
    nbe = length(blade_elements)
    nr = blade_elements[1].num_radial_stations[1]

    # initialize outputs
    Γr = zeros(TF, nbe, nr)
    Σr = zeros(TF, nbe, nr)

    # call in-place function
    return calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf, Vm, Vθ)
end

"""
    calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf [, Vm, Vθ])

In-place version of [`calculate_gamma_sigma`](@ref)
"""
function calculate_gamma_sigma!(Γr, Σr, blade_elements, Vinf, Vm=nothing, Vθ=nothing)

    # problem dimensions
    nbe, nr = size(Γr)

    # loop through rotors
    for i in 1:nbe
        for j = 1:nr

            # extract blade element properties
            B = blade_elements[i].num_blades # number of blades
            c = blade_elements[i].chords[j] # chord length
            twist = blade_elements[i].twists[j] # twist
            r = blade_elements[i].radial_positions[j] # radius
            Ω = blade_elements[i].omega[j] # rotation rate

            # meridional and tangential inflow velocities
            Wm = isnothing(Vm) ? Vinf : Vm[i, j] + Vinf
            Wθ = isnothing(Vθ) ? -Ω*r : Vθ[i, j] - Ω*r
            W = sqrt(Wm^2 + Wθ^2)

            # calculate angle of attack
            alpha = twist - atan(Wm, -Wθ)

            # look up lift and drag data
            cl, cd = search_polars(blade_elements[i].airfoils[j], alpha)

            # calculate vortex strength
            Γr[i, j] = 1/2*W*c*cl

            # calculate source strength
            Σr[i, j] = B/(4*pi*r)*W*c*cd

        end
    end

    return Γr, Σr
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = ccb.afeval(airfoil, alpha, re, ma)
