#=

Functions regarding rotor aerodynamics

Authors: Judd Mehr,

=#

"""
    calculate_induced_velocities_on_rotors(BΓr, Γ_tilde, blade_elements, Ax_br, Ar_br, Γ_b,
        Ax_wr, Ar_wr, Γ_w)

Calculate meridional and tangential induced velocity on the rotors.

The meridional velocity is defined tangent to a streamline in the x-r plane.  Its magnitude
can therefore be computed as `Vm = Vx + Vr`
"""
function calculate_induced_velocities_on_rotors(BΓr, Γ_tilde,
    blade_elements, Ax_br, Ar_br, Γ_b, Ax_wr, Ar_wr, Γ_w)

    # problem dimensions
    nr = length(blade_elements)

    # initialize outputs
    Vm = similar(BΓr) .= 0
    Vθ = similar(BΓr) .= 0

    # loop through each rotor
    for ir in 1:nr

        # add body induced meridional velocities
        Vm[:, ir] .+= Ax_br[ir] * Γ_b
        Vm[:, ir] .+= Ar_br[ir] * Γ_b

        # add wake induced meridional velocities
        for iw in 1:size(Γ_w, 1)
            Vm[:, ir] .+= Ax_wr[iw, ir] * view(Γ_w, iw, :)
            Vm[:, ir] .+= Ar_wr[iw, ir] * view(Γ_w, iw, :)
        end

        # # add rotor induced meridional velocities
        # for iw in 1:size(Σ_r, 1)
        #     Vm[:, ir] .+= Ax_rr[iw, ir] * view(Σ_r, iw, :)
        #     Vm[:, ir] .+= Ar_rr[iw, ir] * view(Σ_r, iw, :)
        # end

        # add self-induced tangential velocity
        Vθ[:, ir] .+= view(BΓr, :, ir) / (4*pi*blade_elements[ir].radial_positions)

        # add tangential velocity from previous rotors
        if ir > 1
            Vθ[:,] .+= view(Γ_tilde, :, ir-1) / (2*pi*blade_elements[ir].radial_positions)
        end

    end

    return Vm, Vθ
end

"""
DONE.HAS TEST. clean up
twist is from plane of rotation
TODO: need to check upstream in the functions that the sign convention for the velocities matches expectations for this function.
"""
function calculate_angle_of_attack(twist, Wm, Wθ)

    inflow =

    # - Calculate Angle of Attack - #

    # if typeof(alpha) == Float64
    #     println("section twist degrees: ", twist * 180 / pi)
    #     println("inflow angle degrees: ", inflow * 180 / pi)
    #     println("angle of attack degrees: ", alpha * 180 / pi)
    # end

    return alpha
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
