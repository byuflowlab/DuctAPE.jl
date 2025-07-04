"""
    struct Imperial

Disptach type used for setting Operating Point Units to Imperial

# Units: feet, pounds, seconds
- Velocity: feet per second
- Temperature: Fahrenheit
- Pressure: pounds per square foot
- Density: slugs per cubic foot
- Dynamic Viscosity: slugs per square foot
"""
struct Imperial end

"""
    struct SI

Disptach type used for setting Operating Point Units to SI.
Note that the default is SI, and is really only used in the backend.

# Units: meters, kilograms, seconds
- Velocity: meters per second
- Temperature: Kelvin
- Pressure: Pascals
- Density: kilograms per meter cubed
- Dynamic Viscosity: Pascal seconds
"""
struct SI end

"""
    OperatingPoint(Vinf, Minf, rhoinf, muinf, asound, Ptot, Ttot, Omega)
    OperatingPoint(
        Vinf, Omega, rhoinf=nothing, muinf=nothing, asound=nothing; altitude=0.0
    )
    OperatingPoint(
        ::Imperial, Vinf, Omega, rhoinf=nothing, muinf=nothing, asound=nothing; altitude=0.0
    )

DuctedRotor operating point information.

Functions that take in `altitude` will populate undefined thermodynamic properties of the freestream using a standard atmosphere model, ideal gas law, and Sutherland's law; defaulting to SI units.
If the `::Imperial` dispatch type is input, then the thermodynamic properties will be converted to Imperial units.

# Fields/Arguments

- `Vinf::AbstractVector{Float}` : Freestream velocity magnitude (which is only in the axial direction).
- `Minf::AbstractVector{Float}` : Freestream Mach number
- `rhoinf::AbstractVector{Float}` : Freestream density
- `muinf::AbstractVector{Float}` : Freestream viscosity
- `asound::AbstractVector{Float}` : Freestream speed of sound
- `Ptot::AbstractVector{Float}` : Freestream total pressure
- `Ttot::AbstractVector{Float}` : Freestream total temperature
- `Omega::AbstractVector{Float}` : Rotor rototation rate(s)
"""
struct OperatingPoint{
    Ta<:AbstractVector,
    TM<:AbstractVector,
    Tm<:AbstractVector,
    To<:AbstractVector,
    Tp<:AbstractVector,
    Tr<:AbstractVector,
    Tt<:AbstractVector,
    Tv<:AbstractVector,
    Tu,
}
    units::Tu
    Vinf::Tv
    Minf::TM
    rhoinf::Tr
    muinf::Tm
    asound::Ta
    Ptot::Tp
    Ttot::Tt
    Omega::To
end

function OperatingPoint(
    units::Imperial,
    Vinf,
    Omega,
    rhoinf=nothing,
    muinf=nothing,
    asound=nothing,
    altitude=0.0,
)

    # Get thermodynamic properties
    Tinf, Pinf, rho_inf, mu_inf = standard_atmosphere(units, altitude)

    # freestream density
    if isnothing(rhoinf)
        rhoinf = rho_inf
    end

    # freestream dynamic viscosity
    if isnothing(muinf)
        muinf = mu_inf
    end

    # freestream speed of sound
    if isnothing(asound)
        asound = speed_of_sound(Pinf, rhoinf)
    end

    # freestream Mach
    Minf = calculate_mach.(Vinf, asound)

    # freestream total pressure and temperature
    Ptot = total_pressure.(Pinf, Minf)
    Ttot = total_temperature.(Tinf, Minf)

    return OperatingPoint(
        Imperial(),
        isscalar(Vinf) ? [Vinf] : Vinf,
        isscalar(Minf) ? [Minf] : Minf,
        isscalar(rhoinf) ? [rhoinf] : rhoinf,
        isscalar(muinf) ? [muinf] : muinf,
        isscalar(asound) ? [asound] : asound,
        isscalar(Ptot) ? [Ptot] : Ptot,
        isscalar(Ttot) ? [Ttot] : Ttot,
        isscalar(Omega) ? [Omega] : Omega,
    )
end

function OperatingPoint(
    Vinf, Omega, rhoinf=nothing, muinf=nothing, asound=nothing; altitude=0.0
)

    # Get thermodynamic properties
    Tinf, Pinf, rho_inf, mu_inf = standard_atmosphere(altitude)

    # freestream density
    if isnothing(rhoinf)
        rhoinf = rho_inf
    end

    # freestream dynamic viscosity
    if isnothing(muinf)
        muinf = mu_inf
    end

    # freestream speed of sound
    if isnothing(asound)
        asound = speed_of_sound(Pinf, rhoinf)
    end

    # freestream Mach
    Minf = calculate_mach.(Vinf, asound)

    # freestream total pressure and temperature
    Ptot = total_pressure.(Pinf, Minf)
    Ttot = total_temperature.(Tinf, Minf)

    return OperatingPoint(
        SI(),
        isscalar(Vinf) ? [Vinf] : Vinf,
        isscalar(Minf) ? [Minf] : Minf,
        isscalar(rhoinf) ? [rhoinf] : rhoinf,
        isscalar(muinf) ? [muinf] : muinf,
        isscalar(asound) ? [asound] : asound,
        isscalar(Ptot) ? [Ptot] : Ptot,
        isscalar(Ttot) ? [Ttot] : Ttot,
        isscalar(Omega) ? [Omega] : Omega,
    )
end

"""
    ReferenceParameters(Vref, Rref)

Reference parameters for post-process non-dimensionalization.

Note that the actual struct requires the inputs to be arrays, but there is a constructor function that will take in scalars and automatically build the array-based struct.

# Arguments

- `Vref::AbstractVector{Float}` : Reference velocity.
- `Rref::AbstractVector{Float}` : Reference rotor tip radius.
"""
struct ReferenceParameters{Tv<:AbstractVector,Tr<:AbstractVector}
    Vref::Tv
    Rref::Tr
end

function ReferenceParameters(Vref, Rref)
    return ReferenceParameters(
        isscalar(Vref) ? [Vref] : Vref, isscalar(Rref) ? [Rref] : Rref
    )
end

"""
    PanelingConstants(
        num_duct_inlet_panels,
        num_center_body_inlet_panels,
        num_panels,
        dte_minus_cbte,
        num_wake_sheets,
        wake_length=1.0,
    )

Constants used in re-paneling geometry.

Note that unlike other input structures, this one, in general, does not define fields as vectors.  This is because these values should not change throughout an optimization, even if the geometry may change.  Otherwise, discontinuities could be experienced.

# Arguments

- `num_duct_inlet_panels::Int` : The number of panels to use for the casing inlet.
- `num_center_body_inlet_panels::Int` : The number of panels to use for the center body inlet.
- `num_panels::AbstractVector{Int}` : A vector containing the number of panels between discrete locations inside the wake. Specifically, the number of panels between the rotors, between the last rotor and the first body trailing edge, between the body trailing edges (if different), and between the last body trailing edge and the end of the wake.  The length of this vector should be N+1 (where N is the number of rotors) if the duct and center_body trailing edges are aligned, and N+2 if not.
- `dte_minus_cbte::Float` : An indicator concerning the hub and duct trailing edge relative locations. Should be set to -1 if the duct trailing edge axial position minus the center_body trailing edge axial position is negative, +1 if positive (though any positive or negative number will suffice), and zero if the trailing edges are aligned.
- `num_wake_sheets::Int` : The number of wake sheets to use. Note this will also be setting the number of blade elements to use.
- `wake_length::Float=1.0` : Non-dimensional (based on the length from the foremost body leading edge and the aftmost body trailing edge) length of the wake extending behind the aftmost body trailing edge.
"""
struct PanelingConstants{TI,TF,TFI}
    num_duct_inlet_panels::TI
    num_center_body_inlet_panels::TI
    num_panels::AbstractVector{TI}
    dte_minus_cbte::TFI
    num_wake_sheets::TI
    wake_length::TF
end

function PanelingConstants(
    num_duct_inlet_panels,
    num_center_body_inlet_panels,
    num_panels,
    dte_minus_cbte,
    num_wake_sheets,
    wake_length=1.0,
)

    # initialize error messages
    throw_error = false
    error_messages = ""
    error_count = 1

    if num_duct_inlet_panels <= 0
        throw_error = true
        error_messages *= "\n\tError $(error_count): `num_duct_inlet_panels` cannot be fewer than 1; must have a non-zero, positive number of panels for the duct inlet."
        error_count += 1
    end
    if num_center_body_inlet_panels <= 0
        throw_error = true
        error_messages *= "\n\tError $(error_count): `num_center_body_inlet_panels` cannot be fewer than 1; must have a non-zero, positive number of panels for the center body inlet"
        error_count += 1
    end
    if any(num_panels .<= 0)
        throw_error = true
        error_messages *= "\n\tError $(error_count): at least one entry in `num_panels` is less than 1; must have non-zero, positive numbers of panels."
        error_count += 1
    end
    if num_wake_sheets < 3
        throw_error = true
        error_messages *= "\n\tError $(error_count): `num_wake_sheets` must be at least 3."
        error_count += 1
    end
    if wake_length < 0
        throw_error = true
        error_messages *= "\n\tError $(error_count): cannont have a negative `wake_length`."
        error_count += 1
    end

    if throw_error
        throw(error_messages)
        return nothing
    else
        return PanelingConstants(
            num_duct_inlet_panels,
            num_center_body_inlet_panels,
            num_panels,
            dte_minus_cbte,
            num_wake_sheets,
            wake_length,
        )
    end
end

"""
    Rotor(
        B, rotor_axial_position, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, is_stator
    )

Composite type containing the rotor(s) geometric properties.

Note that the actual struct requires the inputs to be arrays, but there is a constructor function that will take in scalars and automatically build the array-based struct.

# Arguments
- `B::AbstractVector{Float}` : The number of blades for each rotor. May not be an integer, but usually is.
- `rotor_axial_position::AbstractVector{Float}` : Dimensional, axial position of each rotor.
- `r::AbstractArray{Float}` : Non-dimensional radial locations of each blade element.
- `Rhub::AbstractVector{Float}` : Dimensional hub radius of rotor. (may be changed if it does not match the radial position of the `center_body` geometry at the selected `rotor_axial_position`.
- `Rtip::AbstractVector{Float}` : Dimensional tip radius of rotor. Is used to determine the radial position of the duct if the `autoshiftduct` option is selected.
- `chords::AbstractArray{Float}` : Dimensional chord lengths of the blade elements.
- `twists::AbstractArray{Float}` : Blade element angles, in radians.
- `tip_gap::AbstractVector{Float}` : Currently unused, do not set to anything other than zeros.
- `airfoils::Vector{Vector{AFType}}` : Airfoil types describing the airfoil polars for each rotor and blade element [[rotor 1 airfoils], [rotor 2 airfoils], ...].
- `is_stator::AbstractVector{Bool}` : Flag to indicate if the airfoil lift values should be flipped or not.

# Keyword Arguments
- `i_know_what_im_doing::Bool=false` : if set to false, checks the twist angles, and if greater than 1.75, prints a warning and converts from degrees to radians.
"""
struct Rotor{
    Tb<:AbstractVector,
    TRz<:AbstractVector,
    Tr<:AbstractArray,
    TRh<:AbstractVector,
    TRt<:AbstractVector,
    Tc<:AbstractArray,
    Tt<:AbstractArray,
    TTg<:AbstractVector,
    Taf<:AbstractArray,
    Tf<:AbstractVector,
}
    B::Tb
    rotor_axial_position::TRz
    r::Tr
    Rhub::TRh
    Rtip::TRt
    chords::Tc
    twists::Tt
    tip_gap::TTg
    airfoils::Taf
    is_stator::Tf
end

function Rotor(
    B,
    rotor_axial_position,
    r,
    Rhub,
    Rtip,
    chords,
    twists,
    tip_gap,
    airfoils,
    is_stator,
    i_know_what_im_doing=false,
)

    # - Format Inputs - #
    # check twists
    if !i_know_what_im_doing && length(findall(t -> t > 1.75, twists)) > 2
        @warn "It looks like your input twist angles may be in degrees. Note that the required units for twist are radians. Converting to radians for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic conversion)."
        twists .*= pi / 180.0
    end

    # - Check Errors - #

    # initialize error messages
    throw_error = false
    error_messages = ""
    error_count = 1

    # check that rotors aren't coincident
    if length(unique(rotor_axial_position)) < length(rotor_axial_position)
        throw_error = true
        error_messages *= "\n\tError $(error_count): Cannot place rotors on top of eachother."
        error_count += 1
    end

    # check that tip is larger than hub
    if Rhub > Rtip
        throw_error = true
        error_messages *= "\n\tError $(error_count): `Rtip` must be greater than `Rhub`."
        error_count += 1
    end

    # check that r's defined from hub to tip
    if !all(r[2:end] .> r[1:(end - 1)])
        throw_error = true
        error_messages *= "\n\tError $(error_count): Radial positions, `r`, must be increasing across the blade"
        error_count += 1
    end

    # check that r's aren't negative
    if any(r .<= 0)
        throw_error = true
        error_messages *= "\n\tError $(error_count): Radial positions, `r`, must be postive, non-zero"
        error_count += 1
    end

    if throw_error
        throw(error_messages)
        return nothing
    else
        return Rotor(
            isscalar(B) ? [B] : B,
            if isscalar(rotor_axial_position)
                [rotor_axial_position]
            else
                rotor_axial_position
            end,
            isscalar(r) ? [r] : r,
            isscalar(Rhub) ? [Rhub] : Rhub,
            isscalar(Rtip) ? [Rtip] : Rtip,
            isscalar(chords) ? [chords] : chords,
            isscalar(twists) ? [twists] : twists,
            isscalar(tip_gap) ? [tip_gap] : tip_gap,
            if typeof(airfoils) <: Union{c4b.AFType,c4b.DTCascade,c4b.DFDCairfoil,c4b.ADM}
                [airfoils]
            else
                airfoils
            end,
            isscalar(is_stator) ? [is_stator] : is_stator,
        )
    end
end

"""
    DuctedRotor(
        duct_coordinates, center_body_coordinates, rotor, paneling_constants; kwargs...
    )

# Arguments

- `duct_coordinates::AbstractMatrix` : The [z r] coordinates of the duct geometry beginning at the inner (casing) side trailing edge and proceeding clockwise. Note that the duct geometry absolute radial position does not need to be included here if the `autoshiftduct` option is selected.
- `center_body_coordinates::AbstractMatrix` : The [z r] coordinates of the center_body beginning at the leading edge and ending at the trailing edge. Note that the leading edge is assumed to be placed at a radial distance of 0.0 from the axis of rotation.
- `rotor::Rotor` : Rotor (and possibly stator) geometric paramters.
- `paneling_constants::PanelingConstants` : Constants used in re-paneling the geometry.

# Keyword Arguments
- `i_know_what_im_doing::Bool=false` : if set to false, performs various checks on the inputs and manually adjusts in some cases with a warning, and errors in cases that cannot be automatically adjusted.
- `silence_warnings::Bool=false` : if set to false, prints warnings when automatically adjusting inputs.
"""
struct DuctedRotor{Td<:AbstractMatrix,Tcb<:AbstractMatrix,Trp<:Rotor,Tpc<:PanelingConstants}
    duct_coordinates::Td
    center_body_coordinates::Tcb
    rotor::Trp
    paneling_constants::Tpc
end

function DuctedRotor(
    duct_coordinates,
    center_body_coordinates,
    rotor,
    paneling_constants;
    i_know_what_im_doing=false,
    silence_warnings=false,
)

    ### --- FORMAT INPUTS --- ###

    # check shape of duct coordinates
    if i_know_what_im_doing && size(duct_coordinates, 1) < size(duct_coordinates, 2)
        if !silence_warnings
            @warn "It appears that the duct coordinates have been provided as rows rather than columns. Permuting dimensions for you (this is an allocating action...). Set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing."
        end
        duct_coordinates = permutedims(duct_coordinates)
    end

    # check clockwise duct coordinates
    if i_know_what_im_doing && duct_coordinates[2, 2] > duct_coordinates[end - 1, 2]
        if !silence_warnings
            @warn "It appears that the duct coordinates have been provided counter_clockwise rather than clockwise. Reversing direction for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing)."
        end
        reverse!(duct_coordinates; dims=1)
    end

    # check z and r duct coordinates
    if i_know_what_im_doing &&
        abs(-(extrema(duct_coordinates[:, 2])...)) >
       abs(-(extrema(duct_coordinates[:, 1])...))
        if !silence_warnings
            @warn "It appears that the duct coordinates have been provided as (r,z) rather than (z,r). Reversing direction for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing)."
        end
        reverse!(duct_coordinates; dims=2)
    end

    # check shape of  center body coordinates
    if i_know_what_im_doing &&
        size(center_body_coordinates, 1) < size(center_body_coordinates, 2)
        if !silence_warnings
            @warn "It appears that the center body coordinates have been provided as rows rather than columns. Permuting dimensions for you (this is an allocating action...). Set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing."
        end
        duct_coordinates = permutedims(center_body_coordinates)
    end

    # check clockwise center body coordinates
    if i_know_what_im_doing &&
        center_body_coordinates[1, 1] > center_body_coordinates[end, 1]
        if !silence_warnings
            @warn "It appears that the center body coordinates have been provided counter_clockwise rather than clockwise. Reversing direction for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing)."
        end
        reverse!(center_body_coordinates; dims=1)
    end

    # check z and r duct coordinates
    if i_know_what_im_doing &&
        abs(-(extrema(center_body_coordinates[:, 2])...)) >
       abs(-(extrema(center_body_coordinates[:, 1])...))
        if !silence_warnings
            @warn "It appears that the center body coordinates have been provided as (r,z) rather than (z,r). Reversing direction for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic reversing)."
        end
        reverse!(center_body_coordinates; dims=2)
    end

    ### --- CHECK FOR INPUT ERRORS --- ###

    # initialize error messages
    throw_error = false
    error_messages = ""
    error_count = 1

    # - Bodies - #

    #= TODO: things to check for duct_coordinates
             - no crossover? (how to check this?)
    =#

    # check for negative radial positions in duct_coordinates
    if any(duct_coordinates[:, 2] .< 0.0)
        if !silence_warnings
            @warn "Note that duct radial coordinates of duct must be positive. If you have not set the `autoshiftduct` option in the `Options` input to true, you should see an error shortly."
        end
    end

    # check we have the same lengths for duct coordinates
    if length(duct_coordinates[:, 1]) != length(duct_coordinates[:, 2])
        throw_error = true
        error_messages *= "\n\tError $(error_count): z and r coordinates of duct must have the same length"
        error_count += 1
    end

    # check for negative radial positions in center_body_coordinates
    if any(center_body_coordinates[:, 2] .< 0.0)
        throw_error = true
        error_messages *= "\n\tError $(error_count): Radial coordinates of center body must be positive."
        error_count += 1
    end

    # check we have the same lengths for center_body_coordinates
    if length(center_body_coordinates[:, 1]) != length(center_body_coordinates[:, 2])
        throw_error = true
        error_messages *= "\n\tError $(error_count): z and r coordinates of center body must have the same length"
        error_count += 1
    end

    # check dte_minus_cbte
    if duct_coordinates[1, 1] > center_body_coordinates[end, 1]
        if paneling_constants.dte_minus_cbte <= 0
            throw_error = true
            error_messages *= "\n\tError $(error_count): It appears that the dte_minus_cbte value is incorrect. If the duct trailing edge is behind the center body trailing edge, `dte_minus_cbte` should be positive."
            error_count += 1
        end
    elseif duct_coordinates[1, 1] < center_body_coordinates[end, 1]
        if paneling_constants.dte_minus_cbte >= 0
            throw_error = true
            error_messages *= "\n\tError $(error_count): It appears that the dte_minus_cbte value is incorrect. If the duct trailing edge is ahead of the center body trailing edge, `dte_minus_cbte` should be negative."
            error_count += 1
        end
    else
        if paneling_constants.dte_minus_cbte != 0
            throw_error = true
            error_messages *= "\n\tError $(error_count): It appears that the dte_minus_cbte value is incorrect. If the duct trailing edge is aligned with the center body trailing edge, `dte_minus_cbte` should be zero."
            error_count += 1
        end
    end

    # check number of panels relative to number of rotors and dte_minus_cbte
    if iszero(duct_coordinates[1, 1] == center_body_coordinates[end, 1])
        if length(paneling_constants.num_panels) != length(rotor.rotor_axial_position) + 1
            throw_error = true
            error_messages *= "\n\tError $(error_count): Length of vector `num_panels` should be one more than the length of vector `rotor_axial_position` when the duct and center_body trailing edges align."
            error_count += 1
        end
    else
        if length(paneling_constants.num_panels) != length(rotor.rotor_axial_position) + 2
            throw_error = true
            error_messages *= "\n\tError $(error_count): Length of vector `num_panels` should be two more than the length of vector `rotor_axial_position` when the duct and center_body trailing edges do not align."
            error_count += 1
        end
    end

    # other checks:
    # - rotor location is inside duct
    for rap in rotor.rotor_axial_position
        if rap < minimum(duct_coordinates[1, 1])
            throw_error = true
            error_messages *= "\n\tError $(error_count): Rotor is in front of duct leading edge. Rotor must be within duct."
            error_count += 1
        end
        if rap > duct_coordinates[end, 1]
            throw_error = true
            error_messages *= "\n\tError $(error_count): Rotor is behind duct trailing edge. Rotor must be within duct."
            error_count += 1
        end
        if rap < center_body_coordinates[1, 1]
            throw_error = true
            error_messages *= "\n\tError $(error_count): Rotor is in front of center_body leading edge. Rotor must be attached to center body."
            error_count += 1
        end
        if rap > center_body_coordinates[end, 1]
            throw_error = true
            error_messages *= "\n\tError $(error_count): Rotor is behind center_body trailing edge. Rotor must be attached to center body."
            error_count += 1
        end
    end

    if throw_error
        throw(error_messages)
        return nothing
    else
        return DuctedRotor(
            duct_coordinates, center_body_coordinates, rotor, paneling_constants
        )
    end
end
