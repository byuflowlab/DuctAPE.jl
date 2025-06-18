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

        return new(
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

        return new(
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
- `num_center_body_inlet_panels::Int` : The number of panels to use for the center_body inlet.
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
    wake_length::TF = 1.0

    function PanelingConstants(
        num_duct_inlet_panels::TI,
        num_center_body_inlet_panels::TI,
        num_panels::AbstractVector{TI},
        dte_minus_cbte::TFI,
        num_wake_sheets::TI,
        wake_length::TF=1.0,
    ) where {TI,TF,TFI}

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
            @error error_messages
            return nothing
        else
            return new(
                num_duct_inlet_panels,
                num_center_body_inlet_panels,
                num_panels,
                dte_minus_cbte,
                num_wake_sheets,
                wake_length,
            )
        end
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
- `Rhub::AbstractVector{Float}` : Dimensional hub radius of rotor. (may be changed if it does not match the radial position of the center_body geometry at the selected `rotor_axial_position`.
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

    function Rotor(
        B::Tb,
        rotor_axial_position::TRz,
        r::Tr,
        Rhub::TRh,
        Rtip::TRt,
        chords::Tc,
        twists::Tt,
        tip_gap::TTg,
        airfoils::Taf,
        is_stator::Tf,
        i_know_what_im_doing=false,
    ) where {
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

        # initialize error messages
        throw_error = false
        error_messages = ""
        error_count = 1

        if !i_know_what_im_doing && length(findall(t -> t > 1.75, twists)) > 2
            @warn "It looks like your input twist angles may be in degrees. Note that the required units for twist are radians. Converting to radians for you (set the `i_know_what_im_doing` keyword argument to true to disable automatic conversion)."
            twists .*= pi / 180.0
        end

        if length(unique(rotor_axial_position)) == length(rotor_axial_position)
            throw_error = true
            error_messages *= "\n\tError $(error_count): Cannot place rotors on top of eachother."
            error_count += 1
        end
        if Rhub > Rtip
            throw_error = true
            error_messages *= "\n\tError $(error_count): `Rtip` must be greater than `Rhub`."
            error_count += 1
        end
        if !all(r[2:end] .> r[1:(end - 1)])
            throw_error = true
            error_messages *= "\n\tError $(error_count): Radial positions, `r`, must be increasing across the blade"
            error_count += 1
        end
        if any(r .<= 0)
            throw_error = true
            error_messages *= "\n\tError $(error_count): Radial positions, `r`, must be postive, non-zero"
            error_count += 1
        end

        if throw_error
            @error error_messages
            return nothing
        else
            return new(
                isscalar(B) ? [B] : B,
                isscalar(rotor_axial_position) ? [rotor_axial_position] : rotor_axial_position,
                isscalar(r) ? [r] : r,
                isscalar(Rhub) ? [Rhub] : Rhub,
                isscalar(Rtip) ? [Rtip] : Rtip,
                isscalar(chords) ? [chords] : chords,
                isscalar(twists) ? [twists] : twists,
                isscalar(tip_gap) ? [tip_gap] : tip_gap,
                if typeof(airfoils) <:
                    Union{c4b.AFType,c4b.DTCascade,c4b.DFDCairfoil,c4b.ADM}
                    [airfoils]
                else
                    airfoils
                end,
                isscalar(is_stator) ? [is_stator] : is_stator,
            )
        end
    end
end

"""
    DuctedRotor(duct_coordinates, center_body_coordinates, rotor, paneling_constants)

# Arguments

- `duct_coordinates::AbstractMatrix` : The [z, r] coordinates of the duct geometry beginning at the inner (casing) side trailing edge and proceeding clockwise. Note that the duct geometry absolute radial position does not need to be included here if the `autoshiftduct` option is selected.
- `center_body_coordinates::AbstractMatrix` : The [z, r] coordinates of the center_body beginning at the leading edge and ending at the trailing edge. Note that the leading edge is assumed to be placed at a radial distance of 0.0 from the axis of rotation.
- `rotor::Rotor` : Rotor (and possibly stator) geometric paramters.
- `paneling_constants::PanelingConstants` : Constants used in re-paneling the geometry.
"""
struct DuctedRotor{Td<:AbstractMatrix,Tcb<:AbstractMatrix,Trp<:Rotor,Tpc<:PanelingConstants}
    duct_coordinates::Td
    center_body_coordinates::Tcb
    rotor::Trp
    paneling_constants::Tpc

    function DuctedRotor(
        duct_coordinates::Td,
        center_body_coordinates::Tcb,
        rotor::Trp,
        paneling_constants::Tpc,
    ) where {Td<:AbstractMatrix,Tcb<:AbstractMatrix,Trp<:Rotor,Tpc<:PanelingConstants}

        # initialize error messages
        throw_error = false
        error_messages = ""
        error_count = 1

        #= TODO: things to check for duct_coordinates
                 - correct direction
                 - all positive
                 - no crossover? (how to check this?)
                  @assert
        =#

        if any(duct_coordinates[:, 2] .< 0.0)
            throw_error = true
            error_messages *= "\n\tError $(error_count): Duct Coordinates must be positive."
            error_count += 1
        end
        if length(duct_coordinates[:, 1]) != length(duct_coordinates[:, 2])
            throw_error = true
            error_messages *= "\n\tError $(error_count): z and r coordinates of duct must have the same length"
            error_count += 1
        end

        if length(center_body_coordinates[:, 1]) != length(center_body_coordinates[:, 2])
            throw_error = true
            error_messages *= "\n\tError $(error_count): z and r coordinates of center body must have the same length"
            error_count += 1
        end

        #= TODO: things to check for center_body_coordinates
                 - correct direction
                 - all positive
                  @assert length(x) == length(y) "X and Y vectors must be the same length"
        =#

        #= TODO: things to check for paneling_constants
                 - number of rotors and body alignment vs npanel
                 - dte_minus_cbte vs coordinates
                 if iszero(dte_minus_cbte)
        zd = vcat(rotor_axial_position, cb_tez, cb_tez + wake_length)
        @assert length(num_panels) == length(rotor_axial_position) + 1 "Length of vector `num_panels` should be one more than the length of vector `rotor_axial_position` when the duct and center_body trailing edges align."
        elseif dte_minus_cbte < 0 #duct_tez < cb_tez
        zd = vcat(rotor_axial_position, duct_tez, cb_tez, cb_tez + wake_length)
        @assert length(num_panels) == length(rotor_axial_position) + 2 "Length of vector `num_panels` should be two more than the length of vector `rotor_axial_position` when the duct and center_body trailing edges do not align."
        else #dte_minus_cbte > 0 # duct_tez < cb_tez
        zd = vcat(rotor_axial_position, cb_tez, duct_tez, duct_tez + wake_length)
        @assert length(num_panels) == length(rotor_axial_position) + 2 "Length of vector `num_panels` should be two more than the length of vector `rotor_axial_position` when the duct and center_body trailing edges align."
        end
        =#

        # TODO: go find the various asserts and put them here as appropriate

        # other checks:
        # - rotor location is inside duct
        for rzl in rotor_axial_position
            @assert rzl > duct_lez "Rotor is in front of duct leading edge."
            @assert rzl < duct_tez "Rotor is behind duct trailing edge."
            @assert rzl > cb_lez "Rotor is in front of center_body leading edge."
            @assert rzl < cb_tez "Rotor is behind center_body trailing edge."
        end

        if throw_error
            @error error_messages
            return nothing
        else
            return new(Td, Tcb, Trp, Tpc)
        end
    end
end

