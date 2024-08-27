"""
    OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)

Propulsor operating point information.

# Arguments

- `Vinf::AbstractVector{Float}` : Freestream velocity magnitude (which is only in the axial direction).
- `rhoinf::AbstractVector{Float}` : Freestream density
- `muinf::AbstractVector{Float}` : Freestream viscosity
- `asound::AbstractVector{Float}` : Freestream speed of sound
- `Omega::AbstractVector{Float}` : Rotor rototation rate(s)
"""
struct OperatingPoint{
    Tv<:AbstractVector,
    Tr<:AbstractVector,
    Tm<:AbstractVector,
    Ta<:AbstractVector,
    To<:AbstractVector,
}
    Vinf::Tv
    rhoinf::Tr
    muinf::Tm
    asound::Ta
    Omega::To
end

function OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)
    return OperatingPoint(
        isscalar(Vinf) ? [Vinf] : Vinf,
        isscalar(rhoinf) ? [rhoinf] : rhoinf,
        isscalar(muinf) ? [muinf] : muinf,
        isscalar(asound) ? [asound] : asound,
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
        nduct_inlet,
        ncenterbody_inlet,
        npanels,
        dte_minus_cbte,
        nwake_sheets,
        wake_length=1.0,
    )

Constants used in re-paneling geometry.

Note that unlike other input structures, this one, in general, does not define fields as vectors.  This is because these values should not change throughout an optimization, even if the geometry may change.  Otherwise, discontinuities could be experienced.

# Arguments

- `nduct_inlet::Int` : The number of panels to use for the duct inlet (this number is used for both the casing and nacelle re-paneling)
- `ncenterbody_inlet::Int` : The number of panels to use for the centerbody inlet.
- `npanels::AbstractVector{Int}` : A vector containing the number of panels between discrete locations inside the wake. Specifically, the number of panels between the rotors, between the last rotor and the first body trailing edge, between the body trailing edges (if different), and between the last body trailing edge and the end of the wake.  The length of this vector should be N+1 (where N is the number of rotors) if the duct and centerbody trailing edges are aligned, and N+2 if not.
- `dte_minus_cbte::Float` : An indicator concerning the hub and duct trailing edge relative locations. Should be set to -1 if the duct trailing edge axial position minus the centerbody trailing edge axial position is negative, +1 if positive (though any positive or negative number will suffice), and zero if the trailing edges are aligned.
- `nwake_sheets::Int` : The number of wake sheets to use. Note this will also be setting the number of blade elements to use.
- `wake_length::Float=1.0` : Non-dimensional (based on the length from the foremost body leading edge and the aftmost body trailing edge) length of the wake extending behind the aftmost body trailing edge.
"""
@kwdef struct PanelingConstants{TI,TF,TFI}
    nduct_inlet::TI
    ncenterbody_inlet::TI
    npanels::AbstractVector{TI}
    dte_minus_cbte::TFI
    nwake_sheets::TI
    wake_length::TF = 1.0
end

"""
    RotorStatorParameters(
        B, rotorzloc, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, fliplift
    )

Composite type containing the rotor(s) geometric properties.

Note that the actual struct requires the inputs to be arrays, but there is a constructor function that will take in scalars and automatically build the array-based struct.

# Arguments
- `B::AbstractVector{Float}` : The number of blades for each rotor. May not be an integer, but usually is.
- `rotorzloc::AbstractVector{Float}` : Dimensional, axial position of each rotor.
- `r::AbstractArray{Float}` : Non-dimensional radial locations of each blade element.
- `Rhub::AbstractVector{Float}` : Dimensional hub radius of rotor. (may be changed if it does not match the radial position of the centerbody geometry at the selected `rotorzloc`.
- `Rtip::AbstractVector{Float}` : Dimensional tip radius of rotor. Is used to determine the radial position of the duct if the `autoshiftduct` option is selected.
- `chords::AbstractArray{Float}` : Dimensional chord lengths of the blade elements.
- `twists::AbstractArray{Float}` : Blade element angles, in radians.
- `tip_gap::AbstractVector{Float}` : Currently unused, do not set to anything other than zeros.
- `airfoils::AbstractArray{AFType}` : Airfoil types describing the airfoil polars for each blade element. Currently only fully tested with `C4Blade.DFDCairfoil` types.
- `fliplift::AbstractVector{Bool}` : Flag to indicate if the airfoil lift values should be flipped or not.
"""
struct RotorStatorParameters{
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
    rotorzloc::TRz
    r::Tr
    Rhub::TRh
    Rtip::TRt
    chords::Tc
    twists::Tt
    tip_gap::TTg
    airfoils::Taf
    fliplift::Tf
end

function RotorStatorParameters(
    B, rotorzloc, r, Rhub, Rtip, chords, twists, tip_gap, airfoils, fliplift
)
    return RotorStatorParameters(
        isscalar(B) ? [B] : B,
        isscalar(rotorzloc) ? [rotorzloc] : rotorzloc,
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
        isscalar(fliplift) ? [fliplift] : fliplift,
    )
end

"""
    Propulsor(duct_coordinates, centerbody_coordinates, rotorstator_parameters, operating_point, paneling_constants, reference_parameters)

# Arguments

- `duct_coordinates::AbstractMatrix` : The [z, r] coordinates of the duct geometry beginning at the inner (casing) side trailing edge and proceeding clockwise. Note that the duct geometry absolute radial position does not need to be included here if the `autoshiftduct` option is selected.
- `centerbody_coordinates::AbstractMatrix` : The [z, r] coordinates of the centerbody beginning at the leading edge and ending at the trailing edge. Note that the leading edge is assumed to be placed at a radial distance of 0.0 from the axis of rotation.
- `operating_point::OperatingPoint` : The operating point values.
- `paneling_constants::PanelingConstants` : Constants used in re-paneling the geometry.
- `rotorstator_parameters::RotorStatorParameters` : Rotor (and possibly stator) geometric paramters.
- `reference_parameters::ReferenceParameters` : Reference Parameters.
"""
struct Propulsor{
    Td<:AbstractMatrix,
    Tcb<:AbstractMatrix,
    Top<:OperatingPoint,
    Tpc<:PanelingConstants,
    Trp<:RotorStatorParameters,
    Tref<:ReferenceParameters,
}
    duct_coordinates::Td
    centerbody_coordinates::Tcb
    rotorstator_parameters::Trp
    operating_point::Top
    paneling_constants::Tpc
    reference_parameters::Tref
end

"""
TODO: write this function and have it do all the checks to make sure that the user inputs are going to work.
"""
function verify_input(propulsor)
    # TODO: check number of rotors vs npanel
    # TODO: check rotorzloc is sorted
    # TODO: check dte_minus_cbte vs coordinates
    # TODO: go find all the various asserts and put them here
end
