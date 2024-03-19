#=
This file contains types and functions used for running DuctAPE in a gradient-based optimization setting
=#

# - Functions to check if the input is a scalar
import Base.BroadcastStyle
isscalar(x::T) where {T} = isscalar(T)
isscalar(::Type{T}) where {T} = BroadcastStyle(T) isa Broadcast.DefaultArrayStyle{0}

"""
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
"""
struct ReferenceParameters{Tv<:AbstractVector,Tr<:AbstractVector}
    Vref::Tv
    Rref::Tr
end

function ReferenceParameters(Vref, Rref)
    return ReferenceParameters(
        isscalar(Vref) ? [Vref] : Vref, isscalar(rhoinf) ? [Rref] : Rref
    )
end

"""
"""
struct PanelingConstants{TI,TF}
    nduct_inlet::TI
    ncenterbody_inlet::TI
    npanels::AbstractVector{TI}
    dte_minus_cbte::TI
    nwake_sheets::TI
    wake_length::TF
end

"""
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

function RotorStatorParameters(Vref, Rref)
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
        isscalar(fliplif) ? [fliplift] : fliplift,
    )
end

"""
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

@kwdef struct QuickSolve{TF,TB}
    # Defaults are DFDC hard-coded values
    nosource::TB = true #TODO: remove this and associated functionality
    verbose::TB = false
    maxiter::TF = 1e2
    nrf::TF = 0.4
    bt1::TF = 0.2
    bt2::TF = 0.6
    pf1::TF = 0.4
    pf2::TF = 0.5
    btw::TF = 0.6
    pfw::TF = 1.2
    f_circ::TF = 1e-3
    f_dgamw::TF = 2e-4
end

@kwdef struct QuickWake{TF,TI,TB}
    wake_nlsolve_ftol::TF = 1e-9
    wake_max_iter::TI = 100
    wake_converged::AbstractVector{TB} = [false]
end

@kwdef struct NewtonWake{TF,TI,TB}
    # elliptic grid solve options
    wake_nlsolve_ftol::TF = 1e-14
    wake_max_iter::TI = 100
    wake_converged::AbstractVector{TB} = [false]
    max_wake_relax_iter::TI = 3
    wake_relax_tol::TF = 1e-14
end

@kwdef struct NewtonSolve{TY,TF,TI,TB,Tls,Tlsk}
    # - Options for overall solve - #
    # nlsolve parameters
    nlsolve_method::TY = :newton
    nlsolve_autodiff::TY = :forward
    nlsolve_ftol::TF = 1e-8 #1e-8 is nlsolve default
    nlsolve_iteration_limit::TI = 50 #1000 is nlsolve default
    nlsolve_show_trace::TB = false
    # line search parameters
    nlsolve_linesearch_method::Tls = LineSearches.MoreThuente
    nlsolve_linesearch_kwargs::Tlsk = (;)
    nlsolve_converged::AbstractVector{TB} = [false]
end

# struct QNewtonSolve{}
# end

"""
"""
@kwdef struct Options{TB,TF,TS,Tin,TSo,WS}
    # - General Options - #
    verbose::TB = false
    silence_warnings::TB = true
    # - Geometry Re-interpolation and generation options - #
    finterp::Tin = FLOWMath.akima
    autoshiftduct::TB = true
    # paneling options
    itcpshift::TF = 0.05
    axistol::TF = 1e-15
    tegaptol::TF = 1e1 * eps()
    # - Post-processing Options - #
    write_outputs::TB = false
    outfile::TS = "outputs.jl"
    checkoutfileexists::TB = false
    output_tuple_name::TS = "outs"
    # - Solving Options - #
    wake_options::WS = NewtonWake()
    solve_options::TSo = NewtonSolve()
end

"""
"""
function set_options(; kwargs...)
    return Options(; kwargs...)
end

"""
"""
function quicksolve_options(;
    wake_options=QuickWake(), solve_options=QuickSolve(), kwargs...
)
    return Options(; wake_options=QuickWake(), solve_options=QuickSolve(), kwargs...)
end
