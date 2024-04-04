#=
This file contains types and functions used for running DuctAPE in a gradient-based optimization setting
=#

#---------------------------------#
#          ABSTRACT TYPES         #
#---------------------------------#

# - Types for dispatching CSOR residual and convergence - #
abstract type ConvergenceType end

# - Solver Options - #
abstract type SolverOptionsType end
abstract type ExternalSolverOptions <: SolverOptionsType end

# - Wake Solver Options - #
abstract type GridSolverOptionsType end

# - Quadrature Types - #
abstract type IntegrationMethod end

#---------------------------------#
#           INPUT TYPES           #
#---------------------------------#
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
        isscalar(Vref) ? [Vref] : Vref, isscalar(Rref) ? [Rref] : Rref
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

#---------------------------------#
#         QUADRATURE TYPES        #
#---------------------------------#

@kwdef struct Romberg{TI,TF,TP,TV} <: IntegrationMethod
    max_subdivisions::TI = 10
    atol::TF = 1e-6
end

function set_romberg_options(; max_subdivisions=10, atol=1e-6)
    return Romberg(; max_subdivisions, atol)
end

@kwdef struct GaussKronrod{TI,TF} <: IntegrationMethod
    order::TI = 5
    maxevals::TI = 1000
    atol::TF = 1e-12
end

struct GaussLegendre{TN,TW} <: IntegrationMethod
    sample_points::TN
    weights::TW
end

function GaussLegendre(nsamples=20; silence_warnings=true)
    if silence_warnings && Bool((nsamples) % 2)
        @warn "Must have an even number of GaussLegendre sample points if using for panel self influence"
    end

    nodes, weights = FastGaussQuadrature.gausslegendre(nsamples)

    return GaussLegendre(linear_transform((-1, 1), (0, 1), nodes), weights ./ 2.0)
end

@kwdef struct IntegrationOptions{TN,TS}
    nominal::TN = GaussLegendre(20)
    singular::TS = GaussLegendre(20)
end

#---------------------------------#
#          OPTION TYPES           #
#---------------------------------#

# - Types for CSOR Convergence dispatch - #
struct Relative <: ConvergenceType end
struct Absolute <: ConvergenceType end

# - Solver Options - #
@kwdef struct CSORSolverOptions{TF,TS,TB,TC<:ConvergenceType} <: SolverOptionsType
    # Defaults are DFDC hard-coded values
    verbose::TB = false
    maxiter::TF = 1e2
    nrf::TF = 0.4
    bt1::TF = 0.2
    bt2::TF = 0.6
    pf1::TF = 0.4
    pf2::TF = 0.5
    btw::TF = 0.6
    pfw::TF = 1.2
    relaxation_schedule::TS = [
        reverse!([1e10; 1e-13; 1e-13; 0.0]), reverse!([0.0; 0.0; 1.0; 1.0])
    ]
    f_circ::TF = 1e-3
    f_dgamw::TF = 2e-4
    convergence_type::TC = Relative()
    Vconv::TF = 1.0
    converged::AbstractVector{TB} = [false]
end

# @kwdef struct SolverOptions{TA,TB,TF,TI,TN,TTm,TTr} <: SolverOptionsType
@kwdef struct NonlinearSolveOptions{TA,TB,TF,TI} <: ExternalSolverOptions
    # Algorithm Options
    nlsolve_algorithm::TA = SimpleNonlinearSolve.SimpleDFSane
    # Iteration Controls
    nlsolve_abstol::TF = 1e-10
    nlsolve_maxiters::TI = 100
    converged::AbstractVector{TB} = [false]
end

@kwdef struct NLsolveOptions{TSym,TF,TI,TB,Tls,Tlsk} <: ExternalSolverOptions
    # - Options for overall solve - #
    # TODO: generalize the newton part of this to use NonlinearSolve.jl framework.
    # TODO: consider a tighter default convergence tolerance.
    # nlsolve parameters
    nlsolve_method::TSym = :newton
    nlsolve_ftol::TF = 1e-8 #1e-8 is nlsolve default
    nlsolve_iteration_limit::TI = 20 #1000 is nlsolve default
    nlsolve_show_trace::TB = false
    # line search parameters
    nlsolve_linesearch_method::Tls = LineSearches.MoreThuente
    nlsolve_linesearch_kwargs::Tlsk = (;)
    converged::AbstractVector{TB} = [false]
end

# - Elliptic Grid Solver Options - #
@kwdef struct SLORGridSolverOptions{TF,TI,TB} <: GridSolverOptionsType
    max_wake_relax_iter::TI = 100
    wake_relax_tol::TF = 1e-9
    converged::AbstractVector{TB} = [false]
end

@kwdef struct GridSolverOptions{TSym,TF,TI,TB} <: GridSolverOptionsType
    # elliptic grid solve options
    # TODO: generalize the newton part of this to use NonlinearSolve.jl framework.
    wake_nlsolve_method::TSym = :newton
    wake_nlsolve_autodiff::TSym = :forward
    wake_nlsolve_ftol::TF = 1e-14
    wake_max_iter::TI = 100
    converged::AbstractVector{TB} = [false]
    max_wake_relax_iter::TI = 20
    wake_relax_tol::TF = 1e-9
end

# - Full Option Set - #
"""
"""
@kwdef struct Options{
    TB,TF,TS,Tin,TIo<:IntegrationOptions,TSo<:SolverOptionsType,WS<:GridSolverOptionsType
}
    # - General Options - #
    verbose::TB = false
    silence_warnings::TB = true
    # - Geometry Re-interpolation and generation options - #
    finterp::Tin = FLOWMath.akima
    autoshiftduct::TB = true
    lu_decomp_flag::TB = false
    # paneling options
    itcpshift::TF = 0.05
    axistol::TF = 1e-15
    tegaptol::TF = 1e1 * eps()
    # - Integration Options - #
    integration_options::TIo = IntegrationOptions()
    # - Post-processing Options - #
    write_outputs::TB = false
    outfile::TS = "outputs.jl"
    checkoutfileexists::TB = false
    output_tuple_name::TS = "outs"
    # - Solving Options - #
    grid_solver_options::WS = GridSolverOptions()
    solver_options::TSo = NonlinearSolveOptions()
end

"""
"""
function set_options(; kwargs...)
    return Options(; kwargs...)
end

"""
"""
function quicksolve_options(;
    grid_solver_options=SLORGridSolverOptions(),
    solver_options=CSORSolverOptions(),
    kwargs...,
)
    return Options(;
        grid_solver_options=SLORGridSolverOptions(),
        solver_options=CSORSolverOptions(),
        kwargs...,
    )
end

