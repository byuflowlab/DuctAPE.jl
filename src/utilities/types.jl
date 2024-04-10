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
abstract type MultiSolverOptions <: SolverOptionsType end

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

struct ReferenceParameters{Tv<:AbstractVector,Tr<:AbstractVector}
    Vref::Tv
    Rref::Tr
end

function ReferenceParameters(Vref, Rref)
    return ReferenceParameters(
        isscalar(Vref) ? [Vref] : Vref, isscalar(Rref) ? [Rref] : Rref
    )
end

struct PanelingConstants{TI,TF}
    nduct_inlet::TI
    ncenterbody_inlet::TI
    npanels::AbstractVector{TI}
    dte_minus_cbte::TI
    nwake_sheets::TI
    wake_length::TF
end

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

@kwdef struct Romberg{TI,TF} <: IntegrationMethod
    max_subdivisions::TI = 10
    atol::TF = 1e-6
end

@kwdef struct GaussKronrod{TI,TF} <: IntegrationMethod
    order::TI = 7
    maxevals::TI = 10^7
    atol::TF = 0.0
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

@kwdef struct IntegrationOptions{TN<:IntegrationMethod,TS<:IntegrationMethod}
    nominal::TN = GaussLegendre(20)
    singular::TS = GaussLegendre(20)
end

#---------------------------------#
#          SOLVER TYPES           #
#---------------------------------#

# - Types for CSOR Convergence dispatch - #
struct Relative <: ConvergenceType end
struct Absolute <: ConvergenceType end

##### ----- Fixed Point Solvers ----- #####

# - CSOR Options - #
@kwdef struct CSORSolverOptions{TF,TS,TB,TC<:ConvergenceType} <: SolverOptionsType
    # Defaults are DFDC hard-coded values
    verbose::TB = false
    iteration_limit::TF = 1e2
    nrf::TF = 0.4
    bt1::TF = 0.2
    bt2::TF = 0.6
    pf1::TF = 0.4
    pf2::TF = 0.5
    btw::TF = 0.6
    pfw::TF = 1.2
    relaxation_schedule::TS = [
        reverse!([1e10; 1e-13; 1e-14; 0.0]), reverse!([0.0; 0.0; 1.0; 1.0])
    ]
    f_circ::TF = 1e-3
    f_dgamw::TF = 2e-4
    convergence_type::TC = Relative()
    Vconv::TF = 1.0
    converged::AbstractVector{TB} = [false]
end

@kwdef struct FixedPointOptions{TB,TF,TI} <: ExternalSolverOptions
    iteration_limit::TI = 1000
    vel::TF = 0.9
    ep::TF = 0.01
    atol::TF = 1e-12
    converged::AbstractVector{TB} = [false]
end

@kwdef struct SpeedMappingOptions{TB,TI,TF,TL,TSm,TU} <: ExternalSolverOptions
    orders::AbstractVector{TI} = [3, 2]
    sig_min::TSm = 0 # maybe set to 1?
    stabilize::TB = false # stabilizes before extrapolation
    check_obj::TB = false # checks for inf's and nan's and starts from previous finite point
    atol::TF = 1e-10 # convergence tolerance
    iteration_limit::TI = 1000 # number of "iterations"
    time_limit::TF = Inf
    lower::TL = nothing # box lower bounds
    upper::TU = nothing # box upper bounds
    buffer::TF = 0.01 # if using bounds, buffer brings x inside bounds by buffer amountd
    Lp::TF = Inf # p value for p-norm for convergence criteria
    converged::AbstractVector{TB} = [false]
end

##### ----- Quasi-Newton Solvers ----- #####

@kwdef struct MinpackOptions{TB,TF,TI,TS} <: ExternalSolverOptions
    algorithm::TS = :hybr
    atol::TF = 1e-12
    iteration_limit::TI = 100
    converged::AbstractVector{TB} = [false]
end

@kwdef struct SIAMFANLEOptions{TB,TF,TH,TI,TK} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TH = SIAMFANLEquations.nsoli
    atol::TF = 1e-10
    rtol::TF = 0.0
    iteration_limit::TI = 1000
    linear_iteration_limit::TI = 2
    additional_kwargs::TK = (;)
    # additional_kwargs::TK = (; delta0=1e-3)
    converged::AbstractVector{TB} = [false]
end

##### ----- Newton+ Solvers ----- #####
# NOTE: these also have fixed-point options

@kwdef struct NonlinearSolveOptions{TA,TB,TF,TI,TT} <: ExternalSolverOptions
    # Algorithm Options
    algorithm::TA = SimpleNonlinearSolve.SimpleNewtonRaphson
    additional_kwargs::TT = (;)
    # Iteration Controls
    atol::TF = 1e-10
    iteration_limit::TI = 100
    converged::AbstractVector{TB} = [false]
end

@kwdef struct NLsolveOptions{TSym,TF,TI,TB,Tls,Tlsk} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TSym = :anderson
    atol::TF = 1e-12
    iteration_limit::TI = 100
    # line search parameters
    linesearch_method::Tls = LineSearches.MoreThuente
    linesearch_kwargs::Tlsk = (;)
    converged::AbstractVector{TB} = [false]
end

##### ----- Poly-Algorithm Solvers ----- #####

@kwdef struct CompositeSolverOptions{
    TB,TS<:Union{ExternalSolverOptions,MultiSolverOptions}
} <: MultiSolverOptions
    solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
    ]
    converged::AbstractVector{TB} = [false]
end

@kwdef struct ChainSolverOptions{TB,TS<:Union{ExternalSolverOptions,MultiSolverOptions}} <:
              MultiSolverOptions
    solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
        MinpackOptions(; atol=1e-12),
        NonlinearSolveOptions(;
            algorithm=SimpleNonlinearSolve.SimpleNewtonRaphson,
            atol=1e-12,
            additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoPolyesterForwardDiff()),
        ),
    ]
    converged::AbstractVector{TB} = [false]
end

#---------------------------------#
#   ELLIPTIC GRID SOLVER TYPES    #
#---------------------------------#
@kwdef struct SLORGridSolverOptions{TF,TI,TB} <: GridSolverOptionsType
    relaxation_iteration_limit::TI = 100
    relaxation_atol::TF = 1e-9
    converged::AbstractVector{TB} = [false]
end

@kwdef struct GridSolverOptions{TSym,TF,TI,TB} <: GridSolverOptionsType
    # elliptic grid solve options
    relaxation_iteration_limit::TI = 20
    relaxation_atol::TF = 1e-9
    algorithm::TSym = :newton
    autodiff::TSym = :forward
    atol::TF = 1e-14
    iteration_limit::TI = 10
    converged::AbstractVector{TB} = [false]
end

#---------------------------------#
#         OPTION SET TYPES        #
#---------------------------------#
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
    solver_options::TSo = CompositeSolverOptions()
end

function set_options(; kwargs...)
    return Options(; kwargs...)
end

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
