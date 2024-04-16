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

"""
    OperatingPoint(Vinf, rhoinf, muinf, asound, Omega)

Propulsor operating point information.

Note that the actual struct requires the inputs to be arrays, but there is a constructor function that will take in scalars and automatically build the array-based struct.

Also note that even though each field is required to be a vector, only `Omega` should have more than one entry, and only then if there are more than one rotor.  The purpose behind having vector rather than constant scalar inputs here is for ease of redefinition in an optimization setting when freestream design variables may be present.

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
- `fliplift::AbstractVector{Bool}` : flag to indicate if the airfoil lift values should be flipped or not.
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
- `rotorstator_parameters::OperatingPoint` : The operating point values.
- `operating_point::PanelingConstants` : Constants used in re-paneling the geometry.
- `paneling_constants::RotorStatorParameters` : Rotor (and possibly stator) geometric paramters.
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

function GaussLegendre(nsamples=8; silence_warnings=true)
    if silence_warnings && Bool((nsamples) % 2)
        @warn "Must have an even number of GaussLegendre sample points if using for panel self influence"
    end

    nodes, weights = FastGaussQuadrature.gausslegendre(nsamples)

    return GaussLegendre(linear_transform((-1, 1), (0, 1), nodes), weights ./ 2.0)
end

@kwdef struct IntegrationOptions{TN<:IntegrationMethod,TS<:IntegrationMethod}
    nominal::TN = GaussLegendre(8)
    singular::TS = GaussLegendre(8)
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
            additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoForwardDiff()),
        ),
    ]
    converged::AbstractVector{TB} = [false]
end

function ChainSolverOptions(multipoint)
    lm = length(multipoint)
    return ChainSolverOptions(;
        solvers=[
            NLsolveOptions(; algorithm=:anderson, atol=1e-12, converged=fill(false, lm)),
            MinpackOptions(; atol=1e-12, converged=fill(false, lm)),
            NonlinearSolveOptions(;
                algorithm=SimpleNonlinearSolve.SimpleNewtonRaphson,
                atol=1e-12,
                additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoForwardDiff()),
                converged=fill(false, lm),
            ),
        ],
        converged=fill(false, lm),
    )
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
    TB,
    TBwo,
    TF,
    TI,
    TSf,
    TSt,
    Tin,
    TIo<:IntegrationOptions,
    TSo<:SolverOptionsType,
    WS<:GridSolverOptionsType,
}
    # - General Options - #
    verbose::TB = false
    silence_warnings::TB = true
    multipoint_index::TI = [1]
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
    write_outputs::TBwo = [false]
    outfile::TSf = ["outputs.jl"]
    checkoutfileexists::TB = false
    output_tuple_name::TSt = ["outs"]
    # - Solving Options - #
    grid_solver_options::WS = GridSolverOptions()
    solver_options::TSo = ChainSolverOptions()
end

"""
    set_options(; kwargs...)
    set_options(multipoint; kwargs...)

Set the options for DuctAPE to use.

Note that the vast majority of the available options are defined through keyword arguments.  See the documentation for the various option types for more information.

# Arguments

- `multipoint::AbstractVector{OperatingPoint}` : a vector of operating points to use if running a multi-point analysis.
"""
function set_options(; kwargs...)
    return Options(; kwargs...)
end

function set_options(
    multipoint::AbstractVector{TM};
    write_outputs=nothing,
    outfile=nothing,
    output_tuple_name=nothing,
    kwargs...,
) where {TM<:OperatingPoint}
    lm = length(multipoint)

    if isnothing(outfile)
        outfile = ["outputs" * lpad(i, 3, "0") * ".jl" for i in 1:lm]
    end

    if isnothing(output_tuple_name)
        output_tuple_name = ["outs" for i in 1:lm]
    end

    if isnothing(write_outputs)
        write_outputs = [false for i in 1:lm]
    end

    return set_options(;
        solver_options=ChainSolverOptions(multipoint),
        write_outputs=write_outputs,
        outfile=outfile,
        output_tuple_name=output_tuple_name,
        kwargs...,
    )
end

function DFDC_options(;
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

######################################################################
#                                                                    #
#                        FUNCTION SET HEADER                         #
#                                                                    #
######################################################################

function promote_propulosor_type(p)
    return promote_type(
        eltype(p.duct_coordinates),
        eltype(p.centerbody_coordinates),
        eltype(p.operating_point.Vinf),
        eltype(p.operating_point.rhoinf),
        eltype(p.operating_point.muinf),
        eltype(p.operating_point.asound),
        eltype(p.operating_point.Omega),
        eltype(p.rotorstator_parameters.B),
        eltype(p.rotorstator_parameters.rotorzloc),
        eltype(p.rotorstator_parameters.r),
        eltype(p.rotorstator_parameters.Rhub),
        eltype(p.rotorstator_parameters.Rtip),
        eltype(p.rotorstator_parameters.chords),
        eltype(p.rotorstator_parameters.twists),
    )
end
