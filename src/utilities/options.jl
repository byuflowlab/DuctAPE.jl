#---------------------------------#
#          ABSTRACT TYPES         #
#---------------------------------#

# - Types for dispatching CSOR residual and convergence - #
abstract type ConvergenceType end

# - Solver Options - #
abstract type SolverOptionsType end
abstract type ExternalSolverOptions <: SolverOptionsType end
abstract type PolyAlgorithmOptions <: SolverOptionsType end

# - Wake Solver Options - #
abstract type GridSolverOptionsType end

# - Quadrature Types - #
abstract type IntegrationMethod end


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
    TB,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}
} <: PolyAlgorithmOptions
    solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
    ]
    converged::AbstractVector{TB} = [false]
end

@kwdef struct ChainSolverOptions{TB,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}} <:
              PolyAlgorithmOptions
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
