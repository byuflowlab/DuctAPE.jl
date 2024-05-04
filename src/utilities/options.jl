#---------------------------------#
#          ABSTRACT TYPES         #
#---------------------------------#

# - Types for dispatching CSOR residual and convergence - #
"""
    abstract type ConvergenceType

Used in dispatching the CSOR (controlled successive over relaxation) residual as relative or absolute.
"""
abstract type ConvergenceType end

# - Solver Options - #
"""
    abstract type SolverOptionsType

Used for solver dispatch.
"""
abstract type SolverOptionsType end

"""
    abstract type ExternalSolverOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type ExternalSolverOptions <: SolverOptionsType end

"""
    abstract type PolyAlgorithmOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type PolyAlgorithmOptions <: SolverOptionsType end

# - Wake Solver Options - #
"""
    abstract type GridSolverOptionsType

Used for elliptic grid solver dispatch
"""
abstract type GridSolverOptionsType end

# - Quadrature Types - #
"""
    abstract type IntegrationMethod

Used in integration method dispatch
"""
abstract type IntegrationMethod end

#---------------------------------#
#         QUADRATURE TYPES        #
#---------------------------------#

"""
    struct Romberg{TF,TI} <: IntegrationMethod

Options for Romberg integration method

# Fields
- `max_subdivisions::TI = 10` : maximum number of subdivisions. Note, total number of internvals is 2^N, where N is number of subdivisions.
- `atol::TF = 1e-6` : absolute error tolerance.
"""
@kwdef struct Romberg{TF,TI} <: IntegrationMethod
    max_subdivisions::TI = 10
    atol::TF = 1e-6
end

"""
    struct GaussKronrod{TF,TI} <: IntegrationMethod

Options for Gauss-Kronrod integration method

# Fields
- `order::TI = 7` : order of Legendre polynomial to use on each interval
- `maxevales::TI = 10^7` : maximum number of evaluations in the adaptive method
- `atol::TF = 0.0` : absolute error tolerance. (note, if zero, QuadGK uses sqrt(eps()) relative tolerance).
"""
@kwdef struct GaussKronrod{TF,TI} <: IntegrationMethod
    order::TI = 7
    maxevals::TI = 10^7
    atol::TF = 0.0
end

"""
    struct GaussLegendre{TN,TW} <: IntegrationMethod

Options for Gauss-Legendre integration method

# Fields
- `sample_points::TN` : Sample Points
- `weights::TW` : Gauss weights
"""
struct GaussLegendre{TN,TW} <: IntegrationMethod
    sample_points::TN
    weights::TW
end

"""
    GaussLegendre(nsamples=8; silence_warnings=true)

Helper function that calls `FastGaussQuadrature.jl` to calculate the sample point position and wieght values for Gauss-Legendre quadrature.

# Arguments
- `nsamples::Int=8` : number of points to use on each panel. Note that if the number is not even, then you will have problems with self-induced panels.

# Keyword Arguments
- `silence_warnings::Bool=true` : if false, and you input an odd number of points, a warning will print telling you that's a problem.

# Returns
- `gl_options::GaussLegendre` : A GaussLegendre object defined based on the selected number of samples to use for integration.
"""
function GaussLegendre(nsamples=8; silence_warnings=true)
    if silence_warnings && Bool((nsamples) % 2)
        @warn "Must have an even number of GaussLegendre sample points if using for panel self influence"
    end

    nodes, weights = FastGaussQuadrature.gausslegendre(nsamples)

    return GaussLegendre(linear_transform((-1, 1), (0, 1), nodes), weights ./ 2.0)
end

"""
    struct IntegrationOptions{TN<:IntegrationMethod,TS<:IntegrationMethod}

A struct used to hold the integration options for both the nominal and singular cases.

# Fields
- `nominal::IntegrationMethod=GaussLegendre(8)` : the integration options to use for the nominal case.
- `singular::IntegrationMethod=GaussLegendre(8)` : the integration options to use for the self-induced case.
"""
@kwdef struct IntegrationOptions{TN<:IntegrationMethod,TS<:IntegrationMethod}
    nominal::TN = GaussLegendre(8)
    singular::TS = GaussLegendre(8)
end

#---------------------------------#
#          SOLVER TYPES           #
#---------------------------------#

# - Types for CSOR Convergence dispatch - #
"""
    struct Relative <: ConvergenceType

Used to dispatch the relative residual for CSOR (controlled successive over relaxation) method
"""
struct Relative <: ConvergenceType end

"""
    struct Absolute <: ConvergenceType

Used to dispatch the absolute residual for CSOR (controlled successive over relaxation)  method
"""
struct Absolute <: ConvergenceType end

##### ----- Fixed Point Solvers ----- #####

# - CSOR Options - #
"""
    struct CSORSolverOptions{TB,TC<:ConvergenceType,TF,TS} <: SolverOptionsType

Type containing all the options for the CSOR (controlled successive over relaxation) solver.

Note that the defaults match DFDC with the exception of the relaxation schedule, which is an experimental feature.

# Fields
- `var::type` :
- `verbose::TB = false` : flag to print verbose statements
- `iteration_limit::TF = 1e2` : maximum number of iterations
- `nrf::TF = 0.4` : nominal relaxation factor
- `bt1::TF = 0.2` : backtracking factor 1
- `bt2::TF = 0.6` : backtracking factor 2
- `pf1::TF = 0.4` : press forward factor 1
- `pf2::TF = 0.5` : press forward factor 2
- `btw::TF = 0.6` : backtracking factor for wake
- `pfw::TF = 1.2` : press forward factor for wake
- `relaxation_schedule::TS = [[0.0;1e-14;1e-13;1e10]), [1.0;1.0;0.0;0.0])]` : values used in spline definition for scaling the relaxation factors (second vector) after various convergence values (first vector).
- `f_circ::TF = 1e-3` : convergence tolerance for rotor circulation
- `f_dgamw::TF = 2e-4` : convergence tolerance for wake vortex strength
- `convergence_type::TC = Relative()` : dispatch for relative or absolute convergence criteria.
- `Vconv::AbstractVector{TF} = [1.0]` : velocity used in relative convergence criteria (should be set to Vref).
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct CSORSolverOptions{TB,TC<:ConvergenceType,TF,TI,TS} <: SolverOptionsType
    # Defaults are DFDC hard-coded values
    verbose::TB = false
    iteration_limit::TF = 1000
    nrf::TF = 0.4
    bt1::TF = 0.2
    bt2::TF = 0.6
    pf1::TF = 0.4
    pf2::TF = 0.5
    btw::TF = 0.6
    pfw::TF = 1.2
    relaxation_schedule::TS = [
        reverse!([1e10; 1e-16; 1e-17; 0.0]), reverse!([0.0; 0.0; 1.0; 1.0])
    ]
    f_circ::TF = 1e-3
    f_dgamw::TF = 2e-4
    convergence_type::TC = Relative()
    Vconv::AbstractVector{TF} = [1.0]
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct FixedPointOptions{TB,TF,TI} <: ExternalSolverOptions

Options for the FixedPoint.jl package solver

# Fields
- `iteration_limit::TF = 1000` : maximum number of iterations
- `vel::TF = 0.9` : vel keyword argument, default is package default
- `ep::TF = 0.01` : ep keyword argument, default is package default
- `atol::TF = 1e-12` : absolute convergence tolerance
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct FixedPointOptions{TB,TF,TI} <: ExternalSolverOptions
    iteration_limit::TI = 1000
    vel::TF = 0.9
    ep::TF = 0.01
    atol::TF = 1e-12
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct SpeedMappingOptions{TB,TF,TI,TL,TSm,TU} <: ExternalSolverOptions

Options for the SpeedMapping.jl package solver

# Fields
- `orders::AbstractVector{TI} = [3, 2]
- `sig_min::TSm = 0` : maybe set to 1?
- `stabilize::TB = false` : stabilizes before extrapolation
- `check_obj::TB = false` : checks for inf's and nan's and starts from previous finite point
- `atol::TF = 1e-10` : absolute convergence tolerance
- `iteration_limit::TF = 1000` : maximum number of iterations
- `time_limit::TF = Inf` : time limit in seconds
- `lower::TL = nothing` : box lower bounds
- `upper::TU = nothing` : box upper bounds
- `buffer::TF = 0.01` : if using bounds, buffer brings x inside bounds by buffer amountd
- `Lp::TF = Inf` : p value for p-norm for convergence criteria
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct SpeedMappingOptions{TB,TF,TI,TL,TSm,TU} <: ExternalSolverOptions
    orders::AbstractVector{TI} = [3, 2]
    sig_min::TSm = 0 # maybe set to 1?
    stabilize::TB = false # stabilizes before extrapolation
    check_obj::TB = false # checks for inf's and nan's and starts from previous finite point
    atol::TF = 1e-12 # convergence tolerance
    iteration_limit::TI = 1000 # number of "iterations"
    time_limit::TF = Inf
    lower::TL = nothing # box lower bounds
    upper::TU = nothing # box upper bounds
    buffer::TF = 0.01 # if using bounds, buffer brings x inside bounds by buffer amountd
    Lp::TF = Inf # p value for p-norm for convergence criteria
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

##### ----- Quasi-Newton Solvers ----- #####

"""
    struct MinpackOptions{TB,TF,TI,TSym} <: ExternalSolverOptions

Options for the MINPACK's HYBRJ solver

# Fields
- `algorithm::TSym = :hybr` : algorithm to use in MINPACK.jl (hybr is HYBRJ when the jacobian is provided)
- `atol::TF = 1e-12` : absolute convergence tolerance
- `iteration_limit::TF = 100` : maximum number of iterations
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct MinpackOptions{TB,TF,TI,TSym} <: ExternalSolverOptions
    algorithm::TSym = :hybr
    atol::TF = 1e-12
    iteration_limit::TI = 100
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct SIAMFANLEOptions{TA,TB,TF,TI,TK} <: ExternalSolverOptions

Options for the SIAMFANLEquations pacakge solvers

# Fields
- `algorithm::TA = SIAMFANLEquations.nsoli` : algorithm to use
- `rtol::TF = 0.0` : relative convergence tolerance
- `atol::TF = 1e-10` : absolute convergence tolerance
- `iteration_limit::TF = 1000` : maximum number of iterations
- `linear_iteration_limit::TF = 5` : maximum number of linear solve iterations (GMRES)
- `additional_kwargs::TK = (;)` : any additional keyword arguments for the solver
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct SIAMFANLEOptions{TA,TB,TF,TI,TK} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TA = SIAMFANLEquations.nsoli
    atol::TF = 1e-10
    rtol::TF = 0.0
    iteration_limit::TI = 1000
    linear_iteration_limit::TI = 5
    additional_kwargs::TK = (;)
    # additional_kwargs::TK = (; delta0=1e-3)
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

##### ----- Newton+ Solvers ----- #####
# NOTE: these also have fixed-point options

"""
    struct NonlinearSolveOptions{TA,TB,TF,TI,TT} <: ExternalSolverOptions

Options for the SimpleNonlinearSolve pacakge solvers

# Fields
- `algorithm::TA = SimpleNonlinearSolve.SimpleNewtonRaphson` : algorithm to use
- `additional_kwargs::TK = (;)` : any additional keyword arguments for the solver
- `atol::TF = 1e-12` : absolute convergence tolerance
- `iteration_limit::TF = 25` : maximum number of iterations
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct NonlinearSolveOptions{TA,TB,TF,TI,TK} <: ExternalSolverOptions
    # Algorithm Options
    algorithm::TA = SimpleNonlinearSolve.SimpleNewtonRaphson
    additional_kwargs::TK = (; autodiff=SimpleNonlinearSolve.AutoForwardDiff())
    # Iteration Controls
    atol::TF = 1e-12
    iteration_limit::TI = 25
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct NLsolveOptions{TB,TF,TK,Tls,Tlsk,TSym} <: ExternalSolverOptions

Options for the NLsolve pacakge solvers

# Fields
- `algorithm::TSym = :anderson` : algorithm to use
- `additional_kwargs::TK = (;)` : any additional keyword arguments for the solver
- `atol::TF = 1e-12` : absolute convergence tolerance
- `iteration_limit::TF = 25` : maximum number of iterations
- `linesearch_method::Tls = LineSearches.MoreThuente` : line search method to use
- `linesearch_kwargs::Tlsk = (;)` : any additional lineseach keyword arguments
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct NLsolveOptions{TB,TF,TI,Tls,Tlsk,TSym} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TSym = :anderson
    atol::TF = 1e-12
    iteration_limit::TI = 100
    # line search parameters
    linesearch_method::Tls = LineSearches.MoreThuente
    linesearch_kwargs::Tlsk = (;)
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

##### ----- Poly-Algorithm Solvers ----- #####

"""
    struct CompositeSolverOptions{
        TB,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}
    } <: PolyAlgorithmOptions

Options for Composite Solvers (start with a partial solve of one solve, then finish with another starting where the first left off).

# Fields
- `solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
    ]' : Vector of solver options to use.
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct CompositeSolverOptions{
    TB,TI,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}
} <: PolyAlgorithmOptions
    solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
    ]
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct ChainSolverOptions{TB,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}} <:PolyAlgorithmOptions

Options for Chain Solvers (try one solver, if it doesn't converge, try another)

# Fields
- `solvers::AbstractVector{TS} = [
        NLsolveOptions(; algorithm=:anderson, atol=1e-12),
        MinpackOptions(; atol=1e-12),
        NonlinearSolveOptions(;
            algorithm=SimpleNonlinearSolve.SimpleNewtonRaphson,
            atol=1e-12,
            additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoForwardDiff()),
        ),
    ] : Vector of solver options to use.
- `converged::AbstractVector{TB} = [false]` : flag to track if convergence took place.
"""
@kwdef struct ChainSolverOptions{
    TB,TI,TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}
} <: PolyAlgorithmOptions
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
    iterations::AbstractVector{TI} = [0]
end

"""
    ChainSolverOptions(multipoint)

Convenience function that set's up chain solver options from defaults for a given number of multi-points.

# Arguments
- `multipoint::Vector` : doesn't need to be anything but a vector of the length of multipoints.

# Returns
- `solver_options::ChainSolverOptions` : A ChainSolverOptions object with arrays for the convergence flags in the overall type as well as inside each of the solver options.
"""
function ChainSolverOptions(multipoint)
    lm = length(multipoint)
    return ChainSolverOptions(;
        solvers=[
            NLsolveOptions(;
                algorithm=:anderson,
                atol=1e-12,
                converged=fill(false, lm),
                iterations=zeros(Int, lm),
            ),
            MinpackOptions(;
                atol=1e-12, converged=fill(false, lm), iterations=zeros(Int, lm)
            ),
            NonlinearSolveOptions(;
                algorithm=SimpleNonlinearSolve.SimpleNewtonRaphson,
                atol=1e-12,
                additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoForwardDiff()),
                converged=fill(false, lm),
                iterations=zeros(Int, lm),
            ),
        ],
        converged=fill(false, lm),
        iterations=zeros(Int, lm),
    )
end

#---------------------------------#
#   ELLIPTIC GRID SOLVER TYPES    #
#---------------------------------#
"""
    struct SLORGridSolverOptions{TB,TF,TI} <: GridSolverOptionsType

Options for SLOR (successive line over relaxation) elliptic grid solver.

# Fields
- `iteration_limit::TI = 100` : maximum number of iterations
- `atol::TF = 1e-9` : absolute convergence tolerance
- `converged::AbstractVector{TB} = [false]
"""
@kwdef struct SLORGridSolverOptions{TB,TF,TI} <: GridSolverOptionsType
    iteration_limit::TI = 100
    atol::TF = 1e-9
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

"""
    struct GridSolverOptions{TB,TF,TI,TSym} <: GridSolverOptionsType

Options for Newton elliptic grid solver.

# Fields
- `iteration_limit::TI = 10` : maximum number of iterations
- `atol::TF = 1e-14` : absolute convergence tolerance
- `algorithm::TSym = :newton` : algorithm to use in NLsolve.jl
- `autodiff::TSym = :forward` : differentiation method to use in NLsolve.jl
- `converged::AbstractVector{TB} = [false]
"""
@kwdef struct GridSolverOptions{TB,TF,TI,TSym} <: GridSolverOptionsType
    iteration_limit::TI = 20
    atol::TF = 1e-14
    algorithm::TSym = :newton
    autodiff::TSym = :forward
    converged::AbstractVector{TB} = [false]
    iterations::AbstractVector{TI} = [0]
end

#---------------------------------#
#         OPTION SET TYPES        #
#---------------------------------#

"""
    struct Options{
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

Type containing (nearly) all the available user options.

# Fields
## General Options
- `verbose::TB = false` : flag to print verbose statements
- `silence_warnings::TB = true` : flag to silence warnings
- `multipoint_index::TI = [1]` : holds current index of multi-point solver (no need for user to change this usually)
## Geometry Re-interpolation and generation options` :
- `finterp::Tin = FLOWMath.akima` : interpolation method used for re-paneling bodies
- `autoshiftduct::TB = true` : flag as to whether duct geometry should be shifted based on rotor tip location
- `lu_decomp_flag::TB = false` : flag indicating if panel method LHS matrix factorization was successful
### paneling options
- `itcpshift::TF = 0.05` : factor for internal trailing edge psuedo-panel placement (default is DFDC hard-coded value)
- `axistol::TF = 1e-15` : tolerance for how close the the axis of rotation should be considered on the axis
- `tegaptol::TF = 1e1 * eps()` : tolerance for how large of a trailing edge gap should be considered a gap
## Integration Options
- `integration_options::TIo = IntegrationOptions()` : integration options
## Post-processing Options
- `write_outputs::TBwo = [false]` : Bool for whether to write the outputs of the analysis to an external file (slow)
- `outfile::TSf = ["outputs.jl"]` : External output file name (including path information) for files to write
- `checkoutfileexists::TB = false` : Flag for whether to check if file exists before overwriting
- `output_tuple_name::TSt = ["outs"]` : variable name for named tuple written to out file
## Solving Options
- `grid_solver_options::WS = GridSolverOptions()` : elliptic grid solver options
- `solver_options::TSo = ChainSolverOptions()` : solver options
"""
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

"""
    function DFDC_options(;
        grid_solver_options=SLORGridSolverOptions(),
        solver_options=CSORSolverOptions(),
        kwargs...,
    )

Convenience function to select options used in DFDC.

"""
function DFDC_options(; kwargs...)
    return Options(;
        grid_solver_options=SLORGridSolverOptions(),
        solver_options=CSORSolverOptions(),
        kwargs...,
    )
end

"""
    function DFDC_options(
        multipoint;
        grid_solver_options=SLORGridSolverOptions(),
        solver_options=CSORSolverOptions(),
        kwargs...,
    )

Convenience function to select options used in DFDC and run multipoint analysis.

# Arguments
- `multipoint::Vector` : doesn't need to be anything but a vector of the length of multipoints.
"""
function DFDC_options(
    multipoint, Vconv; write_outputs=nothing, outfile=nothing, output_tuple_name=nothing, kwargs...
)
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

    return Options(;
        grid_solver_options=SLORGridSolverOptions(),
        solver_options=CSORSolverOptions(;
            converged=fill(false, lm),
            iterations=zeros(Int, lm),
            Vconv=Vconv,
        ),
        write_outputs=write_outputs,
        outfile=outfile,
        output_tuple_name=output_tuple_name,
        kwargs...,
    )
end
