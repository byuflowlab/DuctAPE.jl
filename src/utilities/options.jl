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
    abstract type InternalSolverOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type InternalSolverOptions <: SolverOptionsType end

"""
    abstract type ExternalSolverOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type ExternalSolverOptions <: SolverOptionsType end

"""
    abstract type InternalPolyAlgorithmOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type InternalPolyAlgorithmOptions <: InternalSolverOptions end

"""
    abstract type ExternalPolyAlgorithmOptions <: SolverOptionsType

Used for solver dispatch.
"""
abstract type ExternalPolyAlgorithmOptions <: ExternalSolverOptions end

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

"""
    abstract type BoundaryLayerOptions

Used in boundary layer method dispatch
"""
abstract type BoundaryLayerOptions end

"""
    abstract type RK

Used for selecting boundary layer ODE solver
"""
struct RK end

"""
    abstract type DiffEq

Used for selecting boundary layer ODE solver
"""
struct DiffEq end

#---------------------------------#
#         QUADRATURE TYPES        #
#---------------------------------#

"""
    struct Romberg <: IntegrationMethod

Options for Romberg integration method

# Fields
- `max_subdivisions::Int = 10` : maximum number of subdivisions. Note, total number of internvals is 2^N, where N is number of subdivisions.
- `atol::Float = 1e-6` : absolute error tolerance.
"""
@kwdef struct Romberg{TF,TI} <: IntegrationMethod
    max_subdivisions::TI = 10
    atol::TF = 1e-6
end

"""
    struct GaussKronrod <: IntegrationMethod

Options for Gauss-Kronrod integration method

# Fields
- `order::Int = 7` : order of Legendre polynomial to use on each interval
- `maxevales::Int = 10^7` : maximum number of evaluations in the adaptive method
- `atol::Float = 0.0` : absolute error tolerance. (note, if zero, QuadGK uses sqrt(eps()) relative tolerance).
"""
@kwdef struct GaussKronrod{TF,TI} <: IntegrationMethod
    order::TI = 7
    maxevals::TI = 10^7
    atol::TF = 0.0
end

"""
    struct GaussLegendre <: IntegrationMethod

Options for Gauss-Legendre integration method

# Fields
- `sample_points::Int` : Sample Points
- `weights::Int` : Gauss weights
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
    struct IntegrationOptions

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
    struct CSORSolverOptions <: InternalSolverOptions

Type containing all the options for the CSOR (controlled successive over relaxation) solver.

Note that the defaults match DFDC settings.

# Fields
- `verbose::Bool = false` : flag to print verbose statements
- `iteration_limit::Float = 1e2` : maximum number of iterations
- `nrf::Float = 0.4` : nominal relaxation factor
- `bt1::Float = 0.2` : backtracking factor 1
- `bt2::Float = 0.6` : backtracking factor 2
- `pf1::Float = 0.4` : press forward factor 1
- `pf2::Float = 0.5` : press forward factor 2
- `btw::Float = 0.6` : backtracking factor for wake
- `pfw::Float = 1.2` : press forward factor for wake
- `f_circ::Float = 1e-3` : convergence tolerance for rotor circulation
- `f_dgamw::Float = 2e-4` : convergence tolerance for wake vortex strength
- `convergence_type::ConvergenceType = Relative()` : dispatch for relative or absolute convergence criteria.
- `Vconv::AbstractArray{Float} = [1.0]` : velocity used in relative convergence criteria (should be set to Vref).
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct CSORSolverOptions{TB,TC<:ConvergenceType,TF,TI} <: InternalSolverOptions
    # Defaults are DFDC hard-coded values
    verbose::TB = false
    iteration_limit::TI = 1000
    nrf::TF = 0.4
    bt1::TF = 0.2
    bt2::TF = 0.6
    pf1::TF = 0.4
    pf2::TF = 0.5
    btw::TF = 0.6
    pfw::TF = 1.2
    # relaxation_schedule::TS = [
    #     reverse!([1e10; 1e-16; 1e-17; 0.0]), reverse!([0.0; 0.0; 1.0; 1.0])
    # ]
    f_circ::TF = 1e-3
    f_dgamw::TF = 2e-4
    convergence_type::TC = Relative()
    Vconv::AbstractArray{TF} = [1.0]
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    CSORSolverOptions(multipoint; kwargs...)

Convenience function that sets up CSOR solver options from defaults for a given number of multi-points.

# Arguments
- `multipoint::Vector` : doesn't need to be anything but a vector of the length of multipoints.

# Returns
- `solver_options::CSORSolverOptions` : A CSORSolverOptions object with arrays for the convergence flags in the overall type as well as inside each of the solver options.
"""
function CSORSolverOptions(multipoint; kwargs...)
    lm = length(multipoint)

    return CSORSolverOptions(;
        converged=fill(false, (1, lm)),
        iterations=zeros(Int, (1, lm)),
        residuals=-ones(1, lm),
        kwargs...,
    )
end

"""
    struct ModCSORSolverOptions <: InternalSolverOptions

Type containing all the options for the modified CSOR solver.

# Fields
- `verbose::Bool = false` : flag to print verbose statements
- `iteration_limit::Float = 350` : maximum number of iterations
- `relaxation_parameters::NamedTuple` = (;
    - `nrf::Float = 0.4` : nominal relaxation factor
    - `bt1::Float = 0.2` : backtracking factor 1
    - `bt2::Float = 0.6` : backtracking factor 2
    - `pf1::Float = 0.4` : press forward factor 1
    - `pf2::Float = 0.5` : press forward factor 2
    - `btw::Float = 0.6` : backtracking factor for wake
    - `pfw::Float = 1.2` : press forward factor for wake
    ) : parameters for determining relaxation level of states in each iteration.
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct ModCSORSolverOptions{TB,TF,TI,TT} <: InternalSolverOptions
    # Defaults are DFDC hard-coded values
    verbose::TB = false
    iteration_limit::TI = 500
    relaxation_parameters::TT = (;
        nrf=0.4, bt1=0.2, bt2=0.6, pf1=0.4, pf2=0.5, btw=0.6, pfw=1.2
    )
    atol::TF = 5e-10
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    ModCSORSolverOptions(multipoint; kwargs...)

Convenience function that sets up CSOR solver options from defaults for a given number of multi-points.

# Arguments
- `multipoint::Vector` : doesn't need to be anything but a vector of the length of multipoints.

# Returns
- `solver_options::ModCSORSolverOptions` : A ModCSORSolverOptions object with arrays for the convergence flags in the overall type as well as inside each of the solver options.
"""
function ModCSORSolverOptions(multipoint; kwargs...)
    lm = length(multipoint)

    return ModCSORSolverOptions(;
        converged=fill(false, (1, lm)),
        iterations=zeros(Int, (1, lm)),
        residuals=-ones(1, lm),
        kwargs...,
    )
end

"""
    struct FixedPointOptions <: ExternalSolverOptions

Options for the FixedPoint.jl package solver

# Fields
- `iteration_limit::Int = 1000` : maximum number of iterations
- `vel::Float = 0.9` : vel keyword argument, default is package default
- `ep::Float = 0.01` : ep keyword argument, default is package default
- `atol::Float = 1e-12` : absolute convergence tolerance
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct FixedPointOptions{TB,TF,TI} <: ExternalSolverOptions
    iteration_limit::TI = 1000
    vel::TF = 0.9
    ep::TF = 0.01
    atol::TF = 1e-10
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    struct SpeedMappingOptions <: ExternalSolverOptions

Options for the SpeedMapping.jl package solver

# Fields
- `orders::AbstractArray{Int} = [3, 2]`
- `sig_min::Int = 0` : maybe set to 1?
- `stabilize::Bool = false` : stabilizes before extrapolation
- `check_obj::Bool = false` : checks for inf's and nan's and starts from previous finite point
- `atol::Float = 1e-10` : absolute convergence tolerance
- `iteration_limit::Float = 1000` : maximum number of iterations
- `time_limit::Float = Inf` : time limit in seconds
- `lower::Float = nothing` : box lower bounds
- `upper::Float = nothing` : box upper bounds
- `buffer::Float = 0.01` : if using bounds, buffer brings x inside bounds by buffer amountd
- `Lp::Float = Inf` : p value for p-norm for convergence criteria
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct SpeedMappingOptions{TB,TF,TI,TL,TSm,TU} <: ExternalSolverOptions
    orders::AbstractArray{TI} = [3, 2]
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
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

##### ----- Quasi-Newton Solvers ----- #####

"""
    struct MinpackOptions <: ExternalSolverOptions

Options for the MINPACK's HYBRJ solver

# Fields
- `algorithm::Symbol = :hybr` : algorithm to use in MINPACK.jl (hybr is HYBRJ when the jacobian is provided)
- `atol::FLoat = 1e-10` : absolute convergence tolerance
- `iteration_limit::FLoat = 100` : maximum number of iterations
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct MinpackOptions{TB,TF,TI,TSym} <: ExternalSolverOptions
    algorithm::TSym = :hybr
    atol::TF = 1e-10
    iteration_limit::TI = 100
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    struct SIAMFANLEOptions <: ExternalSolverOptions

Options for the SIAMFANLEquations pacakge solvers

# Fields
- `algorithm::SIAMFANLEquations algorithm = SIAMFANLEquations.nsoli` : algorithm to use
- `rtol::Float = 0.0` : relative convergence tolerance
- `atol::Float = 1e-10` : absolute convergence tolerance
- `iteration_limit::Int = 1000` : maximum number of iterations
- `linear_iteration_limit::Float = 5` : maximum number of linear solve iterations (GMRES)
- `additional_kwargs = (;)` : any additional keyword arguments for the solver
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct SIAMFANLEOptions{TA,TB,TF,TI,TK} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TA = SIAMFANLEquations.nsoli
    atol::TF = 1e-10
    rtol::TF = 0.0
    iteration_limit::TI = 1000
    linear_iteration_limit::TI = 5
    additional_kwargs::TK = (;)
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

##### ----- Newton+ Solvers ----- #####
# NOTE: these also have fixed-point options

"""
    struct NonlinearSolveOptions <: ExternalSolverOptions

Options for the SimpleNonlinearSolve pacakge solvers

# Fields
- `algorithm::SimpleNonlinearSolve algorithm = SimpleNonlinearSolve.SimpleNewtonRaphson` : algorithm to use
- `additional_kwargs = (;)` : any additional keyword arguments for the solver
- `atol::Float = 1e-12` : absolute convergence tolerance
- `iteration_limit::Float = 25` : maximum number of iterations
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct NonlinearSolveOptions{TA,TB,TF,TI,TK} <: ExternalSolverOptions
    # Algorithm Options
    algorithm::TA = SimpleNonlinearSolve.SimpleNewtonRaphson
    additional_kwargs::TK = (; autodiff=SimpleNonlinearSolve.AutoForwardDiff())
    # Iteration Controls
    atol::TF = 1e-10
    iteration_limit::TI = 25
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    struct NLsolveOptions <: ExternalSolverOptions

Options for the NLsolve pacakge solvers

# Fields
- `algorithm::Symbol = :anderson` : algorithm to use
- `additional_kwargs = (;)` : any additional keyword arguments for the solver
- `atol::Float = 1e-12` : absolute convergence tolerance
- `iteration_limit::Int = 25` : maximum number of iterations
- `linesearch_method::LineSearches method = LineSearches.MoreThuente` : line search method to use
- `linesearch_kwargs = (;)` : any additional lineseach keyword arguments
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct NLsolveOptions{TB,TF,TI,Tls,Tlsk,TSym} <: ExternalSolverOptions
    # Options for overall solve
    algorithm::TSym = :anderson
    atol::TF = 1e-10
    iteration_limit::TI = 100
    # line search parameters
    linesearch_method::Tls = LineSearches.MoreThuente
    linesearch_kwargs::Tlsk = (;)
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

##### ----- Poly-Algorithm Solvers ----- #####

"""
    struct CompositeSolverOptions <: ExternalPolyAlgorithmOptions

Options for Composite Solvers (start with a partial solve of one solve, then finish with another starting where the first left off).

# Fields
- `solvers::AbstractArray{SolverOptionsType} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-10),
    ]` : Vector of solver options to use.
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct CompositeSolverOptions{
    TB,TI,TF,TS<:Union{ExternalSolverOptions,ExternalPolyAlgorithmOptions}
} <: ExternalPolyAlgorithmOptions
    solvers::AbstractArray{TS} = [
        NLsolveOptions(; algorithm=:newton, iteration_limit=3),
        NLsolveOptions(; algorithm=:anderson, atol=1e-10),
    ]
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    struct ChainSolverOptions <:ExternalPolyAlgorithmOptions

Options for Chain Solvers (try one solver, if it doesn't converge, try another)

# Fields
- `solvers::AbstractArray{SolverOptionsType} = [
        NLsolveOptions(; algorithm=:anderson, atol=1e-10),
        MinpackOptions(; atol=1e-10),
        NonlinearSolveOptions(;
            algorithm=SimpleNonlinearSolve.SimpleNewtonRaphson,
            atol=1e-12,
            additional_kwargs=(; autodiff=SimpleNonlinearSolve.AutoForwardDiff()),
        ),
    ]` : Vector of solver options to use.
- `converged::AbstractArray{Bool} = [false]` : flag to track if convergence took place.
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residuals::AbstractArray{Float} = [-1.0]` : final residual values
"""
@kwdef struct ChainSolverOptions{TB,TI,TF,TS<:ExternalSolverOptions} <:
              ExternalPolyAlgorithmOptions
    solvers::AbstractArray{TS} = [
        NLsolveOptions(; algorithm=:anderson, atol=1e-10, iteration_limit=200),
        NLsolveOptions(; algorithm=:trust_region, atol=1e-10, iteration_limit=100),
        MinpackOptions(; atol=1e-10, iteration_limit=100),
    ]
    converged::AbstractArray{TB} = [false, false, false]
    iterations::AbstractArray{TI} = [0, 0, 0]
    residuals::AbstractArray{TF} = [-1.0]
end

"""
    ChainSolverOptions(multipoint)

Convenience function that set's up chain solver options from defaults for a given number of multi-points.

# Arguments
- `multipoint::Vector` : doesn't need to be anything but a vector of the length of multipoints.

# Returns
- `solver_options::ChainSolverOptions` : A ChainSolverOptions object with arrays for the convergence flags in the overall type as well as inside each of the solver options.
"""
function ChainSolverOptions(multipoint; solvers=nothing)
    lm = length(multipoint)

    if isnothing(solvers)
        solvers = [
            NLsolveOptions(;
                algorithm=:anderson,
                atol=1e-10,
                converged=fill(false, lm),
                iterations=zeros(Int, lm),
            ),
            NLsolveOptions(;
                algorithm=:trust_region,
                atol=1e-10,
                converged=fill(false, lm),
                iterations=zeros(Int, lm),
            ),
            MinpackOptions(;
                atol=1e-10, converged=fill(false, lm), iterations=zeros(Int, lm)
            ),
        ]
    end

    return ChainSolverOptions(;
        solvers=solvers,
        converged=fill(false, (length(solvers), lm)),
        iterations=zeros(Int, (length(solvers), lm)),
        residuals=-ones(1, lm),
    )
end

#---------------------------------#
#   ELLIPTIC GRID SOLVER TYPES    #
#---------------------------------#
"""
    struct SLORGridSolverOptions <: GridSolverOptionsType

Options for SLOR (successive line over relaxation) elliptic grid solver.

# Fields
- `iteration_limit::Int = 100` : maximum number of iterations
- `atol::Float = 1e-9` : absolute convergence tolerance
- `converged::AbstractArray{Bool}` = [false]
- `iterations::AbstractArray{Int} = [0]` : iteration counter
"""
@kwdef struct SLORGridSolverOptions{TB,TF,TI} <: GridSolverOptionsType
    iteration_limit::TI = 200
    atol::TF = eps()
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
end

"""
    struct GridSolverOptions <: GridSolverOptionsType

Options for Newton elliptic grid solver.

# Fields
- `iteration_limit::Int = 20` : maximum number of iterations
- `atol::Float = 3e-10` : absolute convergence tolerance
- `algorithm::Symbol = :newton` : algorithm to use in NLsolve.jl
- `autodiff::Symbol = :forward` : differentiation method to use in NLsolve.jl
- `precondition = false` : flag to precondition with SLOR
- `precondition_max_iterations = 3` : number of precondition iterations
- `converged::AbstractArray{Bool}` = [false]
- `iterations::AbstractArray{Int} = [0]` : iteration counter
- `residual_value::AbstractArray{Int} = [0]` : residual value
"""
@kwdef struct GridSolverOptions{TB,TF,TI,TSym} <: GridSolverOptionsType
    iteration_limit::TI = 30
    atol::TF = 3e-10
    algorithm::TSym = :newton
    autodiff::TSym = :forward
    precondition::TB = false
    precondition_max_iterations::TI = 3
    converged::AbstractArray{TB} = [false]
    iterations::AbstractArray{TI} = [0]
    residual_value::AbstractArray{TF} = [0.0]
end

#---------------------------------#
#       BOUNDARY LAYER TYPES      #
#---------------------------------#
"""
    struct HeadsBoundaryLayerOptions

# Fields:
- `model_drag::Bool=false` : flag to turn on viscous drag approximation
- `termiante::Bool=true` : flag to terminate solver when separation criteria is met.
- `return_last_max_shape_factor::Bool=true` : return the last maximum shape factor to avoid false drops near the trailing edge.
- `cutoff_Hsep::Bool=false` : Cutoff returned Hsep vector at separation criteria
- `n_steps::Int = Int(5e3)` : number of steps to use in boundary layer integration
- `first_step_size::Float = 1e-6` : size of first step in boundary layer integration
- `upper_step_size::Float=nothing` : uses fixed step size rather than total number of steps.
- `lower_step_size::Float=nothing` : uses fixed step size rather than total number of steps.
- `offset::Float = 1e-3` : size of offset for (where to initialize) boundary layer integration
- `solver_type::AbstractODESolverType=RK()` : type of ODE solver (RK() or DiffEq())
- `ode::Function=RK2` : solver to use for boundary layer integration (RadauIIA5, RK4, or RK2 available)
- `separation_criteria::Float=3.0` : value of H12 after which separation should happen.
- `separation_allowance_upper::Int=10` : upper side allowance for how many steps ahead of the trailing edge we'll allow separation without penalty
- `separation_allowance_lower::Int=10` : lower side allowance for how many steps ahead of the trailing edge we'll allow separation without penalty
- `separation_penalty_upper::Float=0.2` : upper side maximum penalty value for separation (at leading edge)
- `separation_penalty_lower::Float=0.2` : lower side maximum penalty value for separation (at leading edge)
- `apply_separation_penalty_to_rotor::Bool=false` : flag to apply separation penalty to rotor performance.
- `dy_eps::Float=0.0` : temporary development parameter.
- `H1_eps::Float=1e-4`: temporary development parameter.
- `H_eps::Float=0.0`: temporary development parameter.
- `verbose::Bool=false` : flag to print verbose statements each iterations (beware; it's a lot)
"""
@kwdef struct HeadsBoundaryLayerOptions{Tb,Te,Tf,Tfun,Ti,To,Tp,Ts,Tsol,Tssl,Tssu} <:
              BoundaryLayerOptions
    model_drag::Tb = false
    terminate::Tb = true
    return_last_max_shape_factor::Tb = true
    cutoff_Hsep::Tb = false
    n_steps::Ti = Int(5e3)
    first_step_size::Tf = 1e-6
    upper_step_size::Tssu = nothing
    lower_step_size::Tssl = nothing
    offset::To = 1e-3
    solver_type::Tsol = RK()
    ode::Tfun = RK2
    separation_criteria::Ts = 3.0
    separation_allowance_upper::Ti = 10
    separation_allowance_lower::Ti = 10
    separation_penalty_upper::Tp = 0.2
    separation_penalty_lower::Tp = 0.2
    apply_separation_penalty_to_rotor::Tb = false
    dy_eps::Te = 0.0
    H1_eps::Te = 1e-4
    H_eps::Te = 0.0
    verbose::Tb = false
end

# """
#     struct GreensBoundaryLayerOptions
#
# NOTE: Green's method is mostly implemented, but there are several bugs still, especially when using Imperial units.
# Known Bugs:
# - Imperial units overestimate momentum thickness.  Likely a unit conversion bug.
# - In some cases of non-separation, the momentum thickens or shape parameter becomes exceedingly large, vastly overestimating the drag coefficient.  Likely the product of one or more of the adjustments to try and make the method more robust.
#
# # Fields:
# - `model_drag::Tb=true` : flag to turn off viscous drag approximation
# - `lambda::Bool = true` : flag to add secondary influences into boundary layer residuals
# - `longitudinal_curvature::Bool = true` : if `lambda`=true, flag to add longitudinal curvature influence into boundary layer residuals
# - `lateral_strain::Bool = true` : if `lambda`=true, flag to add lateral strain influence into boundary layer residuals
# - `dilation::Bool = true` : if `lambda`=true, flag to add dilation influence into boundary layer residuals
# - `n_steps::Int = Int(2e2)` : number of steps to use in boundary layer integration
# - `first_step_size::Float = 1e-3` : size of first step in boundary layer integration
# - `offset::Float = 1e-2` : size of offset for (where to initialize) boundary layer integration
# - `solver_type::Type = DiffEq()` : type of ODE solver (RK() or DiffEq())
# - `ode::Function = RadauIIA5` : solver to use for boundary layer integration (RadauIIA5, RK4, or RK2 available)
# - `separation_allowance_upper::Int=3` : upper side allowance for how many steps ahead of the trailing edge we'll allow separation without penalty
# - `separation_allowance_lower::Int=3` : lower side allowance for how many steps ahead of the trailing edge we'll allow separation without penalty
# - `separation_penalty_upper::Float=0.2` : upper side maximum penalty value for separation (at leading edge)
# - `separation_penalty_lower::Float=0.2` : lower side maximum penalty value for separation (at leading edge)
# """
# @kwdef struct GreensBoundaryLayerOptions{Tb,Tf,Tfun,Ti,To,Tp,Ts,Tsol,Tssl,Tssu} <:
#               BoundaryLayerOptions
#     model_drag::Tb = true
#     terminate::Tb = true
#     lambda::Tb = false
#     longitudinal_curvature::Tb = true
#     lateral_strain::Tb = true
#     dilation::Tb = true
#     n_steps::Ti = Int(2e2)
#     first_step_size::Tf = 1e-6
#     upper_step_size::Tssu = nothing
#     lower_step_size::Tssl = nothing
#     offset::To = 1e-3
#     solver_type::Tsol = DiffEq()
#     ode::Tfun = RadauIIA5
#     separation_criteria::Ts = 0.0
#     separation_allowance_upper::Ti = 25
#     separation_allowance_lower::Ti = 25
#     separation_penalty_upper::Tp = 0.2
#     separation_penalty_lower::Tp = 0.2
#     verbose::Tb = false
# end

#---------------------------------#
#         OPTION SET TYPES        #
#---------------------------------#

"""
    struct Options

Type containing (nearly) all the available user options.

# Fields
## General Options
- `verbose::Bool = false` : flag to print verbose statements
- `silence_warnings::Bool = true` : flag to silence warnings
- `hard_fail::Bool = true` : flag as to whether DuctAPE should return nothing immediately after a failed initialization of the elliptic grid or a failed decomposition of the body influence matrix.  If set to false, DuctAPE will attempt to return objects of the correct size, but with initialized values only.
## Pre-processing Options
### Geometry Interpolation and Generation Options
- `finterp::Interplation Method = FLOWMath.akima` : interpolation method used for re-paneling bodies
- `autoshiftduct::Bool = true` : flag as to whether duct geometry should be shifted based on rotor tip location
### Paneling Options
- `itcpshift::Float = 0.05` : factor for internal trailing edge psuedo-panel placement (default is DFDC hard-coded value)
- `axistol::Float = 1e-15` : tolerance for how close the the axis of rotation should be considered on the axis
- `tegaptol::Float = 1e1 * eps()` : tolerance for how large of a trailing edge gap should be considered a gap
## Integration Options
- `integration_options::IntegrationOptions type = IntegrationOptions()` : integration options
## Post-processing Options
- `boundary_layer_options::BoundaryLayerOptions` : BoundaryLayerOptions object
- `write_outputs::AbstractArray{Bool} = [false]` : Bool for whether to write the outputs of the analysis to an external file (slow)
- `outfile::AbstractArray{String} = ["outputs.jl"]` : External output file name (including path information) for files to write
- `checkoutfileexists::Bool = false` : Flag for whether to check if file exists before overwriting
- `output_tuple_name::AbstractArray{String} = ["outs"]` : variable name for named tuple written to out file
## Solving Options
- `grid_solver_options::GridSolverOptionsType = GridSolverOptions()` : elliptic grid solver options
- `solver_options::SolverOptionsType = ChainSolverOptions()` : solver options
## Bookkeeping Options
- `multipoint_index::Int = [1]` : holds current index of multi-point solver (no need for user to change this usually)
- `lu_decomp_flag::Bool = false` : flag indicating if panel method LHS matrix factorization was successful
"""
@kwdef struct Options{
    TB,
    TBL,
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
    finterp::Tin = (x, y, xp) -> FLOWMath.akima(x, y, xp, 2.0 * eps(), eps())
    autoshiftduct::TB = true
    lu_decomp_flag::TB = false
    # paneling options
    itcpshift::TF = 0.05
    axistol::TF = 1e-15
    tegaptol::TF = 1e1 * eps()
    # - Integration Options - #
    integration_options::TIo = IntegrationOptions()
    # - Post-processing Options - #
    boundary_layer_options::TBL = HeadsBoundaryLayerOptions()
    write_outputs::TBwo = [false]
    outfile::TSf = ["outputs.jl"]
    checkoutfileexists::TB = false
    output_tuple_name::TSt = ["outs"]
    # - Solving Options - #
    grid_solver_options::WS = GridSolverOptions()
    solver_options::TSo = ModCSORSolverOptions()
    # - Failure Options - #
    hard_fail::TB = true
end

"""
    set_options(; kwargs...)
    set_options(multipoint; kwargs...)

Set the options for DuctAPE to use.

Note that the vast majority of the available options are defined through keyword arguments.  See the documentation for the various option types for more information.

# Arguments
- `multipoint::AbstractArray{OperatingPoint}` : a vector of operating points to use if running a multi-point analysis.
"""
function set_options(; kwargs...)
    return Options(; kwargs...)
end

function set_options(
    multipoint::AbstractArray{TM};
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
        solver_options=ModCSORSolverOptions(multipoint),
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
    multipoint,
    Vconv;
    write_outputs=nothing,
    outfile=nothing,
    output_tuple_name=nothing,
    kwargs...,
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
        solver_options=CSORSolverOptions(multipoint; Vconv=Vconv),
        write_outputs=write_outputs,
        outfile=outfile,
        output_tuple_name=output_tuple_name,
        kwargs...,
    )
end
