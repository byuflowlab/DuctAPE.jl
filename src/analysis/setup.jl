"""
    setup_analysis(
        propulsor::Propulsor,
        options::Options=set_options();
        prepost_container_caching=nothing,
        solve_parameter_caching=nothing,
        solve_container_caching=nothing,
    )

Perform pre-processing and cache setup (as needed) for propuslor analysis.

# Arguments
- `propulsor::Propulsor` : Propulsor input object (see docstring for `Propulsor` type)
- `options::Options=set_options()` : Options object (see `set_options` and related functions)

# Keyword Arguments
- `prepost_container_caching=nothing` : Output of `allocate_prepost_container_cache`
- `solve_parameter_caching=nothing` : Output of `allocate_solve_parameter_container_cache`
- `solve_container_caching=nothing` : Output of `allocate_solve_container_cache`

# Returns

- `problem_dimensions::NamedTuple` : Named Tuple contiaining bookkeeping information (problem dimensions)
- `prepost_containers::NamedTuple` : Named Tuple containing reshaped views into the prepost cache
- `solve_parameter_cache_vector::Vector` : Vector containing the relevant typed cache vector of solve parameters
- `solve_parameter_cache_dims::NamedTuple` : Named Tuple containing dimensions used for reshaping the solve parameter cache
- `A_bb_LU::LinearAlgebra.LU` : The LU factorization of the AIC matrix used in the panel method
- `lu_decomp_flag::Bool` : flag indicating if the LU decomposition was successful
- `airfoils::Matrix{AFType}` : Matrix contiaining the blade element airfoil polar objects
- `idmaps::NamedTuple` : Named Tuple containing bookkeeping information (index mappings)
"""
function setup_analysis(
    propulsor::Propulsor,
    options=set_options();
    prepost_container_caching=nothing,
    solve_parameter_caching=nothing,
    solve_container_caching=nothing,
)
    # - Get type to dispatch caches - #
    TF = promote_type(
        eltype(propulsor.duct_coordinates),
        eltype(propulsor.centerbody_coordinates),
        eltype(propulsor.operating_point.Vinf),
        eltype(propulsor.operating_point.Omega),
        eltype(propulsor.operating_point.rhoinf),
        eltype(propulsor.operating_point.muinf),
        eltype(propulsor.operating_point.asound),
        eltype(propulsor.rotorstator_parameters.B),
        eltype(propulsor.rotorstator_parameters.Rhub),
        eltype(propulsor.rotorstator_parameters.Rtip),
        eltype(propulsor.rotorstator_parameters.rotorzloc),
        eltype(propulsor.rotorstator_parameters.chords),
        eltype(propulsor.rotorstator_parameters.twists),
    )

    # - Get Problem Dimensions - #
    problem_dimensions = get_problem_dimensions(propulsor.paneling_constants)

    ##### ----- SET UP CACHES AS NEEDED ----- #####

    # - Set up Pre- and Post-process Cache - #
    # Allocate Cache
    if isnothing(prepost_container_caching)
        prepost_container_caching = allocate_prepost_container_cache(
            propulsor.paneling_constants
        )
    end

    # unpack the caching
    (; prepost_container_cache, prepost_container_cache_dims) = prepost_container_caching

    # Get correct cached types
    prepost_container_cache_vec = @views PreallocationTools.get_tmp(
        prepost_container_cache, TF(1.0)
    )

    # reset cache
    prepost_container_cache_vec .= 0

    # Reshape Cache
    prepost_containers = withdraw_prepost_container_cache(
        prepost_container_cache_vec, prepost_container_cache_dims
    )

    # - Set up Solver Sensitivity Paramter Cache - #

    # Allocate Cache
    if isnothing(solve_parameter_caching)
        solve_parameter_caching = allocate_solve_parameter_cache(
            options.solver_options, propulsor.paneling_constants
        )
    end

    # unpack caching
    (; solve_parameter_cache, solve_parameter_cache_dims) = solve_parameter_caching

    # get correct cache type
    solve_parameter_cache_vector = @views PreallocationTools.get_tmp(
        solve_parameter_cache, TF(1.0)
    )
    # reset cache
    solve_parameter_cache_vector .= 0

    # reshape cache
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # copy over operating point
    for f in fieldnames(typeof(propulsor.operating_point))
        solve_parameter_tuple.operating_point[f] .= getfield(propulsor.operating_point, f)
    end

    # - Do preprocessutations - #
    if options.verbose
        println("Pre-computing Parameters")
    end

    ##### ----- PERFORM PREPROCESSING COMPUTATIONS ----- #####

    # - Preprocess - #
    A_bb_LU, lu_decomp_flag, airfoils, idmaps, _ = precompute_parameters!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        propulsor,
        prepost_containers,
        problem_dimensions;
        grid_solver_options=options.grid_solver_options,
        integration_options=options.integration_options,
        autoshiftduct=options.autoshiftduct,
        itcpshift=options.itcpshift,
        axistol=options.axistol,
        tegaptol=options.tegaptol,
        finterp=options.finterp,
        silence_warnings=options.silence_warnings,
        verbose=options.verbose,
    )

    return problem_dimensions,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    A_bb_LU,
    lu_decomp_flag,
    airfoils,
    idmaps
end