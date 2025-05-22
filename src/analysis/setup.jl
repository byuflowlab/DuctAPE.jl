"""
    setup_analysis(
        ducted_rotor::DuctedRotor,
        operating_point::OperatingPoint,
        options::Options=set_options();
        prepost_container_caching=nothing,
        solve_parameter_caching=nothing,
        solve_container_caching=nothing,
    )

Perform pre-processing and cache setup (as needed) for propuslor analysis.

# Arguments
- `ducted_rotor::DuctedRotor` : DuctedRotor input object (see docstring for `DuctedRotor` type)
- `operating_point::OperatingPoint` : OperatingPoint input object (see docstring for `OperatingPoint` type)
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
- `idmaps::NamedTuple` : Named Tuple containing bookkeeping information (index mappings)
"""
function setup_analysis(
    ducted_rotor::DuctedRotor,
    operating_point::OperatingPoint,
    options=set_options();
    prepost_container_caching=nothing,
    solve_parameter_caching=nothing,
    solve_container_caching=nothing,
)
    # - Get type to dispatch caches - #
    TF = promote_type(
        eltype(ducted_rotor.duct_coordinates),
        eltype(ducted_rotor.centerbody_coordinates),
        eltype(ducted_rotor.rotor.B),
        eltype(ducted_rotor.rotor.Rhub),
        eltype(ducted_rotor.rotor.Rtip),
        eltype(ducted_rotor.rotor.rotorzloc),
        eltype(ducted_rotor.rotor.chords),
        eltype(ducted_rotor.rotor.twists),
        eltype(operating_point.Vinf),
        eltype(operating_point.Omega),
        eltype(operating_point.rhoinf),
        eltype(operating_point.muinf),
        eltype(operating_point.asound),
    )

    # - Get Problem Dimensions - #
    problem_dimensions = get_problem_dimensions(ducted_rotor.paneling_constants)

    ##### ----- SET UP CACHES AS NEEDED ----- #####

    # - Set up Pre- and Post-process Cache - #
    # Allocate Cache
    if isnothing(prepost_container_caching)
        prepost_container_caching = allocate_prepost_container_cache(
            ducted_rotor.paneling_constants
        )
    else
        # reset cache
        prepost_container_caching.prepost_container_cache.du .= 0
        prepost_container_caching.prepost_container_cache.dual_du .= 0
    end

    # unpack the caching
    (; prepost_container_cache, prepost_container_cache_dims) = prepost_container_caching

    # Get correct cached types
    prepost_container_cache_vector = @views PreallocationTools.get_tmp(
        prepost_container_cache, TF(1.0)
    )

    # reset cache
    prepost_container_cache_vector .= 0

    # Reshape Cache
    prepost_containers = withdraw_prepost_container_cache(
        prepost_container_cache_vector, prepost_container_cache_dims
    )

    # - Set up Solver Sensitivity Paramter Cache - #

    # Allocate Cache
    if isnothing(solve_parameter_caching)
        solve_parameter_caching = allocate_solve_parameter_cache(
            options.solver_options,
            ducted_rotor.paneling_constants,
            ducted_rotor.rotor.airfoils,
        )
    else
        # reset cache
        solve_parameter_caching.solve_parameter_cache.du .= 0
        solve_parameter_caching.solve_parameter_cache.dual_du .= 0
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
    for f in fieldnames(typeof(operating_point))
        if f != :units
            solve_parameter_tuple.operating_point[f] .= getfield(operating_point, f)
        end
    end

    ##### ----- PERFORM PREPROCESSING COMPUTATIONS ----- #####
    if options.verbose
        println("Pre-computing Parameters")
    end

    # - Preprocess - #
    A_bb_LU, lu_decomp_flag, idmaps, _ = precompute_parameters!(
        solve_parameter_tuple.ivr,
        solve_parameter_tuple.ivw,
        solve_parameter_tuple.blade_elements,
        solve_parameter_tuple.linsys,
        solve_parameter_tuple.wakeK,
        ducted_rotor,
        operating_point,
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
    idmaps
end
