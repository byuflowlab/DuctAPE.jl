"""
    analyze(
        propulsor::Propulsor,
        options::Options=set_options();
        prepost_container_caching=nothing,
        solve_parameter_caching=nothing,
        solve_container_caching=nothing,
        return_inputs=false,
    )

Analyze propulsor, including preprocessing.

# Arguments
- `propulsor::Propulsor` : Propulsor input object (see docstring for `Propulsor` type)
- `options::Options=set_options()` : Options object (see `set_options` and related functions)

# Keyword Arguments
- `prepost_container_caching=nothing` : Output of `allocate_prepost_container_cache`
- `solve_parameter_caching=nothing` : Output of `allocate_solve_parameter_container_cache`
- `solve_container_caching=nothing` : Output of `allocate_solve_container_cache`
- `return_inputs=false` : flag as to whether or not to return the pre-processed inputs

# Returns
- `outs::NamedTuple` : Named Tuple of various analysis outputs (see docstring for postprocess for details), note, if linear system decomposition fails, no solve is performed and an empty vector is returned.
- `ins::NamedTuple` : Named Tuple of various pre-processed inputs (e.g. panels and body linear system), will only be returned if `return_inputs=true`
- `convergence_flag` : Flag for successful solve convergence
"""
function analyze(
    propulsor::Propulsor,
    options::Options=set_options();
    prepost_container_caching=nothing,
    solve_parameter_caching=nothing,
    solve_container_caching=nothing,
    return_inputs=false,
)

    # - Set Up - #
    problem_dimensions, prepost_containers, solve_parameter_cache_vector, solve_parameter_cache_dims, A_bb_LU, lu_decomp_flag, airfoils, idmaps = setup_analysis(
        propulsor,
        options;
        prepost_container_caching=prepost_container_caching,
        solve_parameter_caching=solve_parameter_caching,
        solve_container_caching=solve_container_caching,
    )

    # - Check that the precomputation went well - #
    #=
      NOTE: If the linear system or wake did not converge, there is likely a serious problem that would lead to an error in the solve, so we will exit here with a fail flag for an optimizer or user
    =#
    if iszero(lu_decomp_flag) || !options.grid_solver_options.converged[1]
        if !options.silence_warnings
            if iszero(lu_decomp_flag)
                @warn "Exiting.  LU decomposition of the LHS matrix for the linear system failed.  Please check your body geometry and ensure that there will be no panels lying directly atop eachother or other similar problematic geometry."
            elseif !options.grid_solver_options.converged[1]
                @warn "Exiting. Wake elliptic grid solve did not converge. Consider a looser convergence tolerance if the geometry looks good."
            end
        end
        #TODO: write a function that returns the same as outs below, but all zeros
        #TODO: probably just call  the post-process function directly and return a reset_container! of the output
        if return_inputs
            return [],#zero_outputs(),
            (; solve_parameter_tuple..., airfoils, idmaps, panels, problem_dimensions),
            false
        else
            return [],#zero_outputs(),
            false
        end
    end

    # - Continue with Analysis - #
    return analyze(
        propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=solve_container_caching,
        return_inputs=return_inputs,
    )
end

"""
    analyze(
        propulsor::Propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options::Options=set_options();
        return_inputs=false,
        solve_container_caching=nothing,
    )

Analyze propulsor, assuming `setup_analysis` has been called and the outputs thereof are being passed in here.

# Arguments
- `propulsor::Propulsor` : Propulsor input object
- `prepost_containers::NamedTuple` : An output from `setup_analysis` containing reshaped views into the prepost cache
- `solve_parameter_cache_vector::Vector` : An output from `setup_analysis` containing the relevant typed cache vector of solve parameters
- `solve_parameter_cache_dims::NamedTuple` : An output from `setup_analysis` containing dimensions used for reshaping the solve parameter cache
- `airfoils::Vector{AFType}` : An output from `setup_analysis` contiaining the blade element airfoil polar objects
- `A_bb_LU::LinearAlgebra.LU` : An output from `setup_analysis` that is the LU decomposition of the AIC matrix used in the panel method
- `idmaps::NamedTuple` : An output from `setup_analysis` containing bookkeeping information (index mappings)
- `problem_dimensions::NamedTuple` : An output from `setup_analysis` contiaining bookkeeping information (problem dimensions)
- `options::Options=set_options()` : Options object

# Keyword Arguments
- `solve_container_caching=nothing` : Output of `allocate_solve_container_cache`
- `return_inputs=false` : flag as to whether or not to return the pre-processed inputs

# Returns
- `outs::NamedTuple` : Named Tuple of various analysis outputs (see docstring for postprocess for details), note, if linear system decomposition fails, no solve is performed and an empty vector is returned.
- `ins::NamedTuple` : Named Tuple of various pre-processed inputs (e.g. panels and body linear system), will only be returned if `return_inputs=true`
- `convergence_flag` : Flag for successful solve convergence
"""
function analyze(
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    idmaps,
    problem_dimensions,
    options::Options=set_options();
    return_inputs=false,
    solve_container_caching=nothing,
)

    # Set up Solve Container Cache
    if isnothing(solve_container_caching)
        solve_container_caching = allocate_solve_container_cache(
            options.solver_options, propulsor.paneling_constants
        )
    end

    # - Process - #
    converged_states = process(
        options.solver_options,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        solve_container_caching,
        idmaps,
        options,
    )

    # - Post-Process - #
    outs = post_process(
        options.solver_options,
        converged_states,
        prepost_containers,
        solve_container_caching,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        propulsor.operating_point,
        propulsor.reference_parameters,
        A_bb_LU,
        airfoils,
        idmaps,
        problem_dimensions;
        write_outputs=options.write_outputs[options.multipoint_index[]],
        outfile=options.outfile[options.multipoint_index[]],
        checkoutfileexists=options.checkoutfileexists,
        output_tuple_name=options.output_tuple_name[options.multipoint_index[]],
        verbose=options.verbose,
    )

    if return_inputs
        solve_parameter_tuple = withdraw_solve_parameter_cache(
            options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
        )

        return outs,
        (;
            prepost_containers.panels,
            prepost_containers.ivb,
            solve_parameter_tuple...,
            blade_elements=(; solve_parameter_tuple.blade_elements..., airfoils...),
            linsys=(; solve_parameter_tuple.linsys..., A_bb_LU),
            idmaps,
            problem_dimensions,
        ),
        options.solver_options.converged[]
    else
        return outs, options.solver_options.converged[options.multipoint_index[]]
    end
end

"""
    analyze(
        multipoint::AbstractVector{OperatingPoint},
        propulsor::Propulsor,
        options::Options=set_options();
        prepost_container_caching=nothing,
        solve_parameter_caching=nothing,
        solve_container_caching=nothing,
        return_inputs=false,
    )

Analyze propulsor, including preprocessing, for a set of operating points.

# Arguments
- `multipoint::AbstractVector{OperatingPoint}` : Vector of Operating Points at which to analyze the propulsor (note that the operating point within the propulsor input will be overwritten with these)
- `propulsor::Propulsor` : Propulsor input object
- `options::Options=set_options()` : Options object

# Keyword Arguments
- `prepost_container_caching=nothing` : Output of `allocate_prepost_container_cache`
- `solve_parameter_caching=nothing` : Output of `allocate_solve_parameter_container_cache`
- `solve_container_caching=nothing` : Output of `allocate_solve_container_cache`
- `return_inputs=false` : flag as to whether or not to return the pre-processed inputs

# Returns
- `outs::Vector{NamedTuple}` : Vector of named tuples of various analysis outputs (see docstring for postprocess for details), note, if linear system decomposition fails, no solve is performed and an empty vector is returned.
- `ins::NamedTuple` : Named Tuple of various pre-processed inputs (e.g. panels and body linear system), will only be returned if `return_inputs=true`
- `convergence_flag` : Flag for successful solve convergence
"""
function analyze(
    multipoint::AbstractVector{TO},
    propulsor::Propulsor,
    options::Options=set_options(multipoint);
    prepost_container_caching=nothing,
    solve_parameter_caching=nothing,
    solve_container_caching=nothing,
    return_inputs=false,
) where {TO<:OperatingPoint}

    # - Set Up - #
    problem_dimensions, prepost_containers, solve_parameter_cache_vector, solve_parameter_cache_dims, A_bb_LU, lu_decomp_flag, airfoils, idmaps = setup_analysis(
        propulsor,
        options;
        prepost_container_caching=prepost_container_caching,
        solve_parameter_caching=solve_parameter_caching,
        solve_container_caching=solve_container_caching,
    )

    # - Check that the precomputation went well - #
    #=
      NOTE: If the linear system or wake did not converge, there is likely a serious problem that would lead to an error in the solve, so we will exit here with a fail flag for an optimizer or user
    =#
    if iszero(lu_decomp_flag) || !options.grid_solver_options.converged[1]
        if !options.silence_warnings
            if iszero(lu_decomp_flag)
                @warn "Exiting.  LU decomposition of the LHS matrix for the linear system failed.  Please check your body geometry and ensure that there will be no panels lying directly atop eachother or other similar problematic geometry."
            elseif !options.grid_solver_options.converged[1]
                @warn "Exiting. Wake elliptic grid solve did not converge. Consider a looser convergence tolerance if the geometry looks good."
            end
        end
        #TODO: write a function that returns the same as outs below, but all zeros
        #TODO: probably just call  the post-process function directly and return a reset_container! of the output
        return [],#zero_outputs(),
        (; solve_parameter_tuple..., ivb, airfoils, idmaps, panels, problem_dimensions),
        false
    end

    return analyze(
        multipoint,
        propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=solve_container_caching,
        return_inputs=return_inputs,
    )
end

"""
    analyze(
        multipoint::Vector{OperatingPoint},
        propulsor::Propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options::Options=set_options();
        return_inputs=false,
        solve_container_caching=nothing,
    )

Analyze propulsor, assuming `setup_analysis` has been called and the inputs are being passed in here.

# Arguments
- `multipoint::AbstractVector{OperatingPoint}` : Vector of Operating Points at which to analyze the propulsor (note that the operating point within the propulsor input will be overwritten with these)
- `propulsor::Propulsor` : Propulsor input object
- `prepost_containers::NamedTuple` : An output from `setup_analysis` containing reshaped views into the prepost cache
- `solve_parameter_cache_vector::Vector` : An output from `setup_analysis` containing the relevant typed cache vector of solve parameters
- `solve_parameter_cache_dims::NamedTuple` : An output from `setup_analysis` containing dimensions used for reshaping the solve parameter cache
- `airfoils::Vector{AFType}` : An output from `setup_analysis` contiaining the blade element airfoil polar objects
- `A_bb_LU::LinearAlgebra.LU` : An output from `setup_analysis` that is the LU decomposition of the AIC matrix used in the panel method
- `idmaps::NamedTuple` : An output from `setup_analysis` containing bookkeeping information (index mappings)
- `problem_dimensions::NamedTuple` : An output from `setup_analysis` contiaining bookkeeping information (problem dimensions)
- `options::Options=set_options()` : Options object

# Keyword Arguments
- `solve_container_caching=nothing` : Output of `allocate_solve_container_cache`
- `return_inputs=false` : flag as to whether or not to return the pre-processed inputs

# Returns
- `outs::Vector{NamedTuple}` : Named Tuple of various analysis outputs (see docstring for postprocess for details), note, if linear system decomposition fails, no solve is performed and an empty vector is returned.
- `ins::NamedTuple` : Named Tuple of various pre-processed inputs (e.g. panels and body linear system), will only be returned if `return_inputs=true`.  Note that some inputs will be overwritten (e.g. the linear system RHS components related to the freestream) and only those associated with the final operating point will be returned.
- `convergence_flag` : Flag for successful solve convergence
"""
function analyze(
    multipoint::AbstractVector{TO},
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    idmaps,
    problem_dimensions,
    options::Options;
    solve_container_caching=nothing,
    return_inputs=false,
) where {TO<:OperatingPoint}
    if options.verbose
        println("\nRunning Multipoint Analysis")
    end

    # Set up Solve Container Cache
    if isnothing(solve_container_caching)
        solve_container_caching = allocate_solve_container_cache(
            options.solver_options, propulsor.paneling_constants
        )
    end

    # reset multipoint index
    options.multipoint_index[] = 0

    outs = [
        analyze_multipoint(
            op,
            propulsor,
            prepost_containers,
            solve_parameter_cache_vector,
            solve_parameter_cache_dims,
            airfoils,
            A_bb_LU,
            idmaps,
            problem_dimensions,
            options;
            solve_container_caching=solve_container_caching,
            return_inputs=false,
        ) for op in multipoint
    ]

    if return_inputs
        solve_parameter_tuple = withdraw_solve_parameter_cache(
            solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
        )

        return outs,
        (;
            prepost_containers.panels,
            prepost_containers.ivb,
            solve_parameter_tuple...,
            blade_elements=(; blade_elements..., airfoils...),
            linsys=(; linsys..., A_bb_LU),
            idmaps,
            problem_dimensions,
        ),
        options.solver_options.converged
    else
        return outs, options.solver_options.converged
    end
end

"""
    analyze_multipoint(
        operating_point::OperatingPoint,
        propulsor::Propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options::Options;
        solve_container_caching=nothing,
        return_inputs=false,
    )

Identical to the single analyze function assuming `setup_analysis` has been called; except here we are running a single operating point for a multipoint analysis, and overwriting the operating point in the propulsor with the explicit operating point input.
"""
function analyze_multipoint(
    operating_point::TO,
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    A_bb_LU,
    idmaps,
    problem_dimensions,
    options::Options;
    solve_container_caching=nothing,
    return_inputs=false,
) where {TO<:OperatingPoint}

    # incremenet operating point to keep track of convergence in multipoint solves
    options.multipoint_index[] += 1

    if options.verbose
        println("\n  Operating Point:")
        for fn in fieldnames(typeof(operating_point))
            println(@sprintf "    %6s = %5.3e" fn getfield(operating_point, fn)[])
        end
    end

    # Reshape solver parameter cache to update multipoint
    solve_parameter_tuple = withdraw_solve_parameter_cache(
        options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
    )

    # - copy over operating point - #
    for f in fieldnames(typeof(operating_point))
        update_operating_point!(propulsor.operating_point, operating_point)
        solve_parameter_tuple.operating_point[f] .= getfield(operating_point, f)
    end

    # - Set up Body Linear System RHS - #
    vinfvec = [operating_point.Vinf[]; 0.0]
    prepost_containers.vdnb .= [
        dot(vinfvec, nhat) for
        nhat in eachcol(prepost_containers.panels.body_vortex_panels.normal)
    ]
    prepost_containers.vdnpcp .= [
        dot(vinfvec, nhat) for
        nhat in eachcol(prepost_containers.panels.body_vortex_panels.itnormal)
    ]
    assemble_rhs_matrix!(
        solve_parameter_tuple.linsys.b_bf,
        prepost_containers.vdnb,
        prepost_containers.vdnpcp,
        prepost_containers.panels.body_vortex_panels.npanel,
        prepost_containers.panels.body_vortex_panels.nnode,
        prepost_containers.panels.body_vortex_panels.totpanel,
        prepost_containers.panels.body_vortex_panels.totnode,
        prepost_containers.panels.body_vortex_panels.prescribednodeidxs,
    )

    return analyze(
        propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        return_inputs=return_inputs,
        solve_container_caching=solve_container_caching,
    )[1]
end