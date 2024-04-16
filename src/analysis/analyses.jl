"""
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
    problem_dimensions, prepost_containers, solve_parameter_cache_vector, solve_parameter_cache_dims, ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps = setup_analysis(
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

    # - Continue with Analysis - #
    return analyze(
        propulsor,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        airfoils,
        ivb,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=solve_container_caching,
        return_inputs=return_inputs,
    )
end

"""
"""
function analyze(
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    ivb,
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
            solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
        )

        return outs,
        (;
            prepost_containers.panels,
            prepost_containers.ivb,
            solve_parameter_tuple...,
            blade_elements=(; blade_elements..., airfoils...),
            linsys=(; linsys..., A_bb_LU),
        ),
        options.solver_options.converged[]
    else
        return outs, options.solver_options.converged[options.multipoint_index[]]
    end
end

"""
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
    problem_dimensions, prepost_containers, solve_parameter_cache_vector, solve_parameter_cache_dims, ivb, A_bb_LU, lu_decomp_flag, airfoils, idmaps = setup_analysis(
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
        ivb,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        solve_container_caching=solve_container_caching,
        return_inputs=return_inputs,
    )
end

function analyze(
    multipoint::AbstractVector{TO},
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    ivb,
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
            ivb,
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
        ),
        options.solver_options.converged
    else
        return outs, options.solver_options.converged
    end
end

function analyze_multipoint(
    operating_point::TO,
    propulsor::Propulsor,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
    airfoils,
    ivb,
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
    println(propulsor.operating_point)
    for f in fieldnames(typeof(operating_point))
        update_operating_point!(propulsor.operating_point, operating_point)
        solve_parameter_tuple.operating_point[f] .= getfield(operating_point, f)
    end
    println(propulsor.operating_point)

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
        ivb,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options;
        return_inputs=return_inputs,
        solve_container_caching=solve_container_caching,
    )[1]
end
