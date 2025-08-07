"""
    analyze_with_warm_start(
        ducted_rotor::DuctedRotor,
        operating_point::Vector{OperatingPoint},
        reference_parameters::ReferenceParameters,
        prepost_containers,
        solve_parameter_cache_vector,
        solve_parameter_cache_dims,
        A_bb_LU,
        idmaps,
        problem_dimensions,
        options::Options=set_options();
        return_inputs=false,
        solve_container_caching=nothing,
    )

Analyze `ducted_rotor`, assuming `setup_analysis` has been called and the inputs are being passed in here.

# Arguments
- `ducted_rotor::DuctedRotor` : DuctedRotor input object
- `operating_point::AbstractVector{OperatingPoint}` : Vector of Operating Points at which to analyze the ducted_rotor
- `reference_parameters::ReferenceParameters` : ReferenceParameters input object (see docstring for `ReferenceParameters` type)
- `prepost_containers::NamedTuple` : An output from `setup_analysis` containing reshaped views into the prepost cache
- `solve_parameter_cache_vector::Vector` : An output from `setup_analysis` containing the relevant typed cache vector of solve parameters
- `solve_parameter_cache_dims::NamedTuple` : An output from `setup_analysis` containing dimensions used for reshaping the solve parameter cache
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
function analyze_with_warm_start(
    ducted_rotor::DuctedRotor,
    operating_point::AbstractVector{TO},
    reference_parameters::ReferenceParameters,
    prepost_containers,
    solve_parameter_cache_vector,
    solve_parameter_cache_dims,
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
            options.solver_options, ducted_rotor.paneling_constants
        )
    end

    # reset multipoint index
    options.multipoint_index[] = 0

    # Run like normal
    outs = [
        analyze_multipoint(
            ducted_rotor,
            op,
            reference_parameters,
            prepost_containers,
            solve_parameter_cache_vector,
            solve_parameter_cache_dims,
            A_bb_LU,
            idmaps,
            problem_dimensions,
            options;
            solve_container_caching=solve_container_caching,
            return_inputs=false,
        ) for op in operating_point
    ]

    nop = length(operating_point)
    run_counts = ones(nop)
    # If all converged, proceed as normal
    while !all(options.solver_options.converged)
        println("WARM STARTING!!!")
        # If all unconverged, check residuals
        if !any(options.solver_options.converged)
            println("  All Failed")
            if all(options.solver_options.residuals .> 1.0)
                println("    None have residuals less than 1.0")
                # If all residuals greater than 1
                ##TODO: set up method to skip solve and return values based on initialization
                options.multipoint_index[] = 0
                options.solver_options.skip_solve .= true
                run_counts .+= 1

                outs = [
                    analyze_multipoint(
                        ducted_rotor,
                        op,
                        reference_parameters,
                        prepost_containers,
                        solve_parameter_cache_vector,
                        solve_parameter_cache_dims,
                        A_bb_LU,
                        idmaps,
                        problem_dimensions,
                        options;
                        solve_container_caching=solve_container_caching,
                        return_inputs=false,
                    ) for op in operating_point
                ]

                break

            elseif all(options.solver_options.residuals .<= 1.0)
                println("    All have residuals less than 1.0")
                # If all residuals less than 1, return as normal
                break
            else
                println("    Some have residuals less than 1.0")
                # If any residuals less than 1, try running cases with residuals greater than 1 and the output of the less than 1 residuals to see if you can do better.

                # determine cases to re-run and closest reasonable case
                re_runs = find_large_and_nearest_small(options.solver_options.residuals)
                println("    re_run check: ", re_runs)

                while !isempty(re_runs)
                    println("      Attempting Re-runs")
                    # loop through re-runs
                    for rr in re_runs
                        if run_counts[rr[1]] < 2

                            # run closest case
                            run_counts[rr[2]] += 1
                            options.multipoint_index[] = rr[2] - 1

                            _ = analyze_multipoint(
                                ducted_rotor,
                                operating_point[rr[2]],
                                reference_parameters,
                                prepost_containers,
                                solve_parameter_cache_vector,
                                solve_parameter_cache_dims,
                                A_bb_LU,
                                idmaps,
                                problem_dimensions,
                                options;
                                solve_container_caching=solve_container_caching,
                                return_inputs=false,
                            )

                            # run unconverged case
                            run_counts[rr[1]] += 1
                            options.multipoint_index[] = rr[1] - 1
                            options.solver_options.warm_start[rr[1]] = true

                            outs[rr[1]] = analyze_multipoint(
                                ducted_rotor,
                                operating_point[rr[1]],
                                reference_parameters,
                                prepost_containers,
                                solve_parameter_cache_vector,
                                solve_parameter_cache_dims,
                                A_bb_LU,
                                idmaps,
                                problem_dimensions,
                                options;
                                solve_container_caching=solve_container_caching,
                                return_inputs=false,
                            )
                        end
                    end

                    # re-check (also updates if closer case now exists)
                    println("      re_run re-check: ", re_runs)
                    re_runs = find_large_and_nearest_small(options.solver_options.residuals)

                    if all(run_counts[getindex.(re_runs, 2)] .>= 2)
                        break
                    end
                end
                # After re-running, if any residuals greater than 1 we should hit the else in the largest scope below and it'll skip solve those cases.

            end

        elseif any(options.solver_options.converged) && all(run_counts .< 2)
            println("  Some Converged")
            # If some converged, re-run unconverged cases starting at nearest converged case

            # determine cases to re-run and closest converged case
            re_runs = find_false_and_nearest_true_sorted(options.solver_options.converged)
            println("  re_run check: ", re_runs)
            re_run_count = [length(re_runs)]

            while re_run_count[1] > 0
                println("    Attempting Re-runs")
                # loop through re-runs
                for rr in re_runs
                    if run_counts[rr[1]] < 2

                        # run closest case
                        run_counts[rr[2]] += 1
                        options.multipoint_index[] = rr[2] - 1

                        _ = analyze_multipoint(
                            ducted_rotor,
                            operating_point[rr[2]],
                            reference_parameters,
                            prepost_containers,
                            solve_parameter_cache_vector,
                            solve_parameter_cache_dims,
                            A_bb_LU,
                            idmaps,
                            problem_dimensions,
                            options;
                            solve_container_caching=solve_container_caching,
                            return_inputs=false,
                        )

                        # run unconverged case
                        run_counts[rr[1]] += 1
                        options.multipoint_index[] = rr[1] - 1
                        options.solver_options.warm_start[rr[1]] = true

                        outs[rr[1]] = analyze_multipoint(
                            ducted_rotor,
                            operating_point[rr[1]],
                            reference_parameters,
                            prepost_containers,
                            solve_parameter_cache_vector,
                            solve_parameter_cache_dims,
                            A_bb_LU,
                            idmaps,
                            problem_dimensions,
                            options;
                            solve_container_caching=solve_container_caching,
                            return_inputs=false,
                        )
                    end
                end

                # re-check (also updates if closer case now exists)
                re_check = find_false_and_nearest_true_sorted(
                    options.solver_options.converged
                )
                println("  re_run re-check: ", re_check)
                re_run_count[1] = min(length(re_check), re_run_count[1])

                if all(run_counts[getindex.(re_runs, 1)] .>= 2)
                    break
                end
            end
            # After re-running, if any residuals greater than 1 we should hit the else in the largest scope below and it'll skip solve those cases.

        else #if any converged and we've tried everything again.
            println("  Tried. Giving up.")
            # We've already tried things once, so find the still unconverged cases and return what is needed
            for n in 1:nop
                if !options.solver_options.converged[n]
                    if options.solver_options.residuals[n] > 1.0
                        # run skip solve if there's no hope.
                        options.multipoint_index[] = n - 1
                        options.solver_options.skip_solve[n] = true
                        run_counts[n] += 1
                        outs[n] = analyze_multipoint(
                            ducted_rotor,
                            operating_point[n],
                            reference_parameters,
                            prepost_containers,
                            solve_parameter_cache_vector,
                            solve_parameter_cache_dims,
                            A_bb_LU,
                            idmaps,
                            problem_dimensions,
                            options;
                            solve_container_caching=solve_container_caching,
                            return_inputs=false,
                        )
                    end
                end
            end
        end
    end

    if return_inputs
        solve_parameter_tuple = withdraw_solve_parameter_cache(
            options.solver_options, solve_parameter_cache_vector, solve_parameter_cache_dims
        )

        return outs,
        (;
            prepost_containers.panels,
            prepost_containers.ivb,
            solve_parameter_tuple...,
            blade_elements=solve_parameter_tuple.blade_elements,
            linsys=(; solve_parameter_tuple.linsys..., A_bb_LU),
            idmaps,
            problem_dimensions,
            reference_parameters,
            multi_point=operating_point,
        ),
        options.solver_options.converged
    else
        return outs, options.solver_options.converged
    end
end
