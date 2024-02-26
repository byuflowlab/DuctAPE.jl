"""
"""
function post_process_iad(Gamr, gamw, sigr, fan_params, cache_and_dims)

    # - Extract fan_params - #
    (;
        duct_coordinates, # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, # Vector of NamedTuples of a bunch of stuff...
        paneling_constants, # NamedTuple of numbers and vectors of numbers
        freestream, # NamedTuple of numbers
        reference_parameters, # NamedTuple of numbers
    ) = fan_params

    # - Extract Cache and Dims - #
    (; cache, dimensions) = cache_and_dims

    # - re-populate cache in this scope so that the derivatives between the fan_params and the post-processing are visible. - #
    # TODO: write this function
    various_postprocessing_parameters = postcompute_parameters_aid!(cache, fan_params)

    # - call solver residual to get the most up-to-date stuff in the cache - #

    # - Do the rest of the post processing - #

    return outs
end

