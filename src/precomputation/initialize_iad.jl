"""
"""
function repanel(propuslor)
    return repanel(
        propulsor.duct_coordinates,
        propulsor.centerbody_coordinates,
        propulsor.rotorstator_parameters,
        propulsor.paneling_constants,
    )
end

"""
    repanel(duct_coordinates, centerbody_coordinates, rotorstator_parameters, paneling_constants)
    repanel(propulsor)

Repanel geometry in preparation for wake definition.
"""
function repanel(
    duct_coordinates, centerbody_coordinates, rotorstator_parameters, paneling_constants
)

    #TODO; what other panels are actually needed? do you need the body wake panels or no?
    return (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
end

"""
"""
function calculate_unit_induced_velocities(cache, cache_dims, panels)
    return (; vz_, vr_)
end

"""
"""
function initialize_linear_system(cache, cache_dims, induced_velocities, body_vortex_panels)
    return (; A_bb, A_bb_LU, lu_decomp_flag=[issuccess(A_bb_LU)], b_bf)
end

"""
"""
function precompute_parameters_iad!(cache, propulsor)

    # - Extract propulsor - #
    (;
        duct_coordinates, # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, # Vector of NamedTuples of a bunch of stuff...
        paneling_constants, # NamedTuple of numbers and vectors of numbers
        freestream, # NamedTuple of numbers
        reference_parameters, # NamedTuple of numbers
    ) = propulsor

    # - Rename for Convenience - #
    rsp = rotorstator_parameters

    # - Panel Everything - #

    # - Compute Influence Matrices - #
    # Note: cache contains containers for these

    # - Assemble Body Parameters - #
    body_parameters = (;)

    # - Assemble Rotor Parameters - #
    rotor_parameters = (;
        Omega=rsp.Omega, B=rsp.B, rotor_panel_center=rotor_source_panels.controlpoint[2, :]
    )

    # - Assemble Wake Parameters - #
    wake_parameters = (;)

    return body_parameters,
    rotor_parameters, wake_parameters, freestream,
    reference_parameters
end

"""
"""
function postcompute_parameters_aid!(cache, propulsor)

    # - Re-compute the precomputed stuff - #
    bp, rp, wp, fp, refp = precompute_parameters_iad!(cache, propulsor)

    # - Compute the rest of the stuff required for post-processing - #

    return nothing
end
