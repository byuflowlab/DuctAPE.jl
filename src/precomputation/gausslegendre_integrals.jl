#---------------------------------#
#             VORTEX              #
#---------------------------------#

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_vortex_panel_integration!(
    integration_options::GaussLegendre,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)

    # - Sample Function - #
    for (s, t) in zip(eachrow(containers.samples), integration_options.sample_points)
        # Define function to integrate
        dt.nominal_vortex_induced_velocity_sample!(
            s, t, node1, node2, influence_length, controlpoint, containers.sample_cache
        )
    end

    # - Do Quadrature - #
    for (i, s) in enumerate(eachcol(containers.samples))
        V[i] = dot(integration_options.weights, s)
    end

    if debug
        return V, err
    else
        return V
    end
end

#---------------------------------#
#             SOURCE              #
#---------------------------------#

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_source_panel_integration!(
    integration_options::GaussLegendre,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)

    # - Sample Function - #
    for (s, t) in zip(eachrow(containers.samples), integration_options.sample_points)
        # Define function to integrate
        # TODO: write this function
        dt.nominal_source_induced_velocity_sample!(
            s, t, node1, node2, influence_length, controlpoint, containers.sample_cache
        )
    end

    # - Do Quadrature - #
    for (i, s) in enumerate(eachcol(containers.samples))
        V[i] = dot(integration_options.weights, s)
    end

    if debug
        return V, err
    else
        return V
    end
end

