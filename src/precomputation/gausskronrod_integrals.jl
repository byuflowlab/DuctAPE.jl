#---------------------------------#
#             VORTEX              #
#---------------------------------#

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_vortex_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    cache_vec;
    nondimrange=[0.0; 1.0],
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return dt.nominal_vortex_induced_velocity_sample(
            t,
            node1,
            node2,
            influence_length,
            controlpoint,
            cache_vec
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
    end
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function self_vortex_panel_integration(
    node1,
    node2,
    influence_length,
    controlpoint,
    cache_vec;
    nondimrange=[0.0; 1.0],
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_vortex_induced_velocity_sample(
            t,
            node1,
            node2,
            influence_length,
            controlpoint,
            cache_vec;
            nondimrange=nondimrange,
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        0.5,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    cache_vec[1], cache_vec[2] = analytically_integrated_vortex_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[1:2] .+= cache_vec[1] / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end

#---------------------------------#
#             SOURCE              #
#---------------------------------#

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_source_panel_integration(
    node1,
    node2,
    influence_length,
    controlpoint,
    cache_vec;
    nondimrange=[0.0; 1.0],
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return nominal_source_induced_velocity_sample(
            t,
            node1,
            node2,
            influence_length,
            controlpoint,
            cache_vec;
            nondimrange=nondimrange,
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
    end
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function self_source_panel_integration(
    node1,
    node2,
    influence_length,
    controlpoint,
    cache_vec;
    nondimrange=[0.0; 1.0],
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_source_induced_velocity_sample(
            t,
            node1,
            node2,
            influence_length,
            controlpoint,
            cache_vec;
            nondimrange=nondimrange,
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        0.5,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    cache_vec[1], cache_vec[2] = analytically_integrated_source_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[3:4] .+= cache_vec[2] / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end
