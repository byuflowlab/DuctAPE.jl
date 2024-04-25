#---------------------------------#
#             VORTEX              #
#---------------------------------#

function nominal_vortex_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return nominal_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache
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

function self_vortex_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
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

    sample_cache[1], sample_cache[2] = analytically_integrated_vortex_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[1:2] .+= sample_cache[1] / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end

#---------------------------------#
#             SOURCE              #
#---------------------------------#

function nominal_source_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return nominal_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
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

function self_source_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
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

    sample_cache[1], sample_cache[2] = analytically_integrated_source_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[3:4] .+= sample_cache[2] / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end
