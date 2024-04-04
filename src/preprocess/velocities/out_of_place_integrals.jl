#---------------------------------#
#             VORTEX              #
#---------------------------------#

"""
"""
function nominal_vortex_panel_integration(
    integration_options,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)
    TF = promote_type(eltype(node1), eltype(node2))
    V = zeros(TF, 2, 2)
    return nominal_vortex_panel_integration!(
        integration_options,
        V,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=debug,
    )
end

function self_vortex_panel_integration(
    integration_options,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)
    TF = promote_type(eltype(node1), eltype(node2))
    V = zeros(TF, 2, 2)
    return self_vortex_panel_integration!(
        integration_options,
        V,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=debug,
    )
end

#---------------------------------#
#             SOURCE              #
#---------------------------------#

function nominal_source_panel_integration(
    integration_options,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)
    TF = promote_type(eltype(node1), eltype(node2))
    V = zeros(TF, 2, 2)
    return nominal_source_panel_integration!(
        integration_options,
        V,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=debug,
    )
end

function self_source_panel_integration(
    integration_options,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)
    TF = promote_type(eltype(node1), eltype(node2))
    V = zeros(TF, 2, 2)

    return self_source_panel_integration!(
        integration_options,
        V,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=debug,
    )
end
