#---------------------------------#
#             VORTEX              #
#---------------------------------#

"""
    nominal_vortex_panel_integration(
        integration_options,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=false,
    )

Integration of vortex panel induced velocity on a point far away.

# Arguments
- `integration_options::IntegrationMethod` : options for itegration methods
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `containers::NamedTuple` : cache for intermediate calculations

# Keyword Arguments
- `debug::Bool=false` : if true, some methods will return the estimation error.

# Returns
- `V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: `[vz_j vr_j; vz_{j+1} vr_{j+1}]`
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

"""
    self_vortex_panel_integration(
        integration_options,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=false,
    )

Integration of linear vortex panel self-induced velocity.

# Arguments
- `integration_options::IntegrationMethod` : options for itegration methods
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `containers::NamedTuple` : cache for intermediate calculations

# Keyword Arguments
- `debug::Bool=false` : if true, some methods will return the estimation error.

# Returns
- `V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: `[vz_j vr_j; vz_{j+1} vr_{j+1}]`
"""
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

"""
    nominal_source_panel_integration(
        integration_options,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=false,
    )

Integration of source panel induced velocity on a point far away.

# Arguments
- `integration_options::IntegrationMethod` : options for itegration methods
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `containers::NamedTuple` : cache for intermediate calculations

# Keyword Arguments
- `debug::Bool=false` : if true, some methods will return the estimation error.

# Returns
- `V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: `[vz_j vr_j; vz_{j+1} vr_{j+1}]`
"""
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

"""
    self_source_panel_integration(
        integration_options,
        node1,
        node2,
        influence_length,
        controlpoint,
        containers;
        debug=false,
    )

Integration of linear source panel self-induced velocity.

# Arguments
- `integration_options::IntegrationMethod` : options for itegration methods
- `node1::Vector{Float}` : first panel node (edge) position.
- `node2::Vector{Float}` : second panel node (edge) position.
- `influence_length::Float` : dimensional length of panel.
- `controlpoint::Vector{Float}` : controlpoint position
- `containers::NamedTuple` : cache for intermediate calculations

# Keyword Arguments
- `debug::Bool=false` : if true, some methods will return the estimation error.

# Returns
- `V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: `[vz_j vr_j; vz_{j+1} vr_{j+1}]`
"""
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
