"""
    allocate_integration_containers(
        integration_options::IntegrationMethod, dispatch_type; cache_size=20
    )

Description

# Arguments
- `integration_options::IntegrationMethod` : options for integration used for dispatch and to size cache
- `dispatch_type::` : an object with `eltype(dispatch_type)` with which to define the type for cache initialization.

# Keyword Arguments
- `cache_size::Int=20` : size needed for intermediate calculations for integration.

# Returns
- `integration_containers::NamedTuple` : A named tuple containing the cache(s) needed for integration.
"""
function allocate_integration_containers(
    integration_options::Romberg, dispatch_type; cache_size=20
)
    TF = eltype(dispatch_type)
    return (;
        sample_cache=zeros(TF, cache_size),
        samples=[(zeros(TF, 4), TF[0.0]) for i in 1:(integration_options.max_subdivisions)],
    )
end

function allocate_integration_containers(
    integration_options::GaussLegendre, dispatch_type; cache_size=20
)
    TF = eltype(dispatch_type)
    return (;
        sample_cache=zeros(TF, cache_size),
        samples=zeros(TF, length(integration_options.sample_points), 4),
    )
end

function allocate_integration_containers(
    integration_options::GaussKronrod, dispatch_type; cache_size=20
)
    TF = eltype(dispatch_type)
    return zeros(TF, cache_size)
end
