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