abstract type IntegrationMethod end

@kwdef struct Romberg{TI,TF,TP,TV} <: IntegrationMethod
    max_subdivisions::TI = 10
    atol::TF = 1e-6
end

function set_romberg_options(; max_subdivisions=10, atol=1e-6)
    return Romberg(; max_subdivisions, atol)
end

@kwdef struct GaussKronrod{TI,TF} <: IntegrationMethod
    order::TI = 5
    maxevals::TI = 1000
    atol::TF = 1e-12
end

struct GaussLegendre{TN,TW} <: IntegrationMethod
    sample_points::TN
    weights::TW
end

function GaussLegendre(nsamples=20; silence_warnings=true)
    if silence_warnings && Bool((nsamples) % 2)
        @warn "Must have an even number of GaussLegendre sample points if using for panel self influence"
    end

    nodes, weights = FastGaussQuadrature.gausslegendre(nsamples)

    return GaussLegendre(linear_transform((-1, 1), (0, 1), nodes), weights ./ 2.0)
end

@kwdef struct IntegrationOptions{TN,TS}
    nominal::TN = GaussLegendre(20)
    singular::TS = GaussLegendre(20)
end
