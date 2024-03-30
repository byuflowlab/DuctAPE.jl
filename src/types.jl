abstract type IntegrationMethod end

@kwdef struct Romberg{TI,TF,TP,TV} <: IntegrationMethod
    nsamples::TI = 1025
    stepsize::TF = 1.0 / 1024
    factors::TP = Primes.factor(1024)
    nfactors::TI = 10
    atol::TF = 1e-6
    sample_points::TV = collect(range(0, 1, 1025))
end

function set_romberg_options(; nsamples=1025, atol=1e-6)
    stepsize = 1.0 / (nsamples - 1)
    factors = Primes.factor(nsamples - 1)
    nfactors = sum(values(factors))

    return Romberg(;
        atol=atol, nsamples=nsamples, factors=factors, nfactors=nfactors, stepsize=stepsize
    )
end

@kwdef struct GaussKronrod{TI,TF} <: IntegrationMethod
    order::TI = 3
    maxevals::TI = 100
    atol::TF = 1e-6
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
    singular::TS = GaussKronrod(3, 1e2, 1e-6)
end
