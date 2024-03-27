abstract type IntegrationMethod end

@kwdef struct Romberg <: IntegrationMethod
    nsamples::Int = 1025
    stepsize::Float64 = 1.0 / 1024
    factors::Primes.Factorization{Int} = Primes.factor(1024)
    nfactors::Int = 10
    atol::Float64 = 1e-6
end

function set_romberg_options(; nsamples=1025, atol=1e-6)
    stepsize = 1.0 / (nsamples - 1)
    factors = Primes.factor(nsamples - 1)
    nfactors = sum(values(factors))

    return Romberg(;
        atol=atol, nsamples=nsamples, factors=factors, nfactors=nfactors, stepsize=stepsize
    )
end

@kwdef struct GaussKronrod <: IntegrationMethod
    order::Int = 3
    maxevals::Int = 1e2
    atol::Float64 = 1e-6
end
