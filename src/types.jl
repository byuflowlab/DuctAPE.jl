abstract type IntegrationMethod end

@kwdef struct Romberg{TI,TF,TP,TV} <: IntegrationMethod
    nsamples::TI = 1025
    stepsize::TF = 1.0 / 1024
    factors::TP = Primes.factor(1024)
    nfactors::TI = 10
    atol::TF = 1e-6
    samplepoints::TV = collect(range(0, 1, 1025))
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
    maxevals::TI = 1e2
    atol::TF = 1e-6
end
