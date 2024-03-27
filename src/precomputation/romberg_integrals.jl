#TODO: decide if dfdc method or julia method is better
# the DFDC method requires the whole extrapolation, but you can limit the number of intervals to evaluate when the integral is easy.
# when the integral is hard though, you end up doing 2x more calculations than you would have if you set your sample size larger to begin with.
# the julia method requires you to evaluate all the sample points up front, but will terminate the extrapolation as soon as you hit your desired tolerance. So you need to have a large number of samples to start with.
# DFDC doesn't actually do a true trapezoidal integration, rather it assumes a linear function and does δx*f(midpoint of iterval), this is technically less accurate than the trapezoidal method, but does avoid sampling at the end and midpoints of panels.  Theoretically, the extrapolation will account for the less accurate integration method.
# perhaps a combo, where you do a DFDC-like integration (not quite trapezoidal), but do the large sample size and early extrapolation termination.


"""
    extrapolate!(fh::AbstractVector; power=1, atol=0, rtol=0, maxeval=typemax(Int), breaktol=Inf)

Similar to `extrapolate(fh)`, performs Richardson extrapolation on an array `fh`
of `(f(h), h)` tuples (in order of decreasing `|h|`), but overwrites the array
`fh` in-place with intermediate calculations.

(Thus, the array `fh` must be a vector of `Tuple{T,H}` values, where `H<:Number` is
the type of `h` and `T` is the type of the extrapolated `f(0)` **result**.  This `T`
should be a floating-point type, i.e. `fh` should contain `float(f(h))` if the
function you are extrapolating is not already floating-point-valued.)
"""
function extrapolate!(
    fh::AbstractVector{<:Tuple{TF1,TF2}};
    power::Int=2,
    atol::Float64=1e-6,
    rtol::Float64=0.0,
    breaktol::Float64=Inf,
    maxeval::Int=typemax(Int),
) where {TF1,TF2}

    # - Check for Problems - #
    (rtol >= 0 && atol >= zero(atol)) ||
        throw(ArgumentError("rtol and atol must be nonnegative"))
    breaktol > 0 || throw(ArgumentError("breaktol must be positive"))
    isempty(fh) && throw(ArgumentError("(f,h) array must be non-empty"))

    # - rename for convenience - #
    (f0, h) = first(fh)

    # - initialize - #
    err::typeof(float(norm(f0))) = Inf
    numeval = 1
    maxeval = min(maxeval, length(fh))

    # - Loop - #
    while numeval < maxeval
        # - increment - #
        numeval += 1

        # - update - #
        (f_prime, h_prime) = fh[numeval + (firstindex(fh) - 1)]
        abs(h) > abs(h_prime) ||
            throw(ArgumentError("|$h_prime| >= |$h| is not decreasing"))
        h = h_prime
        minerr_prime = oftype(err, Inf)
        f_ip1 = f_prime

        # - Inner Loop - #
        for i in (numeval - 1):-1:1
            f_i, h_i = fh[i + (firstindex(fh) - 1)]
            c = (h_i / h_prime)^power
            f_ip1 += (f_ip1 - f_i) / (c - 1)
            fh[i + (firstindex(fh) - 1)] = (f_ip1, h_i)
            err_prime = norm(f_ip1 - f_i)
            minerr_prime = min(minerr_prime, err_prime)
            if err_prime < err
                f0, err = f_ip1, err_prime
            end
        end

        # - stop early if error increases too much - #
        (minerr_prime > breaktol * err || !isfinite(minerr_prime)) && break

        # - converged - #
        err <= max(rtol * norm(f0), atol) && break
    end
    return (f0, err)
end

@views function romberg!(fh, step_size, samples, endsum, factors, numfactors; kwargs...)

    # - Do Trapezoidal Integration - #
    b, e = firstindex(samples), lastindex(samples)
    i = numfactors + 1
    fh[i] = (step_size * (sum(samples[(b + 1):(e - 1)]) + endsum), step_size)
    nstep = 1
    for (j, K) in factors
        @inbounds for k in 1:K
            nstep *= j
            sdx = nstep * step_size
            if i == 2 # last iteration (empty sum)
                fh[1] = (sdx * endsum, sdx)
            else
                fh[i -= 1] = (
                    sdx * (sum(samples[(b + nstep):nstep:(e - nstep)]) + endsum), sdx
                )
            end
        end
    end

    # - Do Richardson Extrapolation until converged - #
    return extrapolate!(fh; power=2, kwargs...)
end

#---------------------------------#
#             VORTEX              #
#---------------------------------#

function allocate_integration_containers(dispatch_type::TF, nsamples, nfactors) where {TF}
    return (;
        err=zeros(TF, 2, 2), #always this size
        sample_cache=zeros(TF, 20), #always this size, this is for the intermediate calcs
        samples=zeros(TF, nsamples, 4), # always 4 wide
        fh=[(TF(1.0), TF(1.0)) for i in 1:(nfactors + 1)],
    )
end

"""
`V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: [vz_j vr_j; vz_{j+1} vr_{j+1}]
"""
function nominal_vortex_panel_integration!(
    integration_options::Romberg,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)

    # Define function to integrate
    for (s, t) in zip(eachrow(containers.samples), integration_options.samplepoints)
        s .= nominal_vortex_induced_velocity_sample(
            t,
            node1,
            node2,
            influence_length,
            controlpoint,
            containers.sample_cache;
            nondimrange=integration_options.nondimrange,
        )
    end

    for (i, s) in enumerate(eachcol(containers.samples))
        V[i], containers.err[i] = romberg!(
            containers.fh,
            integration_options.stepsize,
            s,
            (s[1] + s[end]) / 2.0,
            factors,
            length(factors);
            atol=integration_options.atol,
        )
    end

    if debug
        # return reshape(V, (2, 2)), err
        return V, err
    else
        # return reshape(V, (2, 2))
        return V
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

    V, err = quadgk(fsample, 0.0, 0.5, 1.0; order=3, maxevals=1e2, atol=1e-6)

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

    V, err = quadgk(fsample, 0.0, 1.0; order=3, maxevals=1e2, atol=1e-6)

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

    V, err = quadgk(fsample, 0.0, 0.5, 1.0; order=3, maxevals=1e2, atol=1e-6)

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
