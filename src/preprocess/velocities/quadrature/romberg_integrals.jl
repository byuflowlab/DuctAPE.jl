"""
    extrapolate!(V, err, fh; power=2, atol=1e-6)

Performs Richardson extrapolation on an array `fh` for use in Romberg integration.

# Arguments
- `V::Matrix{Float}` : velocity components due to the jth and j+1th nodes in the format: `[vz_j vr_j; vz_{j+1} vr_{j+1}]`
- `err::Vector{Float}` : estimated errors in velocity approximation
- `fh::Tuple` : `(f(h), h)` tuples (in order of decreasing `|h|`)

"""
function extrapolate!(V, err, fh; power=2, atol=1e-6)

    # - rename for convenience - #
    (f0, h) = first(fh)

    # - initialize - #
    err .= Inf
    numeval = 1
    maxeval = size(fh, 1)

    # - Loop - #
    while numeval < maxeval

        # - increment - #
        numeval += 1

        # - update - #
        (f_prime, h_prime) = fh[numeval + (firstindex(fh) - 1)]
        h = h_prime
        minerr_prime = similar(err) .= Inf
        f_ip1 = f_prime

        # - Inner Loop - #
        for i in (numeval - 1):-1:1
            f_i, h_i = fh[i + (firstindex(fh) - 1)]
            c = (h_i[] / h_prime[])^power
            @. f_ip1 += (f_ip1 - f_i) / (c - 1)
            fh[i + (firstindex(fh) - 1)] = (f_ip1, h_i)
            err_prime = norm.(f_ip1 - f_i)
            minerr_prime = min.(minerr_prime, err_prime)
            if err_prime < err
                V[:], err = f_ip1, err_prime
            end
        end

        # - converged - #
        # all(err .<= atol) && break
    end

    return V, err
end

#---------------------------------#
#             VORTEX              #
#---------------------------------#

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
    reset_containers!(containers)

    # - Loop through number of subdivisions (start with 2) - #
    for ii in 1:(integration_options.max_subdivisions)
        nint = 2^ii

        # - Set step size for integration and extrapolation - #
        dx = 1.0 / nint
        containers.samples[ii][2] .= dx

        # - Loop through intervals for current subdivision level - #
        for i in 1:nint
            # sample at the interval midpoints
            t = (i - 0.5) / nint

            # get sample velocity
            nominal_vortex_induced_velocity_sample!(
                @view(containers.samples[ii])[1][1],
                t,
                node1,
                node2,
                influence_length,
                controlpoint,
                containers.sample_cache,
            )
        end

        # multiply samples by the interval length
        containers.samples[ii][1] .*= dx

        if ii > 1
            # extrapolate once you can
            extrapolate!(
                V,
                @view(containers.sample_cache[1:4]),
                @view(containers.samples[1:ii]);
                atol=integration_options.atol,
            )

            # if you have met the convergence requirements, terminate
            all(@view(containers.sample_cache[1:4]) .<= integration_options.atol) && break
        end
    end

    return V
end

function self_vortex_panel_integration!(
    integration_options::Romberg,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)
    reset_containers!(containers)
    # - Loop through number of subdivisions (start with 2) - #
    for ii in 1:(integration_options.max_subdivisions)
        nint = 2^ii

        # - Set step size for integration and extrapolation - #
        dx = 1.0 / nint
        containers.samples[ii][2] .= dx

        # - Loop through intervals for current subdivision level - #
        for i in 1:nint
            # sample at the interval midpoints
            t = (i - 0.5) / nint

            # get sample velocity
            self_vortex_induced_velocity_sample!(
                @view(containers.samples[ii])[1][1],
                t,
                node1,
                node2,
                influence_length,
                controlpoint,
                containers.sample_cache,
            )
        end

        # multiply samples by the interval length
        containers.samples[ii][1] .*= dx
        # multiply samples by the influence length (not done in self induced unit velocity function)
        containers.samples[ii][1] .*= influence_length

        # - Add in the analytical bits - #
        containers.sample_cache[1], containers.sample_cache[2] = analytically_integrated_vortex_influence(
            controlpoint[2], influence_length
        )
        containers.samples[ii][1][1:2] .+= containers.sample_cache[1] / 2.0

        if ii > 1
            # extrapolate once you can
            extrapolate!(
                V,
                @view(containers.sample_cache[1:4]),
                @view(containers.samples[1:ii]);
                atol=integration_options.atol,
            )

            # if you have met the convergence requirements, terminate
            all(@view(containers.sample_cache[1:4]) .<= integration_options.atol) && break
        end
    end

    if debug
        return V, err
    else
        return V
    end
end

#---------------------------------#
#             SOURCE              #
#---------------------------------#

function nominal_source_panel_integration!(
    integration_options::Romberg,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)
    reset_containers!(containers)

    # - Loop through number of subdivisions (start with 2) - #
    for ii in 1:(integration_options.max_subdivisions)
        nint = 2^ii

        # - Set step size for integration and extrapolation - #
        dx = 1.0 / nint
        containers.samples[ii][2] .= dx

        # - Loop through intervals for current subdivision level - #
        for i in 1:nint
            # sample at the interval midpoints
            t = (i - 0.5) / nint

            # get sample velocity
            nominal_source_induced_velocity_sample!(
                @view(containers.samples[ii])[1][1],
                t,
                node1,
                node2,
                influence_length,
                controlpoint,
                containers.sample_cache,
            )
        end

        # multiply samples by the interval length
        containers.samples[ii][1] .*= dx

        if ii > 1
            # extrapolate once you can
            extrapolate!(
                V,
                @view(containers.sample_cache[1:4]),
                @view(containers.samples[1:ii]);
                atol=integration_options.atol,
            )

            # if you have met the convergence requirements, terminate
            all(@view(containers.sample_cache[1:4]) .<= integration_options.atol) && break
        end
    end

    return V
end

function self_source_panel_integration!(
    integration_options::Romberg,
    V,
    node1,
    node2,
    influence_length,
    controlpoint,
    containers;
    debug=false,
)
    reset_containers!(containers)
    # - Loop through number of subdivisions (start with 2) - #
    for ii in 1:(integration_options.max_subdivisions)
        nint = 2^ii

        # - Set step size for integration and extrapolation - #
        dx = 1.0 / nint
        containers.samples[ii][2] .= dx

        # - Loop through intervals for current subdivision level - #
        for i in 1:nint
            # sample at the interval midpoints
            t = (i - 0.5) / nint

            # get sample velocity
            self_source_induced_velocity_sample!(
                @view(containers.samples[ii])[1][1],
                t,
                node1,
                node2,
                influence_length,
                controlpoint,
                containers.sample_cache,
            )
        end

        # multiply samples by the interval length
        containers.samples[ii][1] .*= dx
        # multiply samples by the influence length (not done in self induced unit velocity function)
        containers.samples[ii][1] .*= influence_length

        # - Add in the analytical bits - #
        containers.sample_cache[1], containers.sample_cache[2] = analytically_integrated_source_influence(
            controlpoint[2], influence_length
        )
        containers.samples[ii][1][3:4] .+= containers.sample_cache[2] / 2.0

        if ii > 1
            # extrapolate once you can
            extrapolate!(
                V,
                @view(containers.sample_cache[1:4]),
                @view(containers.samples[1:ii]);
                atol=integration_options.atol,
            )

            # if you have met the convergence requirements, terminate
            all(@view(containers.sample_cache[1:4]) .<= integration_options.atol) && break
        end
    end

    if debug
        return V, err
    else
        return V
    end
end

