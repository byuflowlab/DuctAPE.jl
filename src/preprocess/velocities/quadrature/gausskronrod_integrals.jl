#---------------------------------#
#             VORTEX              #
#---------------------------------#
"""
    nominal_vortex_panel_integration(
        integration_options::GaussKronrod,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=false,
    ) -> Matrix{Float64} or (Matrix{Float64}, Float64)

Computes the velocity influence matrix at a control point due to a nominal vortex panel, using Gauss-Kronrod quadrature over the unit interval.

# Arguments
- `integration_options::GaussKronrod`: Settings for numerical integration, including order, tolerance, and maximum evaluations.
- `node1`: Start point of the vortex panel (2-element vector).
- `node2`: End point of the vortex panel (2-element vector).
- `influence_length`: Scaling factor representing the effective influence length of the vortex.
- `controlpoint`: The point at which the induced velocity is evaluated (2-element vector).
- `sample_cache`: Workspace array used internally to store intermediate values during integration.

# Keyword Arguments
- `debug::Bool=false`: If `true`, also returns the estimated integration error.

# Returns
- `2×2 Matrix{Float64}`: Induced velocity matrix at the control point.
- If `debug=true`, also returns the integration error `::Float64` as a second return value.
"""
function nominal_vortex_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return nominal_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
    end
end

"""
    self_vortex_panel_integration(
        integration_options::GaussKronrod,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=false,
    ) -> Matrix{Float64} or (Matrix{Float64}, Float64)

Computes the self-induced velocity influence matrix at a control point located on the same vortex panel, using a combination of numerical and analytical integration techniques.

# Description
This function accounts for the singularity in self-influence by applying Gauss-Kronrod quadrature over a split interval `[0.0, 0.5, 1.0]` and adding an analytical correction to improve accuracy near the control point. The result is scaled by the influence length and includes a correction from the analytically integrated near-field effect.

# Arguments
- `integration_options::GaussKronrod`: Configuration for numerical integration (order, max evaluations, tolerance).
- `node1`: Start point of the vortex panel (2-element vector).
- `node2`: End point of the vortex panel (2-element vector).
- `influence_length`: Length scale of the panel, used to scale the final result.
- `controlpoint`: The evaluation point located on the same panel (2-element vector).
- `sample_cache`: Workspace array reused for intermediate quantities, including analytic results.

# Keyword Arguments
- `debug::Bool=false`: If `true`, also returns the estimated integration error.

# Returns
- `2×2 Matrix{Float64}`: The reshaped velocity influence matrix at the control point.
- If `debug=true`, also returns the quadrature error estimate `::Float64`.
"""
function self_vortex_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_vortex_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        0.5,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    sample_cache[1], sample_cache[2] = analytically_integrated_vortex_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[1:2] .+= sample_cache[1] / 2.0

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
    nominal_source_panel_integration(
        integration_options::GaussKronrod,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=false,
    ) -> Matrix{Float64} or (Matrix{Float64}, Float64)

Computes the velocity influence matrix at a control point due to a nominal constant-strength source panel, using Gauss-Kronrod quadrature.

# Description
Numerically integrates the induced velocity from a distributed source along a panel, evaluated at a control point off the panel. The influence is computed over the normalized panel interval `[0, 1]` and reshaped into a `2×2` matrix form.

# Arguments
- `integration_options::GaussKronrod`: Configuration for numerical integration (order, tolerance, max evaluations).
- `node1`: Start point of the panel (2-element vector).
- `node2`: End point of the panel (2-element vector).
- `influence_length`: Effective length used to scale influence computation.
- `controlpoint`: Point where the velocity is being computed (2-element vector).
- `sample_cache`: Temporary storage used during integration to avoid allocation.

# Keyword Arguments
- `debug::Bool=false`: If `true`, also returns the estimated integration error.

# Returns
- `2×2 Matrix{Float64}`: Velocity influence matrix due to the source panel.
- Optionally, the integration error `::Float64` if `debug=true`.
"""
function nominal_source_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return nominal_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    if debug
        return reshape(V, (2, 2)), err
        # return V, err
    else
        return reshape(V, (2, 2))
        # return V
    end
end

"""
    self_source_panel_integration(
        integration_options::GaussKronrod,
        node1,
        node2,
        influence_length,
        controlpoint,
        sample_cache;
        debug=false,
    ) -> Matrix{Float64} or (Matrix{Float64}, Float64)

Computes the self-induced velocity influence matrix at a control point on the same constant-strength source panel, using a combination of numerical and analytical integration techniques.

# Description
This function handles the singularity in self-induced source panel influence by:
- Using Gauss-Kronrod quadrature over a split interval `[0.0, 0.5, 1.0]` for better resolution near singular behavior.
- Adding an analytical correction to the normal component via `analytically_integrated_source_influence`.

The final result is scaled by `influence_length` and returned as a reshaped `2×2` velocity matrix.

# Arguments
- `integration_options::GaussKronrod`: Integration settings including order, tolerance, and maximum evaluations.
- `node1`: Start point of the source panel (2-element vector).
- `node2`: End point of the source panel (2-element vector).
- `influence_length`: Scalar used to scale the integrated influence result.
- `controlpoint`: Point lying on the panel where the velocity is evaluated (2-element vector).
- `sample_cache`: Workspace array used for storing intermediate values, including analytic influence results.

# Keyword Arguments
- `debug::Bool=false`: If `true`, also returns the estimated integration error.

# Returns
- `2×2 Matrix{Float64}`: Reshaped matrix representing the self-induced velocity influence.
- Optionally, the integration error estimate `::Float64` if `debug=true`.
"""
function self_source_panel_integration(
    integration_options::GaussKronrod,
    node1,
    node2,
    influence_length,
    controlpoint,
    sample_cache;
    debug=false,
)

    # Define function to integrate
    function fsample(t)
        return self_source_induced_velocity_sample(
            t, node1, node2, influence_length, controlpoint, sample_cache;
        )
    end

    V, err = quadgk(
        fsample,
        0.0,
        0.5,
        1.0;
        order=integration_options.order,
        maxevals=integration_options.maxevals,
        atol=integration_options.atol,
    )

    sample_cache[1], sample_cache[2] = analytically_integrated_source_influence(
        controlpoint[2], influence_length
    )

    V .*= influence_length
    V[3:4] .+= sample_cache[2] / 2.0

    if debug
        return reshape(V, (2, 2)), err
    else
        return reshape(V, (2, 2))
    end
end
