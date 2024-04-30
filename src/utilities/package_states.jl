"""
    extract_initial_guess(
        solver_options::SolverOptionsType, sensitivity_parameters, state_dims
    )

Extract initial guess from the solve parameters cache vector.

# Arguments
- `solver_options::SolverOptionsType` : used for dispatch
- `sensitivity_parameters::Vector{Float}` : vector form of solve parameter cache passed into the solver.
- `state_dims::NamedTuple` : dimensions and indices of state variables within the solve parameter cache vector

# Returns
- initial_guess::Vector{Float}` : a vector of the solver initial guess
"""
function extract_initial_guess(
    solver_options::TS, sensitivity_parameters, state_dims
) where {TS<:ExternalSolverOptions}
    return view(
        sensitivity_parameters, state_dims.vz_rotor.index[1]:state_dims.Cm_wake.index[end]
    )
end

function extract_initial_guess(
    solver_options::CSORSolverOptions, sensitivity_parameters, state_dims
)
    return view(sensitivity_parameters, state_dims.Gamr.index[1]:state_dims.gamw.index[end])
end

"""
    extract_state_variables(solver_options::SolverOptionsType, vars, dims)

Reshape the state variables from a single vector, to multiple arrays.

# Arguments

# Returns if solver_options <: CSORSolverOptions
- `Gamr::type` : Blade element circulation strengths
- `sigr::type` : Rotor source panel strengths
- `gamw::type` : Wake vortex panel strengths

# Returns if solver_options <: Union{ExternalSolverOptions, PolyAlgorithmOptions}
- `vz_rotor::Vector{Float}` : axial induced rotor velocity state container
- `vtheta_rotor::Vector{Float}` : tangential induced rotor velocity state container
- `Cm_wake::Vector{Float}` : absolute meridional wake control point velocity state container
"""
function extract_state_variables(
    solver_options::TS, vars, dims
) where {TS<:Union{ExternalSolverOptions,PolyAlgorithmOptions}}

    # - Separate out - #
    vz_rotor = @views reshape(vars[dims.vz_rotor.index], dims.vz_rotor.shape)
    vtheta_rotor = @views reshape(vars[dims.vtheta_rotor.index], dims.vtheta_rotor.shape)
    Cm_wake = @views reshape(vars[dims.Cm_wake.index], dims.Cm_wake.shape)

    return vz_rotor, vtheta_rotor, Cm_wake
end

function extract_state_variables(solver_options::CSORSolverOptions, vars, dims)

    # - Separate out - #
    Gamr = @views reshape(vars[dims.Gamr.index], dims.Gamr.shape)
    sigr = @views reshape(vars[dims.sigr.index], dims.sigr.shape)
    gamw = @views reshape(vars[dims.gamw.index], dims.gamw.shape)

    return Gamr, sigr, gamw
end
