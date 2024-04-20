"""
"""
function extract_initial_guess(
    solver_options::TS, sensitivity_parameters, state_dims
) where {TS<:ExternalSolverOptions}
    return view(
        sensitivity_parameters, state_dims.vz_rotor.index[1]:state_dims.Cm_wake.index[end]
    )
end

"""
"""
function extract_initial_guess(
    solver_options::CSORSolverOptions, sensitivity_parameters, state_dims
)
    return view(sensitivity_parameters, state_dims.Gamr.index[1]:state_dims.gamw.index[end])
end

"""
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

"""
"""
function extract_state_variables(solver_options::CSORSolverOptions, vars, dims)

    # - Separate out - #
    Gamr = @views reshape(vars[dims.Gamr.index], dims.Gamr.shape)
    sigr = @views reshape(vars[dims.sigr.index], dims.sigr.shape)
    gamw = @views reshape(vars[dims.gamw.index], dims.gamw.shape)

    return Gamr, sigr, gamw
end
