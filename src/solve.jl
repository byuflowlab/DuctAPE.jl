#=

Overall solver functions.

Authors: Judd Mehr,

=#

"""
This is the function you run to actually solve stuff
"""
function analyze_propulsor(duct_coordinates, hub_coordinates, rotor_parameters, freestream)

    # - Initialize Parameters (general parameters, geometries, meshes, coefficient matrics, etc.) - #
    params = initialize_parameters(
        duct_coordinates, hub_coordinates, rotor_parameters, freestream
    )

    # Solve the no-rotor problem
    body_vortex_strengths = ImplicitAD.implicit_linear(
        params.A_body_to_body, params.bc_freestream_to_body
    )

    # - Calculate the initial guesses for Gamma and Sigma from the inputs. - #
    # Uses the no-rotor solution, setting induced velocities to zero in the wake and on the rotors.
    Gamma_init, Sigma_init = calculate_gamma_sigma(
        params.blade_elements, params.freestream.Vinf
    )

    # - Calculate initial guesses for wake vortex strengths - #
    wake_vortex_strengths = initialize_wake_velocities()

    # - Assemble GammaSigma as a vector - #
    GammaSigma_init = [
        body_vortex_strengths[1:(end - 1)]
        reshape(vcat, wake_vortex_strengths)
        Gamma_init
        # Sigma_init]
    ]

    # - Run solver to find Gamma and Sigma values - #
    GammaSigma = ImplicitAD.implicit(solve!, residual!, GammaSigma_init, params)

    return GammaSigma, params.converged

    #TODO: do the rest.  need to use the converged gamma and sigma values to go through and get all the singularity strengths required for post-processing.

end

"""
This is the function being solved
GammaSigma_init are the inital guess for the gamma and sigma values.
F are the output gamma and sigma values minus the input ones. (this is what we want to drive to zero).
"""
function residual!(F, GammaSigma_init, params)
    GammaSigma_new = deepcopy(GammaSigma_init)

    # - Calculated Updated Gamma and Sigma Values - #
    update_gamma_sigma!(GammaSigma_new, params)

    # - Return Difference in Gamma and Sigma Values - #
    @. F = GammaSigma_init - GammaSigma_new

    return nothing
end

"""
 This function wraps the residual function in order to allow for additional parameters as inputs
 params.converged is updated in place in this function.
 """
function solve!(x_init, params)

    # - Define closure that allows for parameters - #
    rwrap(F, x_init) = residual!(F, x_init, params)

    # - Call NLsolve function using AD for Jacobian - #
    #= res is of type NLsolve.SolverResults.
    The zero field contains the "solution" to the non-linear solve.
    The converged() function tells us if the solver converged.
    =#
    res = NLsolve.nlsolve(
        rwrap,
        x_init;
        autodiff=:forward,
        method=:newton,
        linesearch=BackTracking(; maxstep=1e6),
    )

    # - Overwrite the convergence flag in the parameters - #
    params.converged[1] = converged(res)

    # - Return the values driving the residual to zero
    return res.zero
end

"""
"""
function update_gamma_sigma!(Gamma_Sigma, params)

    #---------------------------------#
    #              SET UP             #
    #---------------------------------#

    # println("\nBegin Iteration\n")

    body_vortex_strengths, wake_vortex_strengths, rotor_circulation_strengths, blade_element_source_strengths = extract_state_variables(
        Gamma_Sigma, params
    )

    rotor_panel_source_strengths =
        (
            blade_element_source_strengths[1:(end - 1)] .+
            blade_element_source_strengths[2:end]
        ) / 2.0

    # - Calculate Enthalpy Jumps - #
    H_tilde = calculate_enthalpy_jumps(rotor_circulation_strengths, params.blade_elements)

    # - Calculate Net Circulation - #
    BGamma, Gamma_tilde = calculate_net_circulation(
        rotor_circulation_strengths, params.blade_elements
    )

    # - Calculate Induced Velocities at Rotors - #
    vi_rotor = calculate_induced_velocities(
        BGamma,
        Gamma_tilde,
        params.blade_elements,
        params.A_body_to_rotor,
        body_vortex_strengths,
        params.A_wake_to_rotor,
        wake_vortex_strengths,
        # params.A_rotor_to_rotor,
        # rotor_panel_source_strengths,
    )

    # - Calculate Meridional Velocities on Wakes - #
    wake_velocities = calculate_wake_velocities(
        A_body_to_wake,
        gamma_body,
        A_wake_to_wake,
        gamma_wake,
        # A_rotor_to_wake,
        # rotor_panel_source_strengths,
    )

    #---------------------------------#
    #             UPDATES             #
    #---------------------------------#

    # - Update Wake Vortex Strengths - #
    calculate_wake_vorticity!(wake_vortex_strengths, wake_velocities)

    # - Update Body Vortex Strengths - #
    calculate_body_vortex_strengths!(
        body_vortex_strengths,
        params.A_body_to_body,
        params.bc_freestream_to_body,
        wake_vortex_strengths,
        params.A_wake_to_body,
    )#, blade_element_source_strengths, params.A_rotor_to_body)

    # - Calculate Updated Circulation and Source Strengths - #
    calculate_gamma_sigma!(
        rotor_circulation_strengths,
        blade_element_source_strengths,
        params.blade_elements,
        params.freestream.Vinf,
        vi_rotor.vm,
        vi_rotor.vtheta,
    )

    return nothing
end

"""
"""
function extract_state_variables(Gamma_Sigma, params)

    # - Rename for Convenience - #
    nbe = length(params.blade_elements[1].radial_positions)
    nr = params.num_rotors
    nx = params.num_wake_x_panels
    nbp = params.num_body_panels
    n_g_bw = nbp + nbe * nx #number of gammas for bodies and wake
    ng = n_g_bw + nbe * nr #number of gammas and Gammas

    #  Body gamma Indices
    body_idx = 1:nbp

    # Wake gamma_theta Indices
    wake_gamma_idx = (nbp + 1):(nbp + nbe * nx)

    # Rotor Gamma Indices
    rotor_Gamma_idx = (n_g_bw + 1):(ng)

    # Rotor Sigma Indices
    rotor_Sigma_idx = (ng + 1):(ng + 1 + nbe * nr)

    # - Extract State Variables - #
    body_vortex_strengths = view(Gamma_Sigma, body_idx)
    wake_vortex_strengths = view(Gamma_Sigma, wake_gamma_idx, (nbe, nx))
    rotor_circulation_strengths = reshape(view(Gamma_Sigma, rotor_Gamma_idx), (nbe, nr))
    blade_element_source_strengths = similar(rotor_circulation_strengths)
    # blade_element_source_strengths = reshape(view(Gamma_Sigma, rotor_Sigma_idx), (nbe, nr))

    return body_vortex_strengths,
    wake_vortex_strengths, rotor_circulation_strengths,
    blade_element_source_strengths
end
