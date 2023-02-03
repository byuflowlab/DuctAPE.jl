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

    # - Assemble GammaSigma as a vector - #
    GammaSigma_init = [body_vortex_strengths[1:(end - 1)]; Gamma_init] #; Sigma_init]

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

    # - Rename for Convenience - #
    nbe = length(params.blade_elements[1].radial_positions)
    nr = params.num_rotors
    nbp = params.num_body_panels
    gamma_idx = (nbp + 1):(nbp + nbe * nr)
    sigma_idx = (nbp + 1 + nbe * nr):(nbp + 2 * nbe * nr)

    #if only using Gamma at first do this
    body_vortex_strengths = view(Gamma_Sigma, 1:nbp)
    Gamma = reshape(view(Gamma_Sigma, gamma_idx), (nbe, nr))
    Sigma = similar(Gamma)
    # Sigma = reshape(view(Gamma_Sigma, sigma_idx), (nbe, nr))

    #---------------------------------#
    #             UPDATES             #
    #---------------------------------#

    # - Calculate Enthalpy Jumps - #
    H_tilde = calculate_enthalpy_jumps(Gamma, params.blade_elements)

    # - Calculate Net Circulation - #
    BGamma, Gamma_tilde = calculate_net_circulation(Gamma, params.blade_elements)

    # - Get Surface Velocity on Duct Inner Surface - #
    #ASSUMES DUCT IS FIRST BODY
    vs = get_surface_velocity(
        body_vortex_strengths, params.body_panels[1], params.wake_panels
    )

    # - Calculate Meridional Velocities - #
    vm = calculate_wake_velocities(
        params.wake_panels[1].panel_center[:, 1],
        vs,
        params.rotoridxs,
        Gamma_tilde,
        H_tilde,
        params.blade_elements,
    )

    # - Calculate Wake Vorticity - #
    wake_gammas = calculate_wake_vorticity(vm)

    # - Calculate Induced Velocities at Rotors - #
    vi = calculate_induced_velocities(
        BGamma,
        Gamma_tilde,
        params.blade_elements,
        params.A_body_to_rotor,
        body_vortex_strengths,
        params.A_wake_to_rotor,
        wake_gammas,
        # params.A_rotor_to_rotor,
        # Sigma,
    )

    # - Solve Full Linear System - #
    update_body_vortex_strengths!(
        body_vortex_strengths,
        params.A_body_to_body,
        params.bc_freestream_to_body,
        wake_gammas,
        params.A_wake_to_body,
    )#, Sigma, params.A_rotor_to_body)

    # - Calculate Updated Circulation and Source Strengths - #
    calculate_gamma_sigma!(
        Gamma, Sigma, params.blade_elements, params.freestream.Vinf, vi.vm, vi.vtheta
    )

    return nothing
end
