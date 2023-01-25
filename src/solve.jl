#=

Overall solver functions.

Authors: Judd Mehr,

=#

"""
This is the function being solved
GammaSigma_init are the inital guess for the gamma and sigma values.
F are the output gamma and sigma values minus the input ones. (this is what we want to drive to zero).
"""
function residual!(F, GammaSigma_init, params)

    # - Calculated Updated Gamma and Sigma Values - #
    GammaSigma_new = update_gamma_sigma(GammaSigma_init, params)

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
    res = NLsolve.nlsolve(rwrap, x_init; autodiff=:forward)

    # - Overwrite the convergence flag in the parameters - #
    params.converged[1] = converged(res)

    # - Return the values driving the residual to zero
    return res.zero
end

"""
This is the function you run to actually solve stuff
"""
function run(duct_coordinates, hub_coordinates, rotor, stator, freestream)

    #TODO: Assemble the parameters that need to be passed into the solver.  for example, all the geometry and coefficient matrices.
    params = initialize_parameters(
        duct_coordinates, hub_coordinates, rotor, stator, freestream
    )
    #TODO: Initialize. Calculate the initial guesses for Gamma and Sigma from the inputs.

    # - Run solver to find Gamma and Sigma values - #
    # params is a tuple
    GammaSigma = ImplicitAD.implicit(solve!, residual!, GammaSigma_init, params)

    return GammaSigma, converged

    #TODO: do the rest.  need to use the converged gamma and sigma values to go through and get all the singularity strengths required for post-processing.

end

"""
"""
function update_gamma_sigma(GammaSigma_init, params)

    # - Calculate Enthalpy Jumps - #
    delta_h = calculate_enthalpy_jumps()

    # - Calculate Net Circulation - #
    Gamma_tilde = calculate_net_circulation()

    # - Calculate Meridional Velocities - #
    vm = calculate_meridional_velocities()

    # - Calculate Wake Vorticity - #
    wake_gammas = calculate_wake_vorticity()

    # - Solve Full Linear System - #
    body_vortex_strengths = solve_linear_system()

    # - Calculate Induced Velocities at Rotors - #
    vi = calculate_induced_velocities()

    # - Calculate Blade Element Angles of Attack - #
    alpha = calculate_angles_of_attack()

    # - Calculate Inflow Velocities at Blade Elements - #
    W = calculate_inflow_velocities()

    # - Look up Blade Element Polar Data - #
    cl, cd = search_polars()

    # - Calculate Updated Circulation and Source Strengths - #
    Gamma, Sigma = calculate_gamma_sigma()

    # - Return Updated Circulation and Source Strengths in Single Vector - #
    return [Gamma; Sigma]
end

