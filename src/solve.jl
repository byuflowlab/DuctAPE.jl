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
function analyze_propulsor(duct_coordinates, hub_coordinates, rotor, stator, freestream)

    # - Initialize Parameters (geometries, meshes, coefficient matrics) - #
    params = initialize_parameters(
        duct_coordinates, hub_coordinates, rotor, stator, freestream
    )

    # - Calculate the initial guesses for Gamma and Sigma from the inputs. - #
    Gamma_rotor_init, Sigma_rotor_init = calculate_gamma_sigma(
        Vinf, params.rotor_blade_elements; vm=0.0, vtheta=0.0
    )
    Gamma_stator_init, Sigma_stator_init = calculate_gamma_sigma(
        Vinf, params.stator_blade_elements; vm=0.0, vtheta=0.0
    )

    # - Assemble GammaSigma as a vector - #
    GammaSigma_init - [Gamma_rotor_init; Gamma_stator_init] #; Sigma_rotor_init; Sigma_stator_init]

    # - Run solver to find Gamma and Sigma values - #
    # params is a tuple
    GammaSigma = ImplicitAD.implicit(solve!, residual!, GammaSigma_init, params)

    return GammaSigma, converged

    #TODO: do the rest.  need to use the converged gamma and sigma values to go through and get all the singularity strengths required for post-processing.

end

"""
"""
function update_gamma_sigma(GammaSigma_init, params)

    # - Rename for Convenience - #
    nbe = length(params.rotor_blade_elements.radial_positions)
    nr = params.num_rotors

    #TODO need to figure out how to do things.  should the params be vectors of stuff rather than splitting out rotor and stator? (yes) need to make updates to the setup structure.
    #TODO: before updating setup structure, figure out how you want THIS function to work, so that you know what the inputs should be and you only have to set it up once.

    #if only using gammas at first do this
    Gammas = reshape(GammaSigma_init, (nbe, nr))

    #TODO: do this when adding in sigmas
    ##first half of the vector are the gammas
    #Gammas = reshape(view(GammaSigma_init, 1:(nbe * nr / 2), 1), (nbe, nr))

    ##second half of the vetor are the sigmas
    #Sigmas = reshape(view(GammaSigma_init, (nbe * nr / 2 + 1):(nbe * nr), 1), (nbe, nr))

    # - Calculate Enthalpy Jumps - #
    delta_h = calculate_enthalpy_jumps(Gammas)

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

