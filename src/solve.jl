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
 """
function solve_wrapper(x_init, params)

    # - Wrap Residual Function to allow for parameters - #
    rwrap(F, x_init) = residual!(F, x_init, params)

    # - Call NLsolve function using AD - #
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
function run(duct_coordinates, hub_coordinates)

    #TODO: Initialize. Calculate the initial guesses for Gamma and Sigma from the inputs.
    #TODO: Assemble the parameters that need to be passed into the solver.  for example, all the geometry and coefficient matrices.

    # - Run solver to find Gamma and Sigma values - #
    # params is a tuple
    GammaSigma = ImplicitAD.implicit(solve_wrapper, residual!, GammaSigma_init, params)

    return GammaSigma, converged

    #TODO: do the rest.  need to use the converged gamma and sigma values to go through and get all the singularity strengths required for post-processing.

end

"""
"""
function update_gamma_sigma(GammaSigma_init, params)

    # - Calculate Enthalpy Jumps - #

    # - Calculate Net Circulation - #

    # - Calculate Meridional Velocities - #

    # - Calculate Wake Vorticity - #

    # - Solve Full Linear System - #

    # - Calculate Induced Velocities at Rotors - #

    # - Calculate Blade Element Angles of Attack - #

    # - Calculate Inflow Velocities at Blade Elements - #

    # - Look up Blade Element Polar Data - #

    # - Calculate Updated Circulation and Source Strengths - #

    # - Return Updated Circulation and Source Strengths in Single Vector - #
    return nothing
end

