#=

Solver Functions

Authors: Judd Mehr,

=#

"""
    residual!(gamma_sigma_residuals, gamma_sigma, problem)

Takes in the vector of Gamma and Sigma inputs, their residuals, and the problem object; calculates the residuals of Gamma and Sigma.

**Arguments:**
- `gamma_sigma_residuals::Array{Float}` : residuals of Gamma and Sigma
- `gamma_sigma::Array{Float}` : Gamma and Sigma input values
- `problem::DuctTAPE.Problem` : object containing all the parameters, geometry, and coefficients required to obtain the residuals

**Returns:**
Nothing, gamma_sigma_residuals and problem are updated as required.
"""
function residual!(gamma_sigma_residuals, gamma_sigma, problem)

    # Unpack Gamma and Sigma inputs
    Gamma = reshape(gamma_sigma[1:(system.ngamma)], system.shapegamma)
    Sigma = reshape(gamma_sigma[(system.ngamma + 1):end], system.shapesigma)

    # Get wake panel strengths
    update_wake!(Gamma, problem)

    # Solve linear system for wall panel strengths
    update_wall_panel_strengths(Gamma, Sigma, problem)

    # Update flow velocities
    update_velocities!(problem)

    # Re-run CCBlade with updated flow inputs
    Ws, cls, cds = run_ccblade(problem.rotors, problem.rotor_velocities)

    #Get Gamma and Sigma Residuals (from eqn 73 and 75)
    gamma_sigma_residuals = calculate_residuals(Gamma, Sigma, Ws, cls, cds, problem.blades)

    return nothing
end

"""
    solve!(problem)

Solves ducted rotor system defined by Problem object.

**Arguments:**
- `problem::DuctTAPE.System` : ducted rotor system problem object

**Returns:**
- `solution::DuctTAPE.System` : updated problem object
- `convergence::Bool` : convergence flag
"""
function solve!(problem)
    gamma_sigma_init = [problem.Gammas; problem.Sigmas]

    function res_wrap!(gamma_sigma_residuals, GS)
        return residual!(gamma_sigma_residuals, GS, problem)
    end

    #use nlsolve to solve newton system.
    solver_results = NLsolve.nlsolve(res_wrap!, gamma_sigma_init; autodiff=:forward)

    #=
    solver_results includes:

    method::String
    initial_x::I
    zero::Z
    residual_norm::rT
    iterations::Int
    x_converged::Bool
    xtol::rT
    f_converged::Bool
    ftol::rT
    trace::SolverTrace
    f_calls::Int
    g_calls::Int
    =#

    #create solution object from problem and solver_results

    #return updated problem and convergence flag
    return solution, NLsolve.converged(solver_results)
end

"""
"""
function calculate_residuals(Gamma, Sigma, Ws, cls, cds, blades)

    #rename for convenience
    #number of elements, number of rotors
    ngel, nr = size(Gamma)
    nsel, _ = size(Sigma)

    #initialize
    gamma_res = [0.0 for i in 1:Gamma[:, 1], j in 1:Gamma[1, :]]
    sigma_res = [0.0 for i in 1:Sigma[:, 1], j in 1:Sigma[1, :]]

    #Loop through radial stations calculating Gamma Values
    for i in 1:ngel
        for j in 1:nr

            #calculate gamma value (function defined in rotors.jl)
            gamma_calc = blade_section_gamma(Ws[i, j], blades[j].cdim[i], cls[i, j])

            #get difference between input and calculated gamma values
            gamma_res[i, j] = gamma_calc - Gamma[i, j]
        end
    end

    #Loop through radial stations, applying averages to calculate Sigma values
    for i in 1:nsel
        for j in 1:nr
            #get averaged values
            W = 0.5 * (Ws[i, j] + Ws[i + 1, j])
            chord = 0.5 * (blades[j].cdim[i] + blades[j].cdim[i + 1])
            cd = 0.5 * (cds[i, j] + cds[i + 1, j])
            r = 0.5 * (blades[j].rdim[i] + blades[j].rdim[i + 1])

            #calculate sigma value (function defined in rotors.jl)
            sigma_calc = blade_section_sigma(W, chord, cd, blades[j].B, r)

            #get difference between input and calculated sigmas
            sigma_res[i, j] = sigma_calc - Sigma[i, j]
        end
    end

    #reformat residuals into single array (matching gamma_sigma input for clarity)
    return [reshape(gamma_res, (:, 1)); reshape(sigma_res, (:, 1))]
end
