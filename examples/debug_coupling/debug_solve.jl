using NLsolve

p = (converged=[false],)
rwrap(r, states) = dt.residual!(r, states, inputs, p)

res = NLsolve.nlsolve(
    rwrap,
    initial_states;
    autodiff=:forward,
    method=:newton,
    iterations=25,
    show_trace=true,
    # linesearch=BackTracking(; maxstep=maximum_linesearch_step_size),
)

strengths = res.zero

gamb, gamw, Gamr, sigr = dt.extract_state_variables(strengths, inputs)
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

plot!(pb, dp, gamd; xlabel="x", ylabel="cp", label="converged duct surface pressure")
# plot!(pb, hp, gamh ; label="iter hub surface pressure")

## -- check rotor circulation and source initial strengths -- ##
plot!(pG, Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="converged")

plot!(ps, sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="converged")

plot!(
    pw,
    gamw,
    inputs.rotor_panel_edges;
    xlabel=L"\gamma_\theta",
    ylabel="r",
    label="converged",
)



include("../rotor_only/run_ccblade.jl")
ccbouts = run_ccblade(Vinf)
plot!(pG, ccbouts.circ, rnondim*Rtip; xlabel=L"\Gamma", ylabel="r", label="CCBlade")


savefig(pb, "examples/debug_coupling/body-strengths-converged.pdf")
savefig(pG, "examples/debug_coupling/circulation-converged.pdf")
savefig(ps, "examples/debug_coupling/sources-converged.pdf")
savefig(pw, "examples/debug_coupling/wake-strengths-converged.pdf")

