project_dir = dirname(dirname(dirname(dirname(@__FILE__))))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE

datapath = project_dir * "/test/data/dfdc_init_iter1/"
savepath = project_dir * "/test/manual_tests/dfdc_comparisons/"

include(project_dir * "/visualize/plots_default.jl")

# read in DFDC data
include(datapath * "reformat_dfdc_data.jl")
include(datapath * "ductape_parameters.jl")
include(datapath * "ductape_formatted_dfdc_geometry.jl")

# generate inputs
inputs = dt.precomputed_inputs(
    system_geometry,
    rotor, #vector of named tuples
    freestream,
    reference_parameters;
    debug=false,
)

# TODO: go get steady cp values instead of velocities and see if that compares easier.
_, _, _, ductinner_zpts, ductouter_zpts, hub_zpts = dt.split_bodies(
    ones(inputs.body_vortex_panels.totpanel), inputs.body_vortex_panels; duct=true, hub=true
)

duct_ctrlpt = reverse([ductinner_zpts; ductouter_zpts])
hub_ctrlpt = reverse(hub_zpts)

# - get DuctAPE body data - #
gamb, wakerhs, rotor1rhs = dt.calculate_body_vortex_strengths!(
    inputs.gamb,
    inputs.A_bb,
    inputs.b_bf,
    -gamw0,
    inputs.A_bw,
    # -waic,
    inputs.A_pw,
    # -waic_pcp,
    sigr0, #note, these are zero for initial and first iteration.
    inputs.A_br,
    inputs.A_pr,
    inputs.RHS;
    post=true,
)

pgamb = plot(
    1:length(gamb), gamb; label="DuctAPE", xlabel="eqn #", ylabel="solution", linewidth=2
)
plot!(pgamb, 1:length(gamb), -gamb0; label="-DFDC")
savefig(pgamb, savepath * "initial-linear-solution-comparison.pdf")

pvinf = plot(
    1:length(wakerhs),
    inputs.b_bf;
    label="DuctAPE",
    xlabel="eqn #",
    ylabel="vinf dot n",
    linewidth=2,
)
plot!(pvinf, 1:length(wakerhs), -vinfrhs0; label="-DFDC")
savefig(pvinf, savepath * "initial-vinf-rhs-comparison.pdf")

prhsw = plot(
    1:length(wakerhs),
    wakerhs;
    label="DuctAPE",
    xlabel="eqn #",
    ylabel="vwake dot n",
    linewidth=2,
)
plot!(prhsw, 1:length(wakerhs), -wakerhs0; label="-DFDC")
savefig(prhsw, savepath * "initial-wake-rhs-comparison.pdf")

# - get body tangential velocities - #
vt = dt.get_body_tangential_velocities(inputs, -gamw0, sigr0)

plot([vt.casing_zpts; vt.nacelle_zpts], vt.Vtot_prejump[1, 1:106])
plot!(
    reverse([vt.casing_zpts; vt.nacelle_zpts]),
    iter0_controlpoint_vels.duct_vel[:, 1];
    color=1,
    label="DFDC",
    linestyle=:dash,
)

# include jump term and normalize
casing_cp = dt.steady_cp(vt.vtan_casing, Vinf, Vref)
nacelle_cp = dt.steady_cp(vt.vtan_nacelle, Vinf, Vref)
centerbody_cp = dt.steady_cp(vt.vtan_centerbody, Vinf, Vref)

# - Plot cp comparison - #
pd = plot(;
    xlabel=L"z",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
pc = plot(;
    xlabel=L"z",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(
    pd, duct_ctrlpt, duct_cpR2; label="DFDC initial", color=2, linestyle=:dash, linewidth=2
)
plot!(pd, vt.casing_zpts, casing_cp; color=1, label="DuctAPE initial")
plot!(pd, vt.nacelle_zpts, nacelle_cp; color=1, label="")
savefig(pd, savepath * "initial_duct_cp_comp.pdf")

plot!(pc, hub_ctrlpt, hub_cpR2; label="DFDC initial", color=2, linestyle=:dash, linewidth=2)
plot!(pc, vt.centerbody_zpts, centerbody_cp; color=1, label="DuctAPE initial")
savefig(pc, savepath * "initial_centerbody_cp_comp.pdf")
