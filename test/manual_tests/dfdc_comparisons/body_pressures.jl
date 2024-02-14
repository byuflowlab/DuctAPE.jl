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
    rotorstator_parameters, #vector of named tuples
    freestream,
    reference_parameters;
    debug=false,
)

_, _, _, ductinner_zpts, ductouter_zpts, hub_zpts = dt.split_bodies(
    ones(inputs.body_vortex_panels.totpanel), inputs.body_vortex_panels; duct=true, hub=true
)

duct_ctrlpt = reverse([ductinner_zpts; ductouter_zpts])
hub_ctrlpt = reverse(hub_zpts)

##### ----- Compare full Cps ----- #####
Gamr = zeros(length(Gamr1), 1) .= copy(Gamr1)
sigr = zeros(length(sigr2), 1) .= copy(sigr2)
gamw = -copy(gamw1)
outs = dt.post_process(Gamr, sigr, gamw, inputs; write_outputs=false)

pcpd = plot(;
    xlabel=L"z",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pcpd, duct_ctrlpt, duct_cpdel2; linewidth=2, linestyle=:dash, color=2, label="DFDC")
plot!(pcpd, outs.bodies.casing_zpts, outs.bodies.cp_casing_out; color=1, label="DuctAPE")
plot!(pcpd, outs.bodies.nacelle_zpts, outs.bodies.cp_nacelle_out; color=1, label="")
savefig(pcpd, savepath * "initial_duct_full_cp_comp.pdf")

pcpc = plot(;
    xlabel=L"z",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(pcpc, hub_ctrlpt, hub_cpdel2; linewidth=2, linestyle=:dash, color=2, label="DFDC")
plot!(
    pcpc, outs.bodies.centerbody_zpts, outs.bodies.cp_centerbody_out; color=1, label="DuctAPE"
)
savefig(pcpc, savepath * "initial_centerbody_full_cp_comp.pdf")

##### ----- Compare Steady Cps ----- #####

# - get stead pressure coefficients manuall - #
# note that gamb inside inputs was updated in the call to the post-process function above.
vt = dt.get_body_tangential_velocities(inputs, -gamw0, sigr0)

casing_cp = dt.steady_cp(vt.vtan_casing_out, Vinf, Vref)
nacelle_cp = dt.steady_cp(vt.vtan_nacelle_out, Vinf, Vref)
centerbody_cp = dt.steady_cp(vt.vtan_centerbody_out, Vinf, Vref)

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
    pd, duct_ctrlpt, duct_cpR0; label="DFDC initial", color=2, linestyle=:dash, linewidth=2
)
plot!(pd, vt.casing_zpts, casing_cp; color=1, label="DuctAPE initial")
plot!(pd, vt.nacelle_zpts, nacelle_cp; color=1, label="")
savefig(pd, savepath * "initial_duct_cp_comp.pdf")

plot!(pc, hub_ctrlpt, hub_cpR0; label="DFDC initial", color=2, linestyle=:dash, linewidth=2)
plot!(pc, vt.centerbody_zpts, centerbody_cp; color=1, label="DuctAPE initial")
savefig(pc, savepath * "initial_centerbody_cp_comp.pdf")
