include("newton_setup.jl")
include("newton_solve_dt.jl")
# include("test/manual_tests/rotor_wake_tests/newton_setup.jl")

states, params = setup_stuff()
states = solve!(states, params)

# - finalize Plots - #
function add_wake_on_rotor_vm_inducement!(
    vx_rotor, vr_rotor, wake_vortex_strengths, vxd_wake_to_rotor, vrd_wake_to_rotor
)

    ##### ----- Wake on Rotor Influence ----- #####
    # - Loop through the wakes - #
    for i in 1:length(vxd_wake_to_rotor)

        # axial direction
        # vx_wake_i_on_rotor = vxd_wake_to_rotor[i] * wake_vortex_strengths[i, :]
        vx_rotor .+= vxd_wake_to_rotor[i] * wake_vortex_strengths[i, :]

        # radial direction
        # vr_wake_i_on_rotor = vrd_wake_to_rotor[i] * wake_vortex_strengths[i, :]
        vr_rotor .+= vrd_wake_to_rotor[i] * wake_vortex_strengths[i, :]
    end

    return nothing
end

"""
TODO: eventually want to add upstream rotor contributions as well.
"""
function add_rotor_self_induced_vtheta!(Wtheta, radial_positions, BGamma)

    # the rotor adds half of its own circulation to the self-induced tangential velocity
    Wtheta .+= 1.0 ./ (4.0 .* pi .* radial_positions) .* BGamma

    return nothing
end

Gamma = states[params.Gammaidx]
display(Gamma)
gamma_theta = states[params.gamma_theta_idx]
display(gamma_theta)

wake_vortex_strengths = repeat(gamma_theta; inner=(1, params.nxwake))

rpc = params.rotor_panel_centers

TF = eltype(Gamma)

# if TF == Float64
#     println("input Gamma")
#     println(Gamma)
# end

vx_rotor = zeros(TF, length(Gamma))
vr_rotor = zeros(TF, length(Gamma))
vtheta_rotor = zeros(TF, length(Gamma))

# - add the wake induced velocities at the rotor plane to Vm - #
add_wake_on_rotor_vm_inducement!(
    vx_rotor,
    vr_rotor,
    wake_vortex_strengths,
    params.vx_rw,
    params.vr_rw,
)

# - add the rotor rotational self-induction directly - #
add_rotor_self_induced_vtheta!(vtheta_rotor, rpc, params.num_blades * Gamma)

# - Get the blade element reference frame total velocity components - #
# the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
Wx_rotor = vx_rotor .+ params.Vinf
# the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
Wtheta_rotor = vtheta_rotor .- params.Omega .* rpc

Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

# - Get the inflow magnitude at the rotor as the combination of all the components - #
Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

# calculate angle of attack
# alpha = dt.calculate_angle_of_attack(params.twist, Wm_rotor, Wtheta_rotor)
alpha = params.twist .- atan.(Wm_rotor, -Wtheta_rotor)

# - Look up lfit and drag data - #
cl = zeros(TF, length(alpha))
cd = zeros(TF, length(alpha))
for a in 1:length(alpha)
    cl[a], cd[a] = dt.search_polars(params.af, alpha[a])
end

# - get average meridional velocity at rotor - #
Vmavg = [
    0.5 * (params.Vinf + Wm_rotor[1])
    0.5 .* (Wm_rotor[1:(end - 1)] .+ Wm_rotor[2:end])
    0.5 * (Wm_rotor[end] + params.Vinf)
]

plot!(
    params.pvm,
    Wm_rotor,
    params.rotor_panel_centers,
    linestyle=:dash,
    label="Vm from wake at panel centers, converged",
)
plot!(
    params.pvm,
    Vmavg,
    params.rotor_panel_edges;
    label="Vm_avg at panel edges, converged",
)
plot!(
    params.pgw, gamma_theta, params.rotor_panel_edges; label="f(1/Vmavg), converged"
)
plot!(params.pg, Gamma, rpc; label="converged")
plot!(params.pcd, cd, rpc; label="converged")
plot!(params.pcl, cl, rpc; label="converged")
plot!(params.pa, alpha, rpc; label="converged")
plot!(params.pW, Wmag_rotor, rpc; label="converged")
plot!(params.pv, vtheta_rotor, rpc; label="converged")
plot!(params.pu, vx_rotor, rpc; label="converged")

savefig(params.pu, "test/manual_tests/rotor_wake_tests/newton_vx.pdf")
savefig(params.pv, "test/manual_tests/rotor_wake_tests/newton_vt.pdf")
savefig(params.pW, "test/manual_tests/rotor_wake_tests/newton_W.pdf")
savefig(params.pa, "test/manual_tests/rotor_wake_tests/newton_aoa.pdf")
savefig(params.pcl, "test/manual_tests/rotor_wake_tests/newton_cl.pdf")
savefig(params.pcd, "test/manual_tests/rotor_wake_tests/newton_cd.pdf")
savefig(params.pg, "test/manual_tests/rotor_wake_tests/newton_circ.pdf")
savefig(params.pgw, "test/manual_tests/rotor_wake_tests/newton_wakegamma.pdf")
savefig(params.pvm, "test/manual_tests/rotor_wake_tests/newton_vms.pdf")

