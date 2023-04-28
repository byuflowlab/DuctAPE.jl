# - Initialize body vortex strengths (rotor-off linear problem) - #
gamb = dt.solve_body_system(inputs.A_bb, inputs.b_bf, inputs.kutta_idxs) # get circulation strengths from solving body to body problem

# - Initialize blade circulation and source strengths (assume open rotor) - #

# get floating point type
TF = promote_type(
    eltype(inputs.blade_elements[1].chords),
    eltype(inputs.blade_elements[1].twists),
    eltype(inputs.blade_elements[1].Omega),
    eltype(gamb),
)

# get problem dimensions (number of radial stations x number of rotors)
nr = length(inputs.blade_elements[1].rbe)
nrotor = length(inputs.blade_elements)

# initialize outputs
Gamr = zeros(TF, nr, nrotor)
sigr = zeros(TF, nr, nrotor)
gamw = zeros(TF, nr + 1, nrotor)
wake_vortex_strengths = repeat(gamw; inner=(1, inputs.num_wake_x_panels))

vx_rotor, vr_rotor, vtheta_rotor = dt.calculate_induced_velocities_on_rotors(
    inputs.blade_elements,
    Gamr,
    inputs.vx_rw,
    inputs.vr_rw,
    wake_vortex_strengths,
    inputs.vx_rr,
    inputs.vr_rr,
    sigr,
    inputs.vx_rb,
    inputs.vr_rb,
    gamb,
)
# vx_rotor = zeros(size(sigr))
# vr_rotor = zeros(size(sigr))
# vtheta_rotor = zeros(size(sigr))

# the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
Wx_rotor = vx_rotor .+ inputs.Vinf

# the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* inputs.rotor_panel_centers

# meridional component
Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

# Get the inflow magnitude at the rotor as the combination of all the components
Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

# initialize circulation and source panel strengths
dt.calculate_gamma_sigma!(
    Gamr, sigr, inputs.blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor
)

# - Calculate net circulation and enthalpy jumps - #
Gamma_tilde = dt.calculate_net_circulation(Gamr, inputs.blade_elements.B)
H_tilde = dt.calculate_enthalpy_jumps(
    Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
)

# - update wake strengths - #
dt.calculate_wake_vortex_strengths!(
    gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
)

# # try initilizing to freestream instead...
# gamw = dt.initialize_wake_vortex_strengths(
#     inputs.Vinf,
#     Gamr,
#     inputs.blade_elements.Omega,
#     inputs.blade_elements.B,
#     inputs.rotor_panel_edges,
# )

# - Combine initial states into one vector - #
initial_states = vcat(
    gamb,               # body vortex panel strengths
    reduce(vcat, gamw), # wake vortex sheet strengths
    reduce(vcat, Gamr), # rotor circulation strengths
    reduce(vcat, sigr), # rotor source panel strengths
)

gamb, gamw, Gamr, sigr = dt.extract_state_variables(initial_states, inputs)

## -- check body surface velocity initialiation -- ##
dp = inputs.body_panels[1].panel_center[:, 1]
# hp = inputs.body_panels[2].panel_center[:, 1]
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

pb = plot(
    dp, gamd; yflip=true, xlabel="x", ylabel="cp", label="initial duct surface pressure"
)
# plot!(pb, hp, gamh ; label="initial hub surface pressure")
savefig(pb, "examples/debug_coupling/initialize-body-vortex-state.pdf")

## -- check rotor circulation and source initial strengths -- ##
pG = plot(Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="initial")
savefig(pG, "examples/debug_coupling/initialize-rotor-circulation-state.pdf")

ps = plot(sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="initial")
savefig(ps, "examples/debug_coupling/initialize-rotor-source-state.pdf")

pw = plot(
    gamw, inputs.rotor_panel_edges; xlabel=L"\gamma_\theta", ylabel="r", label="initial"
)
savefig(pw, "examples/debug_coupling/initialize-wake-vortex-state.pdf")

pvx = plot(vx_rotor, inputs.rotor_panel_centers; xlabel=L"v_x", ylabel="r", label="initial")
savefig("examples/debug_coupling/vx-for-initial-states.pdf")

pvr = plot(vr_rotor, inputs.rotor_panel_centers; xlabel=L"v_r", ylabel="r", label="initial")
savefig("examples/debug_coupling/vr-for-initial-states.pdf")

pvt = plot(
    vtheta_rotor,
    inputs.rotor_panel_centers;
    xlabel=L"v_\theta",
    ylabel="r",
    label="initial",
)
savefig("examples/debug_coupling/vtheta-for-initial-states.pdf")

pwx = plot(Wx_rotor, inputs.rotor_panel_centers; xlabel=L"W_x", ylabel="r", label="initial")
savefig("examples/debug_coupling/wx-for-initial-states.pdf")

pwt = plot(
    Wtheta_rotor,
    inputs.rotor_panel_centers;
    xlabel=L"W_\theta",
    ylabel="r",
    label="initial",
)
savefig("examples/debug_coupling/wtheta-for-initial-states.pdf")

pwm = plot(Wm_rotor, inputs.rotor_panel_centers; xlabel=L"W_m", ylabel="r", label="initial")
savefig("examples/debug_coupling/wm-for-initial-states.pdf")
