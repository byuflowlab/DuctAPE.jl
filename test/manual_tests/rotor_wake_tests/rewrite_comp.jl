# include("newton_setup.jl")
# include("newton_solve_dt.jl")
include("test/manual_tests/rotor_wake_tests/newton_setup.jl")
include("test/manual_tests/rotor_wake_tests/newton_solve_dt.jl")

states, params = setup_stuff()

##### ----- OLD FUNCTIONS ----- #####

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

##### ----- OLD SCRIPT STUFF ----- #####

Gamma = view(states, params.Gammaidx)
gamma_theta = view(states, params.gamma_theta_idx)

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
    vx_rotor, vr_rotor, wake_vortex_strengths, params.vx_rw, params.vr_rw
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

# - calculate circulation and associated jumps - #
Gamma_updated = 0.5 .* Wmag_rotor .* params.chord .* cl

# if TF == Float64
#     println("updated Gamma_updated")
#     println(states[params.Gammaidx])
# end

Gamma_tilde = params.num_blades .* Gamma_updated
H_tilde = params.Omega * params.num_blades * Gamma_updated / (2.0 * pi)

# - get average meridional velocity at rotor - #
Vmavg = [
     ( Wm_rotor[1])
    0.5 .* (Wm_rotor[1:(end - 1)] .+ Wm_rotor[2:end])
     (Wm_rotor[end] )
]

# - assemble components of wake strength expression - #
dGammatilde2 = [
    Gamma_tilde[1]^2
    Gamma_tilde[2:end] .^ 2 .- Gamma_tilde[1:(end - 1)] .^ 2
    -Gamma_tilde[end]^2
]
dHtilde = [H_tilde[1]; H_tilde[2:end] .- H_tilde[1:(end - 1)]; -H_tilde[end]]

# - update wake strengths - #
gamma_theta_updated =
    1.0 ./ (2.0 .* Vmavg) .*
    (-(1.0 ./ (2.0 .* pi .* params.rotor_panel_edges)).^2 .* (dGammatilde2) .+ 2.0 .* (dHtilde))

##### ----- NEW FUNCTIONS ----- #####
function extract_state_variables(states, params)

    # Problem Dimensions
    nrotor = params.num_rotors                             # number of rotors
    nr = length(params.blade_elements[1].radial_positions) # number of radial coordinates

    # State Variable Indices
    iΓr = 1:(nr * nrotor) # rotor circulation strength indices
    iΓw = (iΓr[end] + 1):(iΓr[end] + nrotor * (nr + 1))     # wake vortex strength indices

    # Extract State variables
    Γr = reshape(view(states, iΓr), (nr, nrotor)) # rotor circulation strengths
    Γw = reshape(view(states, iΓw), (nr + 1, nrotor)) # wake circulation strengths

    return Γr, Γw
end

##### ----- NEW STUFF ----- #####
Gamma_new, gamma_theta_new = extract_state_variables(states, params)

wake_vortex_strengths = repeat(gamma_theta_new; inner=(1, params.nxwake))

rpc = params.rotor_panel_centers

TF = eltype(Gamma_new)

vx_rotor_new, vr_rotor_new, vtheta_rotor_new = dt.calculate_induced_velocities_on_rotors(
    params.blade_elements, Gamma_new, params.vx_rw, params.vr_rw, wake_vortex_strengths
)

# the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
Wx_rotor_new = vx_rotor_new .+ params.Vinf
# the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
Wtheta_rotor_new = vtheta_rotor_new .- params.Omega .* rpc

Wm_rotor_new = sqrt.(Wx_rotor_new .^ 2 .+ vr_rotor_new .^ 2)

# - Get the inflow magnitude at the rotor as the combination of all the components - #
Wmag_rotor_new = sqrt.(Wx_rotor_new .^ 2 .+ vr_rotor_new .^ 2 .+ Wtheta_rotor_new .^ 2)

# - Update Gamma_new - #
# TODO: check that this is outputting as expected. make sure gammas and sigmas aren't flipped.
# TODO: need to add/update unit test
dt.calculate_gamma_sigma!(
    Gamma_new,
    similar(states[params.Gammaidx]).=0,
    params.blade_elements,
    Wm_rotor_new,
    Wtheta_rotor_new,
)

#TODO: need to unit test these functions
Gamma_tilde_new = dt.calculate_net_circulation(Gamma_new, params.num_blades)
H_tilde_new = dt.calculate_enthalpy_jumps(Gamma_new, params.Omega, params.num_blades)

# - update wake strengths - #
# TODO: need to unit test this function
dt.calculate_wake_vortex_strengths!(
    gamma_theta_new, params.rotor_panel_edges, Wm_rotor_new, Gamma_tilde_new, H_tilde_new
)
