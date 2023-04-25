# - Extract commonly used items from precomputed inputs - #
blade_elements = inputs.blade_elements
rpc = inputs.rotor_panel_centers
Vinf = inputs.Vinf

# - Fill out wake strengths - #
wake_vortex_strengths = dt.fill_out_wake_strengths(
    gamw, inputs.rotor_indices, inputs.num_wake_x_panels
)

# - Get the induced velocities at the rotor plane - #
vx_rotor, vr_rotor, vtheta_rotor = dt.calculate_induced_velocities_on_rotors(
    blade_elements,
    Gamr,
    inputs.vx_rw,
    inputs.vr_rw,
    wake_vortex_strengths,
    inputs.vx_rr,
    inputs.vr_rr,
    zeros(size(Gamr)),
    # sigr,
    # inputs.vx_rb,
    # inputs.vr_rb,
    # gamb,
)

# the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
Wx_rotor = vx_rotor .+ inputs.Vinf

# the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
Wtheta_rotor = vtheta_rotor .- inputs.blade_elements[1].Omega .* rpc

# meridional component
Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

# Get the inflow magnitude at the rotor as the combination of all the components
Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

#---------------------------------#
#         PLOT AFTER SETUP        #
#---------------------------------#
pvx = plot(vx_rotor, inputs.rotor_panel_centers; xlabel=L"v_x", ylabel="r", label="initial")
savefig("examples/debug_coupling/vx-in-analysis.pdf")

pvr = plot(vr_rotor, inputs.rotor_panel_centers; xlabel=L"v_r", ylabel="r", label="initial")
savefig("examples/debug_coupling/vr-in-analysis.pdf")

pvt = plot(
    vtheta_rotor,
    inputs.rotor_panel_centers;
    xlabel=L"v_\theta",
    ylabel="r",
    label="initial",
)
savefig("examples/debug_coupling/vtheta-in-analysis.pdf")

pwx = plot(Wx_rotor, inputs.rotor_panel_centers; xlabel=L"W_x", ylabel="r", label="initial")
savefig("examples/debug_coupling/wx-in-analysis.pdf")

pwt = plot(
    Wtheta_rotor,
    inputs.rotor_panel_centers;
    xlabel=L"W_\theta",
    ylabel="r",
    label="initial",
)
savefig("examples/debug_coupling/wtheta-in-analysis.pdf")

pwm = plot(Wm_rotor, inputs.rotor_panel_centers; xlabel=L"W_m", ylabel="r", label="initial")
savefig("examples/debug_coupling/wm-in-analysis.pdf")

#---------------------------------#
#         continue analysis       #
#---------------------------------#

# - Calculate body vortex strengths - #
dt.calculate_body_vortex_strengths!(
    gamb,
    inputs.A_bb,
    inputs.b_bf,
    inputs.kutta_idxs,
    inputs.A_bw,
    wake_vortex_strengths,
    inputs.A_br,
    sigr,
)

# dt.calculate_gamma_sigma!(
# Gamr, sigr, inputs.blade_elements, Wm_rotor, Wtheta_rotor, Wmag_rotor
# )

# problem dimensions
nr, nrotor = size(Gamr) #num radial stations, num rotors

### --- Set up for Plotting --- ###
alpha = zeros(nr, nrotor)
phi = zeros(nr, nrotor)
cl = zeros(nr, nrotor)
cd = zeros(nr, nrotor)

# loop through rotors
for irotor in 1:nrotor

    # loop through radial stations
    for ir in 1:nr

        # extract blade element properties
        B = blade_elements[irotor].B # number of blades
        c = blade_elements[irotor].chords[ir] # chord length
        twist = blade_elements[irotor].twists[ir] # twist
        r = blade_elements[irotor].rbe[ir] # radius
        Ω = blade_elements[irotor].Omega # rotation rate

        # calculate angle of attack
        phi[ir, irotor] = atan(Wm_rotor[ir, irotor], -Wtheta_rotor[ir, irotor])
        alpha[ir, irotor] = twist - phi[ir, irotor]

        # look up lift and drag data for the nearest two input sections
        clin, cdin = dt.search_polars(
            blade_elements[irotor].inner_airfoil[ir], alpha[ir, irotor]
        )
        clout, cdout = dt.search_polars(
            blade_elements[irotor].outer_airfoil[ir], alpha[ir, irotor]
        )
        # linearly interpolate between those two values at your blade element location
        cl[ir, irotor] = fm.linear(
            [0.0; 1.0], [clin, clout], blade_elements[irotor].inner_fraction[ir]
        )
        cd[ir, irotor] = fm.linear(
            [0.0; 1.0], [cdin, cdout], blade_elements[irotor].inner_fraction[ir]
        )

        # calculate vortex strength
        Gamr[ir, irotor] = 1 / 2 * Wmag_rotor[ir, irotor] * c * cl[ir, irotor]

        # calculate source strength
        sigr[ir, irotor] = B / (4 * pi * r) * Wmag_rotor[ir, irotor] * c * cd[ir, irotor]
    end
end

#---------------------------------#
# PLOT BLADE ELEMENT PROPERTIES   #
#---------------------------------#

pa = plot(
    alpha[:, 1] * 180.0 / pi, inputs.rotor_panel_centers; xlabel=L"\alpha~(deg)", ylabel="r"
)
savefig(pa, "examples/debug_coupling/angles-of-attack-in-analysis.pdf")
pp = plot(
    phi[:, 1] * 180.0 / pi, inputs.rotor_panel_centers; xlabel=L"\phi~(deg)", ylabel="r"
)
savefig(pp, "examples/debug_coupling/inflow-angles-in-analysis.pdf")
pcl = plot(cl[:, 1], inputs.rotor_panel_centers; xlabel=L"c_\ell", ylabel="r")
savefig(pcl, "examples/debug_coupling/lift-coeffs-in-analysis.pdf")
pcd = plot(cd[:, 1], inputs.rotor_panel_centers; xlabel=L"c_d", ylabel="r")
savefig(pcd, "examples/debug_coupling/drag-coeffs-in-analysis.pdf")

#---------------------------------#
#         continue analysis       #
#---------------------------------#

# - Calculate net circulation and enthalpy jumps - #
Gamma_tilde = dt.calculate_net_circulation(Gamr, blade_elements.B)
H_tilde = dt.calculate_enthalpy_jumps(Gamr, blade_elements.Omega, blade_elements.B)

# - update wake strengths - #
dt.calculate_wake_vortex_strengths!(
    gamw, inputs.rotor_panel_edges, Wm_rotor, Gamma_tilde, H_tilde
)

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
# gamd = gamb[1:length(dp)]
# gamh = gamb[(length(dp) + 1):end]
gamd = 1.0 .- (gamb[1:length(dp)] ./ Vinf) .^ 2
# gamh = 1.0 .- (gamb[(length(dp) + 1):end]./Vinf).^2

plot!(pb, dp, gamd; xlabel="x", ylabel="cp", label="iter duct surface pressure")
# plot!(pb, hp, gamh ; label="iter hub surface pressure")

## -- check rotor circulation and source initial strengths -- ##
plot!(pG, Gamr, inputs.rotor_panel_centers; xlabel=L"\Gamma", ylabel="r", label="iter")

plot!(ps, sigr, inputs.rotor_panel_centers; xlabel=L"\sigma", ylabel="r", label="iter")

plot!(pw, gamw, inputs.rotor_panel_edges; xlabel=L"\gamma_\theta", ylabel="r", label="iter")
savefig(pb, "examples/debug_coupling/body-strengths-in-analysis.pdf")
savefig(pG, "examples/debug_coupling/circulation-in-analysis.pdf")
savefig(ps, "examples/debug_coupling/source-strengths-in-analysis.pdf")
savefig(pw, "examples/debug_coupling/wake-strengths-in-analysis.pdf")
