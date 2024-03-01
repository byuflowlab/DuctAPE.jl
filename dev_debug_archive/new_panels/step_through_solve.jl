project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/dev_debug_archive/dfdc_archive/dfdc_testing/"

include(datapath * "straightwake.case.jl")

include(project_dir * "/visualize/plots_default.jl")

using DuctAPE
const dt = DuctAPE

using FLOWMath
const fm = FLOWMath

# initialize various inputs used in analysis
inputs = dt.precomputed_inputs(
    duct_coordinates,
    hub_coordinates,
    paneling_constants,
    rotor_parameters,
    freestream,
    reference_parameters;
)

(;
    converged,
    #freestream
    freestream,
    #reference for post process
    reference_parameters,
    # - Panels - #
    body_doublet_panels,
    rotor_source_panels,
    wake_vortex_panels,
    wakeK, #constant based on radial location that is used to define the wake panel strength
    duct_wake_panels,
    hub_wake_panels,
    # - rotors - #
    blade_elements, # blade elements
    num_rotors,
    rotor_panel_edges,
    rotor_panel_centers,
    # - Book Keeping - #
    num_wake_x_panels, # number of wake panels in the axial direction
    num_body_panels,
    # ductTE_index=tip_gaps[1] == 0.0 ? ductTE_index : nothing,
    # hubTE_index=!nohub ? hubTE_index : nothing,
    ductidsaftofrotors,
    hubidsaftofrotors,
    ductwakeinterfaceid, # wake panel indices that lie on top of duct wall
    hubwakeinterfaceid, # wake panel indices that lie on top of hub wall
    rotorwakenodeid, # [rotor panel edge index, and closest forward rotor id] for each wake panel
    # - Linear System - #
    body_system_matrices, # includes the various LHS and RHS matrics and vectors for solving the linear system for the body
    # - Influence Matrices - #
    A_bb, # body to body
    b_bf, # freestream contribution to body boundary conditions
    prescribedpanels, # prescribed panels
    A_br, # rotor to body (total)
    A_bw, # wake to body (total)
    v_rb, # body to rotor
    v_rbte,
    v_wb, # body to wake
    v_wbte,
    v_wr, # rotor to wake
    v_ww, # wake to wake
    v_dwb, # body to duct wake
    v_dwr, # rotor to duct wake
    v_dww, # wake to duct wake
    v_hwb, # body to hub wake
    v_hwr, # rotor to hub wake
    v_hww, # wake to hub wake
    vz_rb, # body to rotor (x-direction)
    vr_rb, # body to rotor (r-direction)
    vz_rbte, # bodyTE to rotor (x-direction)
    vr_rbte, # bodyTE to rotor (r-direction)
    vz_rr, # rotor to rotor (x-direction)
    vr_rr, # rotor to rotor ( r-direction)
    vz_rw, # wake to rotor (x-direction)
    vr_rw, # wake to rotor ( r-direction)
    vz_wb, # body to wake (x-direction)
    vr_wb, # body to wake ( r-direction)
    vz_wbte, # bodyTE to wake (x-direction)
    vr_wbte, # bodyTE to wake ( r-direction)
    vz_wr, # rotor to wake (x-direction)
    vr_wr, # rotor to wake ( r-direction)
    vz_ww, # wake to wake (x-direction)
    vr_ww, # wake to wake ( r-direction)
    vz_dwb, # body to duct wake (x-direction)
    vr_dwb, # body to duct wake ( r-direction)
    vz_dwbte, # bodyTE to duct wake (x-direction)
    vr_dwbte, # bodyTE to duct wake ( r-direction)
    vz_dwr, # rotor to duct wake (x-direction)
    vr_dwr, # rotor to duct wake ( r-direction)
    vz_dww, # wake to duct wake (x-direction)
    vr_dww, # wake to duct wake ( r-direction)
    vz_hwb, # body to hub wake (x-direction)
    vr_hwb, # body to hub wake ( r-direction)
    vz_hwbte, # bodyTE to hub wake (x-direction)
    vr_hwbte, # bodyTE to hub wake ( r-direction)
    vz_hwr, # rotor to hub wake (x-direction)
    vr_hwr, # rotor to hub wake ( r-direction)
    vz_hww, # wake to hub wake (x-direction)
    vr_hww, # wake to hub wake ( r-direction)
    # operating conditions
    Vinf, # freestream parameters
    # - Debugging/Plotting
    duct_coordinates,
    hub_coordinates,
    isduct,
    ishub,
    wakexgrid,
    wakergrid,
) = inputs

states = dt.initialize_states(inputs)

# - Extract states - #
mub, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)

###############################################################################

debug = true
nrotor = num_rotors
TEidxs = (p -> p.idx).(body_doublet_panels.TEnodes)
### --- Get Velocities Before Updating States --- ###
vz_rotor = similar(Gamr) .= 0 # axial induced velocity
vr_rotor = similar(Gamr) .= 0 # radial induced velocity
vtheta_rotor = similar(Gamr) .= 0 # tangential induced velocity

if debug

    # initialize outputs
    vxb_rotor = similar(Gamr) .= 0 # axial induced velocity
    vrb_rotor = similar(Gamr) .= 0 # radial induced velocity

    vxw_rotor = similar(Gamr) .= 0 # axial induced velocity
    vrw_rotor = similar(Gamr) .= 0 # radial induced velocity

    vxr_rotor = similar(Gamr) .= 0 # axial induced velocity
    vrr_rotor = similar(Gamr) .= 0 # radial induced velocity
end

# loop through each affected rotor
for irotor in 1:nrotor

    # add body induced velocities
    if mub != nothing
        @views vz_rotor[:, irotor] .+= vz_rb[irotor] * mub
        @views vr_rotor[:, irotor] .+= vr_rb[irotor] * mub

        if debug
            @views vxb_rotor[:, irotor] .+= vz_rb[irotor] * mub
            @views vrb_rotor[:, irotor] .+= vr_rb[irotor] * mub
        end

        #note: v?_rbte[irotor] should have been defined as zeros for endpoints that are not coincident.
        @views vz_rotor .+= vz_rbte[irotor] * mub[TEidxs]
        @views vr_rotor .+= vr_rbte[irotor] * mub[TEidxs]

        if debug
            @views vxb_rotor .+= vz_rbte[irotor] * mub[TEidxs]
            @views vrb_rotor .+= vr_rbte[irotor] * mub[TEidxs]
        end
    end

    # add wake induced velocities
    @views vz_rotor[:, irotor] .+= vz_rw[irotor] * gamw
    @views vr_rotor[:, irotor] .+= vr_rw[irotor] * gamw

    if debug
        @views vxw_rotor[:, irotor] .+= vz_rw[irotor] * gamw[:]
        @views vrw_rotor[:, irotor] .+= vr_rw[irotor] * gamw[:]
    end

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vz_rotor[:, irotor] .+= vz_rr[irotor, jrotor] * sigr[:, jrotor]
        @views vr_rotor[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]

        if debug
            @views vxr_rotor[:, irotor] .+= vz_rr[irotor, jrotor] * sigr[:, jrotor]
            @views vrr_rotor[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]
        end
    end

    # add self-induced tangential velocity
    B = blade_elements[irotor].B
    r = blade_elements[irotor].rbe
    @views vtheta_rotor[:, irotor] .+= B .* Gamr[:, irotor] ./ (4 * pi * r)

    # add induced tangential velocity from upstream rotors
    for jrotor in 1:(irotor - 1)
        B = blade_elements[jrotor].B
        r = blade_elements[jrotor].rbe
        @views vtheta_rotor[:, irotor] .+= B .* Gamr[:, jrotor] ./ (2 * pi * r)
    end
end

# the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
Wx_rotor = vz_rotor .+ Vinf

# the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
Wtheta_rotor = similar(vtheta_rotor) .= vtheta_rotor

for i in 1:length(Omega)
    Wtheta_rotor[:, i] .-= Omega[i] .* rotor_panel_centers[:, i]
end

# meridional component
Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

# Get the inflow magnitude at the rotor as the combination of all the components
Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

#######################################################################################
# initialize outputs
vz_wake = similar(gamw) .= 0 # axial induced velocity
vr_wake = similar(gamw) .= 0 # radial induced velocity

if debug
    # initialize outputs
    vxb_wake = similar(gamw) .= 0 # axial induced velocity
    vrb_wake = similar(gamw) .= 0 # radial induced velocity
    vxr_wake = similar(gamw) .= 0 # axial induced velocity
    vrr_wake = similar(gamw) .= 0 # radial induced velocity
    vxw_wake = similar(gamw) .= 0 # axial induced velocity
    vrw_wake = similar(gamw) .= 0 # radial induced velocity
end

# add body induced velocities
if mub != nothing
    @views vz_wake .+= vz_wb * mub
    @views vr_wake .+= vr_wb * mub
    if debug
        @views vxb_wake .+= vz_wb * mub
        @views vrb_wake .+= vr_wb * mub
    end

    @views vz_wake .+= vz_wbte * mub[TEidxs]
    @views vr_wake .+= vr_wbte * mub[TEidxs]
    if debug
        @views vxb_wake .+= vz_wbte * mub[TEidxs]
        @views vrb_wake .+= vr_wbte * mub[TEidxs]
    end
end

# add rotor induced velocities
for jrotor in 1:nrotor
    @views vz_wake .+= vz_wr[jrotor] * sigr[:, jrotor]
    @views vr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
    if debug
        @views vxr_wake .+= vz_wr[jrotor] * sigr[:, jrotor]
        @views vrr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
    end
end

# add wake induced velocities
@views vz_wake .+= vz_ww * gamw
@views vr_wake .+= vr_ww * gamw
if debug
    @views vxw_wake .+= vz_ww * gamw
    @views vrw_wake .+= vr_ww * gamw
end

Vx_wake = vz_wake .+ Vinf

Wm_wake = sqrt.(Vx_wake .^ 2 .+ vr_wake .^ 2)

########################################################################################
LHS = A_bb

# - Generate raw RHS, viz. velocities on body, (before updating state dependencies) - #
RHS = zeros(size(b_bf, 1))
# start with freestream contributions to right hand side
# note: the negative was already included in the precomputation for the freestream.
RHS .+= b_bf

# add wake vortex sheet contributions to right hand side
# note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract here
RHS .-= A_bw * gamw

# add rotor source sheet contributions to right hand side
# note: the subtraction was not included in the other coefficients in the precomputation, so we need to subract here
for jrotor in 1:length(A_br)
    # get induced velocity in the x-direction
    RHS .-= A_br[jrotor] * view(sigr, :, jrotor)
end

LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

mured = LHSlsq \ RHSlsq

dt.mured2mu!(mub, mured, prescribedpanels)

########################################################################################

nw = length(gamw)

# get net circulation of upstream rotors
Gamma_tilde = dt.calculate_net_circulation(Gamr, inputs.blade_elements.B)

# get enthalpy jump across disks
H_tilde = dt.calculate_enthalpy_jumps(
    Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
)

# get the circulation squared and enthalpy jumps across the wake sheets
deltaGamma2, deltaH = dt.get_sheet_jumps(Gamma_tilde, H_tilde)

for iw in 1:nw

    # calculate the wake vortex strength
    if Wm_wake[iw] <= 0.0
        # avoid division by zero
        gamw[iw] = 0.0
    else

        # wake strength density taken from rotor to next rotor constant along streamlines
        gamw[iw] =
            (
                inputs.wakeK[iw] *
                deltaGamma2[inputs.rotorwakenodeid[iw, 1], inputs.rotorwakenodeid[iw, 2]] +
                deltaH[inputs.rotorwakenodeid[iw], inputs.rotorwakenodeid[iw, 2]]
            ) / Wm_wake[iw]
    end
end

# - Wake-Body Interface Treatment - #
if inputs.ductwakeinterfaceid != nothing
    gamw[inputs.ductwakeinterfaceid] .= 0.0

    # gamw[end, 1:(inputs.ductTE_index)] = range(
    #     0.0, gamw[end, inputs.ductTE_index]; length=inputs.ductTE_index
    # )
end
if inputs.hubwakeinterfaceid != nothing
    gamw[inputs.hubwakeinterfaceid] .= 0.0

    # gamw[1, 1:(inputs.hubTE_index)] = range(
    #     0.0, gamw[1, inputs.hubTE_index]; length=inputs.hubTE_index
    # )
end

#########################################################################################

# problem dimensions
nr, nrotor = size(Gamr) #num radial stations, num rotors

if debug
    # initialize extra outputs
    phidb = zeros(eltype(Gamr), size(Gamr))
    alphadb = zeros(eltype(Gamr), size(Gamr))
    cldb = zeros(eltype(Gamr), size(Gamr))
    cddb = zeros(eltype(Gamr), size(Gamr))
end

# loop through rotors
for irotor in 1:nrotor

    # loop through radial stations
    for ir in 1:nr

        # extract blade element properties
        B = blade_elements[irotor].B # number of blades
        c = blade_elements[irotor].chords[ir] # chord length
        twist = blade_elements[irotor].twists[ir] # twist
        #stagger is twist angle but from axis
        # stagger = 0.5 * pi - twist
        stagger = blade_elements[irotor].stagger[ir]
        r = blade_elements[irotor].rbe[ir] # radius
        Ω = blade_elements[irotor].Omega # rotation rate
        solidity = blade_elements[irotor].solidity[ir]

        # calculate angle of attack
        phi, alpha = dt.calculate_inflow_angles(
            Wm_rotor[ir, irotor], Wtheta_rotor[ir, irotor], twist
        )

        # printval("Wm_rotor[$(ir),$(irotor)]: ", Wm_rotor[ir,irotor])
        # printval("Wtheta[$(ir),$(irotor)]: ", Wtheta_rotor[ir,irotor])
        # printval("twist: ", twist*180.0/pi)
        # printval("phi: ", phi*180.0/pi)
        # printval("alpha: ", alpha*180.0/pi)
        # printval("Wmag: ", Wmag_rotor[ir, irotor])
        # printval("Mach: ", Wmag_rotor[ir,irotor] / freestream.asound)

        # look up lift and drag data for the nearest two input sections
        if typeof(blade_elements[irotor].inner_airfoil[ir]) <: dt.DFDCairfoil

            # get local reynolds number
            reynolds =
                c * abs(Wmag_rotor[ir, irotor]) * freestream.rhoinf / freestream.muinf

            # printval("Re: ", reynolds)

            #get inner values
            clin, cdin, _ = dt.dfdceval(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity,
                stagger,
                alpha,
                blade_elements[irotor].inner_airfoil[ir],
                freestream.asound;
                verbose=false,
            )
            # get outer values
            clout, cdout, _ = dt.dfdceval(
                Wmag_rotor[ir, irotor],
                reynolds,
                solidity,
                stagger,
                alpha,
                blade_elements[irotor].outer_airfoil[ir],
                freestream.asound;
                verbose=false,
            )

        else
            clin, cdin = dt.search_polars(blade_elements[irotor].inner_airfoil[ir], alpha)
            clout, cdout = dt.search_polars(blade_elements[irotor].outer_airfoil[ir], alpha)
            # linearly interpolate between those two values at your blade element location
        end

        # interpolate inner and outer values
        cl = fm.linear([0.0; 1.0], [clin, clout], blade_elements[irotor].inner_fraction[ir])
        cd = fm.linear([0.0; 1.0], [cdin, cdout], blade_elements[irotor].inner_fraction[ir])

        # printval("cl in: ", clin)
        # printval("cl interp: ", cl)

        dt.gamma_sigma_from_coeffs!(
            view(Gamr, ir, irotor),
            view(sigr, ir, irotor),
            Wmag_rotor[ir, irotor],
            B,
            c,
            r,
            cl,
            cd,
        )

        if debug
            phidb[ir, irotor] = phi
            alphadb[ir, irotor] = alpha
            cldb[ir, irotor] = cl
            cddb[ir, irotor] = cd
        end
    end
end
