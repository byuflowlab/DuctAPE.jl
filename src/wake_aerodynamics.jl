function radially_average_velocity(Vmr, nx)
    Vmavg = [Vmr[1]; [0.5 * (Vmr[i] + Vmr[i + 1]) for i in 1:(length(Vmr) - 1)]; Vmr[end]]

    if nx > 1
        Vms = repeat(Vmavg; inner=(1, nx))
        return Vms
    else
        return Vmavg
    end
end

# function calculate_wake_on_wake_average_velocities(vx_ww, vr_ww, gamw)

#     # - Initialize - #
#     # get floating point type
#     TF = eltype(gamw)
#     # get dimensions
#     nx = length(gamw[1, :])
#     nr = length(gamw[:, 1]) - 1
#     vxw = zeros(TF, nr, nx) # axial induced velocity
#     vrw = zeros(TF, nr, nx) # radial induced velocity
#     vxbar = similar(gamw) .= 0.0
#     vrbar = similar(gamw) .= 0.0

#     # - Loop through affected wake "rotor" planes - #
#     for iplane in 1:nx

#         # - Loop through wake vortex sheets - #
#         # add wake induced velocities
#         for jwake in 1:nr
#             @views vxw[:, iplane] .+= vx_ww[iplane,jwake] * gamw[jwake, :]
#             @views vrw[:, iplane] .+= vr_ww[iplane,jwake] * gamw[jwake, :]
#         end

#     end

#     # - Average velocities - #
#     for iplane in 1:nx
#         vxbar[:, iplane] = radially_average_velocity(vxw[:, iplane], 1)
#         vrbar[:, iplane] = radially_average_velocity(vrw[:, iplane], 1)
#     end

#     return vxbar, vrbar
# end

"""
"""
function calculate_induced_velocities_on_wakes(
    vx_ww,
    vr_ww,
    gamw,
    vx_wr,
    vr_wr,
    sigr,
    vx_wb=nothing,
    vr_wb=nothing,
    gamb=nothing;
    debug=false,
)

    # get number of wakes
    nwake, _ = size(gamw)
    # get number of rotors
    _, nrotor = size(sigr)

    # initialize outputs
    vx = similar(gamw) .= 0 # axial induced velocity
    vr = similar(gamw) .= 0 # radial induced velocity

    if debug
        # initialize outputs
        vxb = similar(gamw) .= 0 # axial induced velocity
        vrb = similar(gamw) .= 0 # radial induced velocity
        vxr = similar(gamw) .= 0 # axial induced velocity
        vrr = similar(gamw) .= 0 # radial induced velocity
        vxw = similar(gamw) .= 0 # axial induced velocity
        vrw = similar(gamw) .= 0 # radial induced velocity
    end

    # loop through each wake sheet
    for iwake in 1:nwake

        # add body induced velocities
        if gamb != nothing
            @views vx[iwake, :] .+= vx_wb[iwake] * gamb
            @views vr[iwake, :] .+= vr_wb[iwake] * gamb
            if debug
                @views vxb[iwake, :] .+= vx_wb[irotor] * gamb
                @views vrb[iwake, :] .+= vr_wb[irotor] * gamb
            end
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vx[iwake, :] .+= vx_wr[iwake, jrotor] * sigr[:, jrotor]
            @views vr[iwake, :] .+= vr_wr[iwake, jrotor] * sigr[:, jrotor]
            if debug
                @views vxr[iwake, :] .+= vx_wr[iwake, jrotor] * sigr[:, jrotor]
                @views vrr[iwake, :] .+= vr_wr[iwake, jrotor] * sigr[:, jrotor]
            end
        end

        # add wake induced velocities
        for jwake in 1:nwake
            @views vx[iwake, :] .+= vx_ww[iwake, jwake] * gamw[jwake, :]
            @views vr[iwake, :] .+= vr_ww[iwake, jwake] * gamw[jwake, :]
            if debug
                @views vxw[iwake, :] .+= vx_ww[iwake, jwake] * gamw[jwake, :]
                @views vrw[iwake, :] .+= vr_ww[iwake, jwake] * gamw[jwake, :]
            end
        end
    end

    # return raw induced velocities
    if debug
        return vx, vr, vxb, vrb, vxr, vrr, vxw, vrw
    else
        return vx, vr
    end
end

function reframe_wake_velocities(vx_wake, vr_wake, Vinf)
    #add freestream to induced axial velocity
    Vx = vx_wake .+ Vinf

    # return meridional velocities
    return sqrt.(Vx .^ 2 .+ vr_wake .^ 2)
end

function get_sheet_jumps(Gamma_tilde, H_tilde)
    # get problem size
    nr, nrotor = size(Gamma_tilde)
    # number of wakes is one more than number or blade elements
    nw = nr + 1

    # get floating point type
    TF = eltype(Gamma_tilde)

    deltaGamma2 = zeros(TF, nw, nrotor)
    deltaH = zeros(TF, nw, nrotor)

    for irotor in 1:nrotor

        # loop through each radial position
        for iw in 1:nw
            # net circulation squared and enthalpy jumps across sheet
            if iw == 1
                #minus zero at hub
                deltaGamma2[iw, irotor] = Gamma_tilde[iw, irotor]^2
                # deltaH[iw, irotor] = H_tilde[iw, irotor]
                #NOTE: DFDC sets this to zero and comments out their equivalent to the line above this note.  It doesn't seem to affect much in the solve though.
                deltaH[iw, irotor] = 0.0
            elseif iw == nw
                # zero minus at duct
                deltaGamma2[iw, irotor] = -Gamma_tilde[iw - 1, irotor]^2
                deltaH[iw, irotor] = -H_tilde[iw - 1, irotor]
            else
                # difference between
                deltaGamma2[iw, irotor] =
                    Gamma_tilde[iw, irotor]^2 - Gamma_tilde[iw - 1, irotor]^2
                deltaH[iw, irotor] = H_tilde[iw, irotor] - H_tilde[iw - 1, irotor]
            end
        end
    end

    return deltaGamma2, deltaH
end

function get_wake_k(wake_vortex_panels)

    # get floating point type
    TF = eltype(wake_vortex_panels[1].panel_center)

    # initialize
    K = zeros(
        TF, length(wake_vortex_panels), length(wake_vortex_panels[1].panel_center[:, 1])
    )
    nw, np = size(K)

    for iw in 1:nw
        for ip in 1:np
            K[iw, ip] = -1.0 / (8.0 * pi^2 * wake_vortex_panels[iw].panel_center[ip, 2]^2)
        end
    end

    return K
end

"""

Calculate wake vortex strengths

"""
function calculate_wake_vortex_strengths!(Gamr, gamw, sigr, gamb, inputs; debug=false)

    # - get problem sizes - #
    # number of streamwise and radial panels
    nw, nx = size(gamw)

    # NOTE: this now includes wake
    # get induced velocities from body and rotor on wake panels.
    vx_wake, vr_wake = calculate_induced_velocities_on_wakes(
        inputs.vx_ww,
        inputs.vr_ww,
        gamw,
        inputs.vx_wr,
        inputs.vr_wr,
        sigr,
        inputs.vx_wb,
        inputs.vr_wb,
        gamb;
        debug=false,
    )

    ## TODO: remove this patch in favor of actual wake-on-wake contributions
    ## TODO: begin stuff to be removed
    ########################################
    ##### patch: use wake induced velocity on first rotor plane
    # # get induced velocities at rotor plane
    #vxr, vrr, _, vxb, vrb, vxw, vrw, _, _ = calculate_induced_velocities_on_rotors(
    #    inputs.blade_elements,
    #    Gamr,
    #    inputs.vx_rw,
    #    inputs.vr_rw,
    #    gamw,
    #    inputs.vx_rr,
    #    inputs.vr_rr,
    #    sigr,
    #    inputs.vx_rb,
    #    inputs.vr_rb,
    #    gamb;
    #    debug=true,
    #)

    ## # - average rotor plane velocities
    ## # TODO: this is part of the patch that needs to be removed
    #vxavg = radially_average_velocity(vxw, nx)
    #vravg = radially_average_velocity(vrw, nx)

    ###TODO: also part that needs to be removed
    ### - add rotor plane velocities to wake (this is to replace the wake-on-wake interactions that you don't have. it won't be right, but hopefully it won't be terribly wrong)
    #vx_wake .+= vxavg
    #vr_wake .+= vravg

    # TODO: end stuff to be removed
    #######################################

    # reframe velocities to get meridional velocity on wake panels
    Vm = reframe_wake_velocities(vx_wake, vr_wake, inputs.Vinf)

    # get net circulation of upstream rotors
    Gamma_tilde = calculate_net_circulation(Gamr, inputs.blade_elements.B)

    # get enthalpy jump across disks
    H_tilde = calculate_enthalpy_jumps(
        Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
    )

    # get the circulation squared and enthalpy jumps across the wake sheets
    deltaGamma2, deltaH = get_sheet_jumps(Gamma_tilde, H_tilde)

    # calculate radius dependent "constant" for wake strength calcualtion
    K = get_wake_k(inputs.wake_vortex_panels)

    #loop  through rotors where wake strengths are generated
    for iw in 1:nw

        # loop through each radial position
        for ix in 1:nx
            #get index of last rotor before current x location
            irotor = searchsortedlast(inputs.rotor_indices_in_wake, ix)

            # calculate the wake vortex strength
            if Vm[iw, ix] <= 0.0
                # avoid division by zero
                gamw[iw, ix] = 0.0
            else

                # wake strength density taken from rotor to next rotor constant along streamlines
                gamw[iw, ix] =
                    (K[iw, ix] * deltaGamma2[iw, irotor] + deltaH[iw, irotor]) / Vm[iw, ix]
            end
        end
    end

    # - Wake-Body Interface Treatment - #
    if inputs.ductTE_index != nothing
        gamw[end, 1:(inputs.ductTE_index)] .= 0.0

        # gamw[end, 1:(inputs.ductTE_index)] = range(
        #     0.0, gamw[end, inputs.ductTE_index]; length=inputs.ductTE_index
        # )
    end
    if inputs.hubTE_index != nothing
        gamw[1, 1:(inputs.hubTE_index)] .= 0.0

        # gamw[1, 1:(inputs.hubTE_index)] = range(
        #     0.0, gamw[1, inputs.hubTE_index]; length=inputs.hubTE_index
        # )
    end

    if debug
        return gamw, deltaGamma2, deltaH, K
    else
        return gamw
    end
end

function initialize_wake_vortex_strengths(Vinf, Gamr, Omega, B, rotor_panel_edges)

    # get net circulations
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # get enthalpy jumps
    H_tilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # initialize vm with freestream and rotation assuming open rotor
    vm = marched_vm(Vinf, Gamma_tilde, H_tilde, rotor_panel_edges)

    # return initialized wake vortex strengths
    gwhub = vm[1, :] .- Vinf
    gw = vm[2:end, :] .- vm[1:(end - 1), :]
    gwduct = Vinf .- vm[end, :]
    return [gwhub'; gw; gwduct']
end

"""
get meridional velocities using the marching method for when Vinf outside of wake is known.
"""
function marched_vm(Vinf, Gamma_tilde, H_tilde, rotor_panel_edges)

    #get floating point type
    TF = eltype(Gamma_tilde)

    #rename for convenience
    nr, nrotor = size(Gamma_tilde)

    #initialize
    Vm = zeros(TF, nr, nrotor)

    # Loop through rotors
    for irotor in 1:nrotor

        # Loop through radial stations
        #march from just outside the tip to the hub
        for ir in nr:-1:1
            if ir == nr

                #this is the change at the tip, so we need to set Vm2 to Vinf, and set Gamma dn H 2's to zeros
                radical =
                    Vinf^2 +
                    (1.0 / (2.0 * pi * rotor_panel_edges[ir, irotor]))^2 *
                    (-Gamma_tilde[ir, irotor]^2) - 2.0 * (-H_tilde[ir, irotor])
            else

                #otherwise, we just take the differences inside as-is
                radical =
                    Vm[ir + 1, irotor]^2 +
                    (1.0 / (2.0 * pi * rotor_panel_edges[ir, irotor]))^2 *
                    (Gamma_tilde[ir + 1, irotor]^2 - Gamma_tilde[ir, irotor]^2) -
                    2.0 * (H_tilde[ir + 1, irotor] - H_tilde[ir, irotor])
            end

            # compute the meridional velocity value, avoiding negative square roots and divisions by zero later
            if radical > 0.0
                Vm[ir, irotor] = sqrt(radical)
            else
                Vm[ir, irotor] = 0.1 * Vinf
            end
        end
    end

    #Vm should now be in the order of hub to tip
    return Vm
end
