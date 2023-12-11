function radially_average_velocity(Vmr, nx)
    Vmavg = [Vmr[1]; [0.5 * (Vmr[i] + Vmr[i + 1]) for i in 1:(length(Vmr) - 1)]; Vmr[end]]

    if nx > 1
        Vms = repeat(Vmavg; inner=(1, nx))
        return Vms
    else
        return Vmavg
    end
end

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
    gamb=nothing,
    vx_wbte=nothing,
    vr_wbte=nothing,
    TEidxs=nothing;
    debug=false,
)

    # get number of rotors
    nrotor = size(sigr, 2)

    # initialize outputs
    vx_wake = similar(gamw, size(vx_ww, 1), 1) .= 0 # axial induced velocity
    vr_wake = similar(gamw, size(vr_ww, 1), 1) .= 0 # axial induced velocity

    if debug
        # initialize outputs
        vxb_wake = similar(vx_wake) .= 0 # axial induced velocity
        vrb_wake = similar(vx_wake) .= 0 # radial induced velocity
        vxr_wake = similar(vx_wake) .= 0 # axial induced velocity
        vrr_wake = similar(vx_wake) .= 0 # radial induced velocity
        vxw_wake = similar(vx_wake) .= 0 # axial induced velocity
        vrw_wake = similar(vx_wake) .= 0 # radial induced velocity
    end

    # add body induced velocities
    if gamb != nothing
        @views vx_wake .+= vx_wb * gamb
        @views vr_wake .+= vr_wb * gamb
        if debug
            @views vxb_wake .+= vx_wb * gamb
            @views vrb_wake .+= vr_wb * gamb
        end

        @views vx_wake .+= vx_wbte * gamb[TEidxs]
        @views vr_wake .+= vr_wbte * gamb[TEidxs]
        if debug
            @views vxb_wake .+= vx_wbte * gamb[TEidxs]
            @views vrb_wake .+= vr_wbte * gamb[TEidxs]
        end
    end

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vx_wake .+= vx_wr[jrotor] * sigr[:, jrotor]
        @views vr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
        if debug
            @views vxr_wake .+= vx_wr[jrotor] * sigr[:, jrotor]
            @views vrr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
        end
    end

    # add wake induced velocities
    @views vx_wake .+= vx_ww * gamw
    @views vr_wake .+= vr_ww * gamw
    if debug
        @views vxw_wake .+= vx_ww * gamw
        @views vrw_wake .+= vr_ww * gamw
    end

    # return raw induced velocities
    if debug
        return vx_wake, vr_wake, vxb_wake, vrb_wake, vxr_wake, vrr_wake, vxw_wake, vrw_wake
    else
        return vx_wake, vr_wake
    end
end

function reframe_wake_velocities(vx_wake, vr_wake, Vinf)
    #add freestream to induced axial velocity
    Wx_wake = vx_wake .+ Vinf

    # return meridional velocities
    return sqrt.(Wx_wake .^ 2 .+ vr_wake .^ 2)
end

"""
"""
function calculate_wake_velocities(gamw, sigr, gamb, inputs)

    # - Get induced velocities on wake - #
    vx_wake, vr_wake = calculate_induced_velocities_on_wakes(
        inputs.vx_ww,
        inputs.vr_ww,
        gamw,
        inputs.vx_wr,
        inputs.vr_wr,
        sigr,
        inputs.vx_wb,
        inputs.vr_wb,
        gamb,
        inputs.vx_wbte,
        inputs.vr_wbte,
        (p -> p.idx).(inputs.body_vortex_panels.TEnodes),
    )

    # - Reframe rotor velocities into blade element frames
    return reframe_wake_velocities(vx_wake, vr_wake, inputs.freestream.Vinf)
end

function get_sheet_jumps(Gamma_tilde, H_tilde)
    # get problem size
    nr, nrotor = size(Gamma_tilde)
    # number of wakes is one more than number of blade elements
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

"""

Calculate wake vortex strengths

"""
function calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs; debug=false)

    # get net circulation of upstream rotors
    Gamma_tilde = calculate_net_circulation(Gamr, inputs.blade_elements.B)

    # get enthalpy jump across disks
    H_tilde = calculate_enthalpy_jumps(
        Gamr, inputs.blade_elements.Omega, inputs.blade_elements.B
    )

    # get the circulation squared and enthalpy jumps across the wake sheets
    deltaGamma2, deltaH = get_sheet_jumps(Gamma_tilde, H_tilde)

    for (iw, (gw, K, sheetid, rotorid)) in enumerate(
        zip(eachrow(gamw), inputs.wakeK, inputs.rotorwakeid[:, 1], inputs.rotorwakeid[:, 2])
    )

        # Get average meridional velocity at node
        if iw in inputs.wake_vortex_panels.endnodeidxs[1, :]
            Wm_avg = Wm_wake[searchsortedfirst(inputs.wake_vortex_panels.nodemap[1, :], iw)]
        elseif iw in inputs.wake_vortex_panels.endnodeidxs[2, :]
            Wm_avg = Wm_wake[searchsortedfirst(inputs.wake_vortex_panels.nodemap[2, :], iw)]
        else
            Wm_avg =
                0.5 * (
                    Wm_wake[searchsortedfirst(
                        inputs.wake_vortex_panels.nodemap[1, :], iw
                    )] +
                    Wm_wake[searchsortedfirst(inputs.wake_vortex_panels.nodemap[2, :], iw)]
                )
        end

        # calculate the wake vortex strength
        if Wm_avg <= 0.0
            # avoid division by zero
            gw[1] = 0.0
        else

            # wake strength density taken from rotor to next rotor constant along streamlines
            gw[1] = (K * deltaGamma2[sheetid, rotorid] + deltaH[sheetid, rotorid]) / Wm_avg
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

    if debug
        return gamw, deltaGamma2, deltaH
    else
        return gamw
    end
end

function initialize_wake_vortex_strengths(Vinf, Gamr, Omega, B, rotor_panel_edges, nxwake)

    # get net circulations
    Gamma_tilde = calculate_net_circulation(Gamr, B)

    # get enthalpy jumps
    H_tilde = calculate_enthalpy_jumps(Gamr, Omega, B)

    # initialize Wm_wake with freestream and rotation assuming open rotor
    Wm_wake = marched_vm(Vinf, Gamma_tilde, H_tilde, rotor_panel_edges)

    # return initialized wake vortex strengths
    gwhub = Wm_wake[1, :] .- Vinf
    gw = Wm_wake[2:end, :] .- Wm_wake[1:(end - 1), :]
    gwduct = Vinf .- Wm_wake[end, :]
    planegw = [gwhub'; gw; gwduct']
    gamw = repeat(planegw; inner=(1, nxwake))
    return reduce(vcat, gamw')
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
    Wm_wake = zeros(TF, nr, nrotor)

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
                    Wm_wake[ir + 1, irotor]^2 +
                    (1.0 / (2.0 * pi * rotor_panel_edges[ir, irotor]))^2 *
                    (Gamma_tilde[ir + 1, irotor]^2 - Gamma_tilde[ir, irotor]^2) -
                    2.0 * (H_tilde[ir + 1, irotor] - H_tilde[ir, irotor])
            end

            # compute the meridional velocity value, avoiding negative square roots and divisions by zero later
            if radical > 0.0
                Wm_wake[ir, irotor] = sqrt(radical)
            else
                Wm_wake[ir, irotor] = 0.1 * Vinf
            end
        end
    end

    #Wm_wake should now be in the order of hub to tip
    return Wm_wake
end
