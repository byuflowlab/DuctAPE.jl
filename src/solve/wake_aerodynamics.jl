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
    vz_ww,
    vr_ww,
    gamw,
    vz_wr,
    vr_wr,
    sigr,
    vz_wb=nothing,
    vr_wb=nothing,
    gamb=nothing;
    post=false,
)

    # get number of rotors
    nrotor = size(sigr, 2)

    # initialize outputs
    vz_wake = similar(gamw, size(vz_ww, 1), 1) .= 0 # axial induced velocity
    vr_wake = similar(gamw, size(vr_ww, 1), 1) .= 0 # axial induced velocity

    if post
        # initialize outputs
        vzb_wake = similar(vz_wake) .= 0 # axial induced velocity
        vrb_wake = similar(vz_wake) .= 0 # radial induced velocity
        vzr_wake = similar(vz_wake) .= 0 # axial induced velocity
        vrr_wake = similar(vz_wake) .= 0 # radial induced velocity
        vzw_wake = similar(vz_wake) .= 0 # axial induced velocity
        vrw_wake = similar(vz_wake) .= 0 # radial induced velocity
    end

    # add body induced velocities
    if gamb != nothing
        #TODO: look into signs here, may need to flip them
        @views vz_wake .+= vz_wb * gamb
        @views vr_wake .+= vr_wb * gamb

        if post
            @views vzb_wake .+= vz_wb * gamb
            @views vrb_wake .+= vr_wb * gamb
        end
    end

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vz_wake .+= vz_wr[jrotor] * sigr[:, jrotor]
        @views vr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
        if post
            @views vzr_wake .+= vz_wr[jrotor] * sigr[:, jrotor]
            @views vrr_wake .+= vr_wr[jrotor] * sigr[:, jrotor]
        end
    end

    # add wake induced velocities
    @views vz_wake .+= vz_ww * gamw
    @views vr_wake .+= vr_ww * gamw

    #  return raw induced velocities
    if post
        return vz_wake, vr_wake, vzb_wake, vrb_wake, vzr_wake, vrr_wake, vzw_wake, vrw_wake
    else
        return vz_wake, vr_wake
    end
end

function reframe_wake_velocities(vz_wake, vr_wake, Vinf; post=false)
    #add freestream to induced axial velocity
    Wz_wake = vz_wake .+ Vinf

    if post
        return Wz_wake, sqrt.(Wz_wake .^ 2 .+ vr_wake .^ 2)
    else
        # return meridional velocities
        return sqrt.(Wz_wake .^ 2 .+ vr_wake .^ 2)
    end
end

"""
"""
function calculate_wake_velocities(gamw, sigr, inputs)

    # - Get induced velocities on wake - #
    vz_wake, vr_wake = calculate_induced_velocities_on_wakes(
        inputs.vz_ww, inputs.vr_ww, gamw, inputs.vz_wr, inputs.vr_wr, sigr
    )

    # - Reframe rotor velocities into blade element frames
    return reframe_wake_velocities(vz_wake, vr_wake, inputs.freestream.Vinf)
end

"""
"""
function calculate_wake_velocities(gamw, sigr, gamb, inputs)

    # - Get induced velocities on wake - #
    vz_wake, vr_wake = calculate_induced_velocities_on_wakes(
        inputs.vz_ww,
        inputs.vr_ww,
        gamw,
        inputs.vz_wr,
        inputs.vr_wr,
        sigr,
        inputs.vz_wb,
        inputs.vr_wb,
        gamb,
    )

    # - Reframe rotor velocities into blade element frames
    return reframe_wake_velocities(vz_wake, vr_wake, inputs.freestream.Vinf)
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
function calculate_wake_vortex_strengths!(gamw, Gamr, Wm_wake, inputs; post=false)

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
        if Wm_avg < eps()
            # avoid division by zero
            gw[1] = 0.0
        else

            # wake strength density taken from rotor to next rotor constant along streamlines
            gw[1] = -(K * deltaGamma2[sheetid, rotorid] + deltaH[sheetid, rotorid]) / Wm_avg
        end
    end

    # - Wake-Body Interface Treatment - #
    if inputs.ductwakeinterfaceid != nothing
        # gamw[inputs.ductwakeinterfaceid] .= 0.0

        #TODO: double check this is working properly.
        gamw[inputs.ductwakeinterfacenodeid] =
            gamw[inputs.ductwakeinterfacenodeid[end]+1] * (
                1.0 .-
                (inputs.ductwakeinterfacenodeid[end] .- inputs.ductwakeinterfacenodeid) / (
                    inputs.ductwakeinterfacenodeid[end] -
                    inputs.ductwakeinterfacenodeid[1] + 1
                )
            )
    end
    if inputs.hubwakeinterfaceid != nothing
        # gamw[inputs.hubwakeinterfaceid] .= 0.0

        #TODO: double check this is working properly.
        println(gamw[inputs.hubwakeinterfacenodeid[end]+1])
        gamw[inputs.hubwakeinterfacenodeid] =
            gamw[inputs.hubwakeinterfacenodeid[end]+1] * (
                1.0 .-
                (reverse(inputs.hubwakeinterfacenodeid) .- inputs.hubwakeinterfacenodeid[1]) /
                length(inputs.hubwakeinterfacenodeid)
            )
    end

    if post
        return gamw, deltaGamma2, deltaH
    else
        return gamw
    end
end
