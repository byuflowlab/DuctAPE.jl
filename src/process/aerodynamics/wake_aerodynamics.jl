"""
    calculate_induced_velocities_on_wakes(
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
    ) -> Tuple{Matrix{Float64}, Matrix{Float64}}

Compute the induced velocities on the wake panels due to contributions from wakes, rotors, and optionally bodies.

# Arguments
- `vz_ww::Matrix{Float}` : Axial influence of wakes on wakes.
- `vr_ww::Matrix{Float}` : Radial influence of wakes on wakes.
- `gamw::Vector{Float}` : Circulation strength for each wake panel.
- `vz_wr::Matrix{Float}` : Axial influence of rotor vortex sheets on wakes.
- `vr_wr::Matrix{Float}` : Radial influence of rotor vortex sheets on wakes.
- `sigr::Matrix{Float}` : Rotor panel source strengths.
- `vz_wb::Matrix{Float} (optional)` : Axial influence of body vortices on wakes. Defaults to `nothing`.
- `vr_wb::Matrix{Float} (optional)` : Radial influence of body vortices on wakes. Defaults to `nothing`.
- `gamb::Vector{Float} (optional)` : Body vortex strengths. Defaults to `nothing`.

# Keyword Arguments
- `post::Bool=false` : If true, additional diagnostic terms may be returned (from the in-place version).

# Returns
- `vz_wake::Matrix{Float}` : Total induced axial velocity on each wake panel.
- `vr_wake::Matrix{Float}` : Total induced radial velocity on each wake panel.

# Notes
- This function initializes output matrices and wraps the in-place computation performed by `calculate_induced_velocities_on_wakes!`.
- Dimensions must be consistent across input influence matrices and strength vectors.
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

    # initialize outputs
    vz_wake = similar(gamw, size(vz_ww, 1), 1) .= 0 # axial induced velocity
    vr_wake = similar(gamw, size(vr_ww, 1), 1) .= 0 # axial induced velocity

    return calculate_induced_velocities_on_wakes!(
        vz_wake,
        vr_wake,
        vz_ww,
        vr_ww,
        gamw,
        vz_wr,
        vr_wr,
        sigr,
        vz_wb,
        vr_wb,
        gamb;
        post=post,
    )
end

"""
    calculate_induced_velocities_on_wakes!(
        vz_wake,
        vr_wake,
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

In place version of calculate_induced_velocities_on_wakes.
"""
function calculate_induced_velocities_on_wakes!(
    vz_wake,
    vr_wake,
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

    # reset induced velocities on wake
    vz_wake .= 0
    vr_wake .= 0

    # get number of rotors
    nbe, nrotor = size(sigr)

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
    if !isnothing(gamb)
        @views vz_wake .+= vz_wb * gamb
        @views vr_wake .+= vr_wb * gamb

        if post
            @views vzb_wake .+= vz_wb * gamb
            @views vrb_wake .+= vr_wb * gamb
        end
    end

    # add rotor induced velocities
    for jrotor in 1:nrotor
        jrotorrange = (nbe * (jrotor - 1) + 1):(nbe * jrotor)
        @views vz_wake .+= vz_wr[:, jrotorrange] * sigr[:, jrotor]
        @views vr_wake .+= vr_wr[:, jrotorrange] * sigr[:, jrotor]
        if post
            @views vzr_wake .+= vz_wr[:, jrotorrange] * sigr[:, jrotor]
            @views vrr_wake .+= vr_wr[:, jrotorrange] * sigr[:, jrotor]
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

"""
    reframe_wake_velocities!(Cm_wake, vz_wake, vr_wake, Vinf; post=false)

In place version of reframe_wake_velocities.
"""
function reframe_wake_velocities!(Cm_wake, vz_wake, vr_wake, Vinf; post=false)
    Cm_wake .= sqrt.((vz_wake .+ Vinf) .^ 2 .+ vr_wake .^ 2)

    if post
        return vz_wake .+ Vinf, Cm_wake
    else
        # return meridional velocities
        return Cm_wake
    end
end

"""
    calculate_wake_velocities!(Cm_wake, vz_wake, vr_wake, gamw, sigr, gamb, ivw, Vinf)

In place version of calculate_wake_velocities.
"""
function calculate_wake_velocities!(Cm_wake, vz_wake, vr_wake, gamw, sigr, gamb, ivw, Vinf)
    # - Get induced velocities on wake - #
    calculate_induced_velocities_on_wakes!(
        vz_wake,
        vr_wake,
        @view(ivw.v_ww[:, :, 1]),
        @view(ivw.v_ww[:, :, 2]),
        gamw,
        @view(ivw.v_wr[:, :, 1]),
        @view(ivw.v_wr[:, :, 2]),
        sigr,
        @view(ivw.v_wb[:, :, 1]),
        @view(ivw.v_wb[:, :, 2]),
        gamb,
    )

    # - Reframe rotor velocities into blade element frames
    return reframe_wake_velocities!(Cm_wake, vz_wake, vr_wake, Vinf)
end

"""
    get_sheet_jumps(Gamma_tilde, H_tilde) -> Tuple{Matrix{Float64}, Matrix{Float64}}

Compute jumps in circulation and enthalpy across wake sheets shed by the rotors.

This function calculates the differences in cumulative rotor circulation (`Γ̃`) and enthalpy jump (`H̃`) between adjacent radial stations. These jumps are used to define wake sheet strengths.

# Arguments
- `Gamma_tilde::Matrix{Float}` : Net cumulative circulation at each blade element and rotor. Dimensions: (number of blade elements, number of rotors).
- `H_tilde::Matrix{Float}` : Cumulative enthalpy jump at each blade element and rotor. Dimensions: (number of blade elements, number of rotors).

# Returns
- `deltaGamma2::Matrix{Float}` : Circulation difference (jump) between successive elements. Dimensions: (number of wake sheets, number of rotors).
- `deltaH::Matrix{Float}` : Enthalpy difference (jump) between successive elements. Same dimensions as `deltaGamma2`.
"""
function get_sheet_jumps(Gamma_tilde, H_tilde)

    # get problem size
    nr, nrotor = size(Gamma_tilde)
    # number of wakes is one more than number of blade elements
    nw = nr + 1

    # get floating point type
    TF = eltype(Gamma_tilde)

    deltaGamma2 = zeros(TF, nw, nrotor)
    deltaH = zeros(TF, nw, nrotor)

    get_sheet_jumps!(deltaGamma2, deltaH, Gamma_tilde, H_tilde)

    return deltaGamma2, deltaH
end

"""
    get_sheet_jumps!(deltaGamma2, deltaH, Gamma_tilde, H_tilde)

In place function of get_sheet_jumps.
"""
function get_sheet_jumps!(deltaGamma2, deltaH, Gamma_tilde, H_tilde)
    # get problem size
    nw, nrotor = size(deltaGamma2)

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
    average_wake_velocities!(Cm_avg, Cm_wake, nodemap, endnodeidxs)

In place version of average_wake_velocities.
"""
function average_wake_velocities!(Cm_avg, Cm_wake, nodemap, endnodeidxs)

    # Loop through each node, averaging velocities from panel centers
    for iw in eachindex(Cm_avg)

        # Get average meridional velocity at node
        if iw in @view(endnodeidxs[1, :])
            Cm_avg[iw] = Cm_wake[searchsortedfirst(@view(nodemap[1, :]), iw)]
        elseif iw in @view(endnodeidxs[2, :])
            Cm_avg[iw] = Cm_wake[searchsortedfirst(@view(nodemap[2, :]), iw)]
        else
            Cm_avg[iw] =
                0.5 * (
                    Cm_wake[searchsortedfirst(@view(nodemap[1, :]), iw)] +
                    Cm_wake[searchsortedfirst(@view(nodemap[2, :]), iw)]
                )
        end
    end

    return Cm_avg
end

"""
    calculate_wake_vortex_strengths!(
        gamw,
        Gamma_tilde,
        H_tilde,
        deltaGamma2,
        deltaH,
        Gamr,
        Cm_avg,
        B,
        Omega,
        wakeK,
        wake_node_sheet_be_map,
        wake_node_ids_along_casing_wake_interface,
        wake_node_ids_along_center_body_wake_interface;
        post=false,
    )

In place version of calculate_wake_vortex_strengths.
"""
function calculate_wake_vortex_strengths!(
    gamw,
    Gamma_tilde,
    H_tilde,
    deltaGamma2,
    deltaH,
    Gamr,
    Cm_avg,
    B,
    Omega,
    wakeK,
    wake_node_sheet_be_map,
    wake_node_ids_along_casing_wake_interface,
    wake_node_ids_along_center_body_wake_interface;
    post=false,
)

    # get net circulation of upstream rotors
    calculate_net_circulation!(Gamma_tilde, Gamr, B)

    # get enthalpy jump across disks
    calculate_enthalpy_jumps!(H_tilde, Gamr, Omega, B)

    # get the circulation squared and enthalpy jumps across the wake sheets
    get_sheet_jumps!(deltaGamma2, deltaH, Gamma_tilde, H_tilde)

    for (Cm, gw, K, sheetid, rotorid) in zip(
        Cm_avg,
        eachrow(gamw),
        wakeK,
        @view(wake_node_sheet_be_map[:, 1]),
        @view(wake_node_sheet_be_map[:, 2]),
    )

        # calculate the wake vortex strength
        if abs(Cm) < eps()
            # avoid division by zero
            gw[1] = 0.0
        else

            # calculate wake node strength
            gw[1] = -(K * deltaGamma2[sheetid, rotorid] + deltaH[sheetid, rotorid]) / Cm
        end
    end

    # - Wake-Body Interface Treatment - #
    if !isnothing(wake_node_ids_along_casing_wake_interface)
        # gamw[wake_node_ids_along_casing_wake_interface] .= 0.0

        gamw[wake_node_ids_along_casing_wake_interface] =
            gamw[wake_node_ids_along_casing_wake_interface[end] + 1] * (
                1.0 .-
                (
                    wake_node_ids_along_casing_wake_interface[end] .-
                    wake_node_ids_along_casing_wake_interface
                ) / (
                    wake_node_ids_along_casing_wake_interface[end] -
                    wake_node_ids_along_casing_wake_interface[1] + 1
                )
            )
    end
    if !isnothing(wake_node_ids_along_center_body_wake_interface)
        # gamw[wake_node_ids_along_center_body_wake_interface] .= 0.0

        gamw[wake_node_ids_along_center_body_wake_interface] =
            gamw[wake_node_ids_along_center_body_wake_interface[end] + 1] * (
                1.0 .-
                (
                    reverse(wake_node_ids_along_center_body_wake_interface) .-
                    wake_node_ids_along_center_body_wake_interface[1]
                ) / length(wake_node_ids_along_center_body_wake_interface)
            )
    end

    if post
        return gamw, deltaGamma2, deltaH
    else
        return gamw
    end
end
