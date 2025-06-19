"""
    get_body_tangential_velocities(
        gamb,
        gamw,
        sigr,
        ivb,
        Vinf,
        totnode,
        totpanel,
        nnode,
        npanel,
        tangent,
        controlpoints,
        endpanelidxs,
        wake_panel_ids_along_center_body_wake_interface,
        wake_panel_ids_along_casing_wake_interface,
        center_body_panel_ids_along_center_body_wake_interface,
        duct_panel_ids_along_casing_wake_interface,
        num_casing_panels,
    )

Get the tangential velocities along the body surfaces.

# Arguments
- `gamb::Vector{Float}` : the body panel strengths
- `gamw::Vector{Float}` : the wake panel strengths
- `sigr::Vector{Float}` : the rotor panel strengths
- `ivb::NamedTuple` : the unit induced velocities on the bodies
- `Vinf::Vector{Float}` : one element vector containing the freestream magnitude
- `totnode::Int` : total number of nodes between all bodies
- `totpanel::Int` : total number of panels between all bodies
- `nnode::Vector{Int}` : number of nodes in each body
- `npanel::Vector{Int}` : number of panels in each body.
- `tangent::Matrix{Float}` : unit tangent vectors for each panel
- `controlpoints::Matrix{Float}` : control point locations for each panel
- `endpanelidxs::Matrix{Int}` : the indices of the first and last panels for each body
- `wake_panel_ids_along_center_body_wake_interface::Vector{Int}` : the indices of the wake panels coincident with the center_body panels
- `wake_panel_ids_along_casing_wake_interface::Vector{Int}` : the indices of the wake panels coincident with the duct casing (inner surface) panels
- `center_body_panel_ids_along_center_body_wake_interface::Vector{Int}` : the indices of the center_body panels coincident with the wake panels
- `duct_panel_ids_along_casing_wake_interface::Vector{Int}` : the indices of the duct panels coincident with the wake panels
- `num_casing_panels::Int` : the number of panels between the leading and trailing edge of the duct on the duct inner side (casing)

# Returns
- `vtan_tuple::NamedTuple` : a named tuple containing the body tangential surface velocities and various useful breakdowns thereof.
"""
function get_body_tangential_velocities(
    gamb,
    gamw,
    sigr,
    ivb,
    Vinf,
    totnode,
    totpanel,
    nnode,
    npanel,
    tangent,
    controlpoints,
    endpanelidxs,
    wake_panel_ids_along_center_body_wake_interface,
    wake_panel_ids_along_casing_wake_interface,
    center_body_panel_ids_along_center_body_wake_interface,
    duct_panel_ids_along_casing_wake_interface,
    num_casing_panels,
)
    if isnothing(gamw)
        TF = promote_type(eltype(gamb))
    else
        TF = promote_type(eltype(gamb), eltype(gamw), eltype(sigr))
    end

    # - initialize total velocity - #
    Vtot = zeros(TF, 2, totpanel)
    Vtot_in = similar(Vtot) .= 0
    Vtot_out = similar(Vtot) .= 0
    Vtan_in = similar(Vtot, totpanel) .= 0
    Vtan_out = similar(Vtot, totpanel) .= 0
    Vtot_prejump = similar(Vtot) .= 0
    vtot_body = similar(Vtot) .= 0
    duct_jump = similar(Vtot, (npanel[1],))
    center_body_jump = similar(Vtot, (npanel[2],))
    body_jump_term = similar(Vtot, totpanel) .= 0.0
    vtot_jump = similar(Vtot) .= 0.0
    vtot_wake = similar(Vtot) .= 0
    vtot_rotors = similar(Vtot) .= 0.0
    casing_zpts = zeros(TF, num_casing_panels)
    vtan_casing_in = similar(casing_zpts) .= 0
    vtan_casing_out = similar(casing_zpts) .= 0
    nacelle_zpts = zeros(TF, npanel[1] - num_casing_panels)
    vtan_nacelle_in = similar(nacelle_zpts) .= 0
    vtan_nacelle_out = similar(nacelle_zpts) .= 0
    center_body_zpts = zeros(TF, npanel[2])
    vtan_center_body_in = similar(center_body_zpts) .= 0
    vtan_center_body_out = similar(center_body_zpts) .= 0

    vtan_tuple = (;
        # Totals and Components:
        Vtan_in,
        Vtot_in,
        Vtan_out,
        Vtot_out,
        Vtot_prejump,
        vtot_body,
        duct_jump,
        center_body_jump,
        body_jump_term,
        vtot_jump,
        vtot_wake,
        vtot_rotors,
        # Splits:
        vtan_casing_in,
        vtan_casing_out,
        vtan_nacelle_in,
        vtan_nacelle_out,
        vtan_center_body_in,
        vtan_center_body_out,
    )

    return get_body_tangential_velocities!(
        vtan_tuple,
        gamb,
        gamw,
        sigr,
        ivb,
        Vinf,
        totnode,
        totpanel,
        nnode,
        npanel,
        tangent,
        controlpoints,
        endpanelidxs,
        wake_panel_ids_along_center_body_wake_interface,
        wake_panel_ids_along_casing_wake_interface,
        center_body_panel_ids_along_center_body_wake_interface,
        duct_panel_ids_along_casing_wake_interface,
        (; casing_zpts, nacelle_zpts, center_body_zpts),
    )
end

"""
function get_body_tangential_velocities!(
    vtan_tuple,
    gamb,
    gamw,
    sigr,
    ivb,
    Vinf,
    totnode,
    totpanel,
    nnode,
    npanel,
    tangent,
    controlpoints,
    endpanelidxs,
    wake_panel_ids_along_center_body_wake_interface,
    wake_panel_ids_along_casing_wake_interface,
    center_body_panel_ids_along_center_body_wake_interface,
    duct_panel_ids_along_casing_wake_interface,
    zpts,
)

In-place version of `get_body_tangential_velocities`.

# Additional Arguments
- `zpts::NamedTuple` : a named tuple containing the z-coordinates of the control points of the duct casing, duct nacelle, and center_body.
"""
function get_body_tangential_velocities!(
    vtan_tuple,
    gamb,
    gamw,
    sigr,
    ivb,
    Vinf,
    totnode,
    totpanel,
    nnode,
    npanel,
    tangent,
    controlpoints,
    endpanelidxs,
    wake_panel_ids_along_center_body_wake_interface,
    wake_panel_ids_along_casing_wake_interface,
    center_body_panel_ids_along_center_body_wake_interface,
    duct_panel_ids_along_casing_wake_interface,
    zpts,
)

    # - setup - #
    if !isnothing(sigr)
        nws, nrotor = size(sigr)
        (; v_bb, v_br, v_bw) = ivb
    else
        nws = nrotor = nothing
        (; v_bb) = ivb
    end

    # rename for convenience
    hwi = wake_panel_ids_along_center_body_wake_interface
    dwi = wake_panel_ids_along_casing_wake_interface
    whi = center_body_panel_ids_along_center_body_wake_interface
    wdi = duct_panel_ids_along_casing_wake_interface

    (;
        # Totals and Components:
        Vtan_in,
        Vtot_in,
        Vtan_out,
        Vtot_out,
        Vtot_prejump,
        vtot_body,
        duct_jump,
        center_body_jump,
        body_jump_term,
        vtot_jump,
        vtot_wake,
        vtot_rotors,
        # Splits:
        vtan_casing_in,
        vtan_casing_out,
        vtan_nacelle_in,
        vtan_nacelle_out,
        vtan_center_body_in,
        vtan_center_body_out,
    ) = vtan_tuple

    (; casing_zpts, nacelle_zpts, center_body_zpts) = zpts

    # TODO also consider including the body wakes here as well.

    # - Velocity Contributions from body - #
    for (i, vt) in enumerate(eachrow(vtot_body))
        vt .= @views v_bb[:, :, i] * gamb[1:size(v_bb, 2)]
    end
    Vtot_out .+= vtot_body

    # - Velocity Contributions from wake - #
    if !isnothing(gamw)
        for (i, vt) in enumerate(eachrow(vtot_wake))
            vt .= @views v_bw[:, :, i] * gamw
        end
        Vtot_out .+= vtot_wake # opposite sign from linear solve
    end

    if !isnothing(sigr)
        # - Velocity Contributions from rotors - #
        for jrotor in 1:nrotor
            rotorrange = (nws * (jrotor - 1) + 1):(nws * jrotor)
            for (i, vt) in enumerate(eachrow(vtot_rotors))
                vt .+= @views v_br[:, rotorrange, i] * sigr[rotorrange]
            end
        end
        Vtot_out .+= vtot_rotors # opposite sign from linear solve
    end

    # - Influence from Freestream - #
    Vtot_out[1, :] .+= Vinf # opposite sign from linear solve
    Vtot_prejump .= copy(Vtot_out)

    # - Add in Jump Term - #
    # duct
    duct_jump .= @views (gamb[1:(npanel[1])] + gamb[2:(nnode[1])]) / 2

    # wake panels interfacing with duct
    if !isnothing(gamw)
        duct_jump[wdi] .+= @views (
            gamw[dwi[1]:(dwi[end] - 1)] + gamw[(dwi[1] + 1):dwi[end]]
        ) / 2.0
    end

    # center body panels
    center_body_jump .= @views (
        gamb[(nnode[1] + 1):(totnode - 1)] + gamb[(nnode[1] + 2):(totnode)]
    ) / 2.0

    # wake panels interfacing with center body
    if !isnothing(gamw)
        center_body_jump[whi] .+= @views (
            gamw[hwi[1]:(hwi[end] - 1)] + gamw[(hwi[1] + 1):hwi[end]]
        ) / 2.0
    end

    body_jump_term[1:length(duct_jump)] .= duct_jump
    body_jump_term[(length(duct_jump) + 1):end] .= center_body_jump

    for (vt, tan) in zip(eachrow(vtot_jump), eachrow(tangent))
        vt .+= body_jump_term .* tan ./ 2.0
    end

    # assign velocities to each side of the panel
    Vtot_in .= Vtot_out .+ vtot_jump # inner side of boundary
    Vtot_out .-= vtot_jump # outer side of boundary

    # Get the magnitude of the sum of the velocities and this is the surface velocity since the body velocities have been solved to eliminate the normal components in the summed velocities
    Vtan_out .= @views sqrt.(Vtot_out[1, :] .^ 2 .+ Vtot_out[2, :] .^ 2)
    Vtan_in .= @views sqrt.(Vtot_in[1, :] .^ 2 .+ Vtot_in[2, :] .^ 2)

    # - Split Velocities associates with inner and outer duct and hub - #
    # total tangential velocities
    split_bodies!(
        vtan_casing_out,
        vtan_nacelle_out,
        vtan_center_body_out,
        casing_zpts,
        nacelle_zpts,
        center_body_zpts,
        Vtan_out,
        controlpoints,
    )
    split_bodies!(
        vtan_casing_in,
        vtan_nacelle_in,
        vtan_center_body_in,
        casing_zpts,
        nacelle_zpts,
        center_body_zpts,
        Vtan_in,
        controlpoints,
    )

    return vtan_tuple
end

"""
    calculate_vtheta(Gamma_tilde, r)

Calculate tangential velocity for a given net circulation and radial location

# Arguments
- `Gamma_tilde::Matrix{Float}` : Sum of upstream circulation values
- `r::Matrix{Float}` : blade element radial positions

# Returns
- `vtheta::Vector{Float}` : Tangential velocity at each radial location.
"""
function calculate_vtheta(Gamma_tilde, r)
    TF = promote_type(eltype(Gamma_tilde), eltype(r))
    vtheta = zeros(TF, length(r))

    for i in 1:length(r)
        if isapprox(r[i], 0.0)
            vtheta[i] = 0.0
        else
            vtheta[i] = Gamma_tilde ./ (2.0 * pi * r[i])
        end
    end

    return vtheta
end

"""
    calculate_induced_velocities_on_bodywake(
        vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
    )

Calculate the induced velocities on one of the body wakes (unit velocity inputs determine which one)

# Arguments
- `vz_w::Matrix{Float}` : unit axial induced velocity of the wake onto the body wake
- `vr_w::Matrix{Float}` : unit radial induced velocity of the wake onto the body wake
- `gamw::Vector{Float}` : wake panel strengths
- `vz_r::Matrix{Float}` : unit axial induced velocity of the rotor onto the body wake
- `vr_r::Matrix{Float}` : unit radial induced velocity of the rotor onto the body wake
- `sigr::Vector{Float}` : rotor panel strengths
- `vz_b::Matrix{Float}` : unit axial induced velocity of the bodies onto the body wake
- `vr_b::Matrix{Float}` : unit radial induced velocity of the bodies onto the body wake
- `gamb::Vector{Float}` : body panel strengths
- `Vinf::Vector{Float}` : one element vector containing the velocity magnitude

# Returns
- `vz::Vector{Float}` : Total axial velocity at the body wake control points.
- `vr::Vector{Float}` : Total radial velocity at the body wake control points.
"""
function calculate_induced_velocities_on_bodywake(
    vz_w, vr_w, gamw, vz_r, vr_r, sigr, vz_b, vr_b, gamb, Vinf
)

    # problem dimensions
    nrotor = size(sigr, 2) # number of rotors
    np = size(vz_b, 1) # number of panels in bodywake

    # initialize outputs
    vz = Vinf * ones(eltype(gamw), np) # axial induced velocity
    vr = zeros(eltype(gamw), np) # radial induced velocity

    # add body induced velocities
    @views vz[:] .+= vz_b * gamb
    @views vr[:] .+= vr_b * gamb

    # add wake induced velocities
    @views vz[:] .+= vz_w * gamw
    @views vr[:] .+= vr_w * gamw

    # add rotor induced velocities
    for jrotor in 1:nrotor
        @views vz[:] .+= vz_r[jrotor] * sigr[:, jrotor]
        @views vr[:] .+= vr_r[jrotor] * sigr[:, jrotor]
    end

    # return raw induced velocities
    return vz, vr
end
