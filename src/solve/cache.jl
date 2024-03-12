
"""
length from size
move to utilities
"""
function lfs(shape)
    if length(shape) == 1
        return shape[1]
    else
        return *(shape...)
    end
end

"""
note: containers must be Arrays, structs of arrays, or tuples of arrays
move to utilities
"""
function reset_containers!(c)
    if typeof(c) <: AbstractArray
        #do nothing if it's a string
        (eltype(c) == String) || (c .= 0)
    else
        for p in propertynames(c)
            cp = getfield(c, p)
            if typeof(cp) <: AbstractArray
                #do nothing if it's a string
                (eltype(cp) == String) || (cp .= 0)
            else
                reset_containers!(cp)
            end
        end
    end

    return c
end

######################################################################
#                                                                    #
#                         CACHE ALLOCATIONS                          #
#                                                                    #
######################################################################

"""
"""
function get_problem_dimensions(paneling_constants)

    # - Extract Paneling Constants - #
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

    # number of rotors is one less than the length of npanels if the duct and hub trailing edges line up, and is two less if they don't
    nrotor = iszero(dte_minus_cbte) ? length(npanels) - 1 : length(npanels) - 2

    # number of wake sheets (blade nodes)
    nws = nwake_sheets

    # number of blade elements (panels)
    nbe = nws - 1

    # number of body panels
    # TODO: need number of duct and centerbody nodes separately as well
    ndp = nduct_inlet * 2
    ncbp = ncenterbody_inlet
    # add rest of panels mutual between centerbody and duct
    if iszero(dte_minus_cbte)
        ndp += sum(npanels[1:(end - 1)]) * 2
        ncbp += sum(npanels[1:(end - 1)])
    else
        ndp += sum(npanels[1:(end - 2)]) * 2
        ncbp += sum(npanels[1:(end - 2)])
    end

    # add additional duct or centerbody panels if one extends further back
    if dte_minus_cbte > 0
        ndp += npanels[end - 1] * 2
    elseif dte_minus_cbte < 0
        ncbp += npanels[end - 1]
    end

    # duct and center body nodes are 1 more than number of panels
    ndn = ndp + 1
    ncbn = ncbp + 1

    # number of body panels is sum of duct and centerbody panels
    nbp = ndp + ncbp
    # number of body nodes is number of panels + number of bodies
    nbn = ndn + ncbn

    # number of nodes in each wake sheet
    nwsp = sum(npanels)
    nwsn = nwsp + 1

    # number of wake panels is the total number of npanels times the number of wake sheets
    nwp = sum(npanels) * nwake_sheets

    # number of wake nodes is one more than the number of panels for each wake sheet
    nwn = nwp + nwake_sheets

    # number of duct-wake and centerbody-wake interface nodes
    if iszero(dte_minus_cbte)
        ndwin = sum(npanels[1:(end - 1)]) + 1
        ncbwin = sum(npanels[end - 1]) + 1
    elseif dte_minus_cbte < 0
        ndwin = sum(npanels[1:(end - 2)]) + 1
        ncbwin = sum(npanels[end - 1]) + 1
    else
        ndwin = sum(npanels[1:(end - 1)]) + 1
        ncbwin = sum(npanels[end - 2]) + 1
    end

    return (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes in each wake sheet
        nwsp,   # number of panels in each wake sheet
        ndwin,  # number of duct-wake interfacing nodes
        ncbwin, # number of centerbody-wake interfacing nodes
    )
end

"""
"""
function allocate_precomp_container_cache(paneling_constants)
    pd = get_problem_dimensions(paneling_constants)

    (;
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nwsn,   # number of nodes per wake sheet
    ) = pd

    # - initialize - #
    total_length = 0

    # - Set up Dimensions - #
    dcshape = (2, ndn)
    dclength = lfs(dcshape)
    rp_duct_coordinates = (;
        index=(total_length + 1):(total_length + dclength), shape=dcshape
    )
    total_length += dclength

    cbcshape = (2, ncbn)
    cbclength = lfs(cbcshape)
    rp_centerbody_coordinates = (;
        index=(total_length + 1):(total_length + cbclength), shape=cbcshape
    )
    total_length += cbclength

    wgshape = (2, nwsn, nws)
    wglength = lfs(wgshape)
    wake_grid = (; index=(total_length + 1):(total_length + wglength), shape=wgshape)
    total_length += wglength

    aicnshape = (nbp, nbn)
    aicnlength = lfs(aicnshape)
    AICn = (; index=(total_length + 1):(total_length + aicnlength), shape=aicnshape)
    total_length += aicnlength

    aicpcpshape = (nbp, nbn)
    aicpcplength = lfs(aicpcpshape)
    AICpcp = (; index=(total_length + 1):(total_length + aicpcplength), shape=aicpcpshape)
    total_length += aicpcplength

    # return tuple of initialized cache and associated dimensions
    return (;
        precomp_container_cache=PreallocationTools.DiffCache(zeros(total_length)),
        precomp_container_cache_dims=(;
            rp_duct_coordinates, rp_centerbody_coordinates, wake_grid, AICn, AICpcp
        ),
    )
end

"""
"""
function allocate_solve_parameter_cache(paneling_constants)
    pd = get_problem_dimensions(paneling_constants)

    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes per wake sheet
        ndwin,  # number of duct-wake interfacing nodes
        ncbwin, # number of centerbody-wake interfacing nodes
    ) = pd

    # - initialize - #
    total_length = 0

    ############################################################
    ############################################################
    ##### ----- TODO; YOUR ARE HERE ----- ######################
    ############################################################
    ############################################################

    # - Induced Velocities on Rotors - #
    v_rb
    v_rr
    v_rw

    # - Induced Velocities on Wakes - #
    v_wb
    v_wr
    v_ww

    # - Induced Velocities on Bodies - #
    v_bb
    v_br
    v_bw

    # - Linear System - #
    A_bb
    A_bb_LU
    b_bf
    A_bw
    A_pw
    A_br
    A_pr

    # - Blade Elements - #
    B
    Omega
    fliplift
    chords
    twists
    stagger
    solidity
    rotor_panel_center
    inner_fraction
    # TODO: WHAT TO DO ABOUT AIRFOILS?

    # - Index Maps - #
    wake_node_ids_along_casing_wake_interface
    wake_node_ids_along_centerbody_wake_interface
    rotorwakenodeid
    wake_nodemap
    wake_endnodeidxs
    rotor_indices_in_wake
    body_totnodes

    #TODO; decide if this can be all floats or if you actually need to use the diff cache stuff.
    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(zeros(total_length)),
        solve_parameter_cache_dims=(;
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(;),
            ivb=(;),
            linsys=(; A_bb, A_bb_LU, b_bf, A_bw, A_pw, A_br, A_pr),
            blade_elements=(;),
            idmaps=(;),
        ),
    )
end

"""
"""
function allocate_solve_container_cache(paneling_constants)
    pd = get_problem_dimensions(paneling_constants)

    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbe,    # number of blade elements (also rotor panels)
    ) = pd

    # - initialize - #
    total_length = 0

    # Strengths
    # TODO: is this always going to be 2?, may want to add nb for nbodies to problem dims and then do +nb
    s = (nbn + 2,)
    l = lfs(s)
    gamb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    Gamr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe + 1, nrotor)
    l = lfs(s)
    sigr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    gamw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Blade Element Values
    s = (nbe, nrotor)
    l = lfs(s)
    beta1 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Ctheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Cmag_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cl = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    cd = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    alpha = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    reynolds = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    mach = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Circulation
    Gamma_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    H_tilde = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltaGamma2 = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    deltaH = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # Wake Velocities
    s = (nwp,)
    l = lfs(s)
    vz_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vr_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # State estimates
    s = (nbe, nrotor)
    l = lfs(s)
    vz_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtheta_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp,)
    l = lfs(s)
    Cm_est = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwn,)
    l = lfs(s)
    Cm_avg = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    return (;
        solve_container_cache=PreallocationTools.DiffCache(zeros(total_length)),
        solve_container_cache_dims=(;
            gamb,
            Gamr,
            sigr,
            gamw,
            beta1,
            alpha,
            reynolds,
            mach,
            Cz_rotor,
            Ctheta_rotor,
            Cmag_rotor,
            cl,
            cd,
            vz_wake,
            vr_wake,
            vz_est,
            vtheta_est,
            Cm_est,
            Cm_avg,
            Gamma_tilde,
            H_tilde,
            deltaGamma2,
            deltaH,
        ),
    )
end

"""
"""
function allocate_post_cache(paneling_constants)
    pd = get_problem_dimensions(paneling_constants)

    (;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ndn,    # number of duct nodes
        ncbn,   # number of centerbody nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes per wake sheet
        ndwin,  # number of duct-wake interfacing nodes
        ncbwin, # number of centerbody-wake interfacing nodes
    ) = pd

    # - initialize - #
    total_length = 0

    return (;
        post_cache=PreallocationTools.DiffCache(zeros(total_length)), post_cache_dims=(;)
    )
end

######################################################################
#                                                                    #
#                          CACHE RESHAPING                           #
#                                                                    #
######################################################################

"""
"""
function withdraw_precomp_container_cache(vec, dims)
    rp_duct_coordinates = reshape(
        @view(vec[dims.rp_duct_coordinates.index]), dims.rp_duct_coordinates.shape
    )
    rp_centerbody_coordinates = reshape(
        @view(vec[dims.rp_centerbody_coordinates.index]),
        dims.rp_centerbody_coordinates.shape,
    )
    wake_grid = reshape(@view(vec[dims.wake_grid.index]), dims.wake_grid.shape)
    AICn = reshape(@view(vec[dims.AICn.index]), dims.AICn.shape)
    AICpcp = reshape(@view(vec[dims.AICpcp.index]), dims.AICpcp.shape)

    return (; rp_duct_coordinates, rp_centerbody_coordinates, wake_grid, AICn, AICpcp)
end

"""
"""
function withdraw_solve_parameter_cache(vec, dims)
    # - induced velocities on rotor - #
    ivr = (;
        v_rb=reshape(@view(vec[dims.ivr.v_rb.index]), dims.ivr.v_rb.shape),
        v_rr=reshape(@view(vec[dims.ivr.v_rr.index]), dims.ivr.v_rr.shape),
        v_rw=reshape(@view(vec[dims.ivr.v_rw.index]), dims.ivr.v_rw.shape),
    )

    # - induced velocities on wake - #
    ivw = (;
        v_wb=reshape(@view(vec[dims.ivr.v_wb.index]), dims.ivr.v_wb.shape),
        v_wr=reshape(@view(vec[dims.ivr.v_wr.index]), dims.ivr.v_wr.shape),
        v_ww=reshape(@view(vec[dims.ivr.v_ww.index]), dims.ivr.v_ww.shape),
    )

    # - induced velocities on body - #
    ivb = (;
        v_bb=reshape(@view(vec[dims.ivr.v_bb.index]), dims.ivr.v_bb.shape),
        v_br=reshape(@view(vec[dims.ivr.v_br.index]), dims.ivr.v_br.shape),
        v_bw=reshape(@view(vec[dims.ivr.v_bw.index]), dims.ivr.v_bw.shape),
    )

    # - linear system - #
    linsys = (;
        A_bb=reshape(@view(vec[dims.linsys.A_bb.index]), dims.linsys.A_bb.shape),
        A_bb_LU=LinearAlgebra.LU(
            reshape(@view(vec[dims.linsys.A_bb_LU.index]), dims.linsys.A_bb_LU.shape),
            [i for i in 1:dims.linsys.A_bb_LU.shape[1]],
            0,
        ),
        b_bf=reshape(@view(vec[dims.linsys.b_bf.index]), dims.linsys.b_bf.shape),
        A_bw=reshape(@view(vec[dims.linsys.A_bw.index]), dims.linsys.A_bw.shape),
        A_pw=reshape(@view(vec[dims.linsys.A_pw.index]), dims.linsys.A_pw.shape),
        A_br=reshape(@view(vec[dims.linsys.A_br.index]), dims.linsys.A_br.shape),
        A_pr=reshape(@view(vec[dims.linsys.A_pr.index]), dims.linsys.A_pr.shape),
    )

    # - blade element geometry - #
    blade_elements(;
        B=reshape(@view(vec[dims.blade_elements.B.index]), dims.blade_elements.B.shape),
        fliplift=reshape(
            @view(vec[dims.blade_elements.fliplift.index]),
            dims.blade_elements.fliplift.shape,
        ),
        chords=reshape(
            @view(vec[dims.blade_elements.chords.index]), dims.blade_elements.chords.shape
        ),
        twists=reshape(
            @view(vec[dims.blade_elements.twists.index]), dims.blade_elements.twists.shape
        ),
        stagger=reshape(
            @view(vec[dims.blade_elements.stagger.index]), dims.blade_elements.stagger.shape
        ),
        solidity=reshape(
            @view(vec[dims.blade_elements.solidity.index]),
            dims.blade_elements.solidity.shape,
        ),
        rotor_panel_center=reshape(
            @view(vec[dims.blade_elements.rotor_panel_center.index]),
            dims.blade_elements.rotor_panel_center.shape,
        ),
        inner_fraction=reshape(
            @view(vec[dims.blade_elements.inner_fraction.index]),
            dims.blade_elements.inner_fraction.shape,
        ),
    )

    # - index maps - #
    idmaps = (;
        wake_node_ids_along_casing_wake_interface=reshape(
            @view(vec[dims.idmaps.wake_node_ids_along_casing_wake_interface.index]),
            dims.idmaps.wake_node_ids_along_casing_wake_interface.shape,
        ),
        wake_node_ids_along_centerbody_wake_interface=reshape(
            @view(vec[dims.idmaps.wake_node_ids_along_centerbody_wake_interface.index]),
            dims.idmaps.wake_node_ids_along_centerbody_wake_interface.shape,
        ),
        id_of_first_casing_panel_aft_of_each_rotor=reshape(
            @view(vec[dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.index]),
            dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.shape,
        ),
        id_of_first_centerbody_panel_aft_of_each_rotor=reshape(
            @view(vec[dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.index]),
            dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.shape,
        ),
        rotorwakenodeid=reshape(
            @view(vec[dims.idmaps.rotorwakenodeid.index]), dims.idmaps.rotorwakenodeid.shape
        ),
        wake_nodemap=reshape(
            @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
        ),
        wake_endnodeidxs=reshape(
            @view(vec[dims.idmaps.wake_endnodeidxs.index]),
            dims.idmaps.wake_endnodeidxs.shape,
        ),
        rotor_indices_in_wake=reshape(
            @view(vec[dims.idmaps.rotor_indices_in_wake.index]),
            dims.idmaps.rotor_indices_in_wake.shape,
        ),
        body_totnodes=reshape(
            @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
        ),
    )

    return (; ivr, ivw, linsys, blade_elements, idmaps)
end

"""
"""
function withdraw_solve_container_cache(vec, dims)
    return (;
        # Strengths
        gamb=reshape(@view(vec[dims.gamb.index]), dims.gamb.shape),
        Gamr=reshape(@view(vec[dims.Gamr.index]), dims.Gamr.shape),
        sigr=reshape(@view(vec[dims.sigr.index]), dims.sigr.shape),
        gamw=reshape(@view(vec[dims.gamw.index]), dims.gamw.shape),

        # Blade Element Values
        Cz_rotor=reshape(@view(vec[dims.Cz_rotor.index]), dims.Cz_rotor.shape),
        Ctheta_rotor=reshape(@view(vec[dims.Ctheta_rotor.index]), dims.Ctheta_rotor.shape),
        Cmag_rotor=reshape(@view(vec[dims.Cmag_rotor.index]), dims.Cmag_rotor.shape),
        cl=reshape(@view(vec[dims.cl.index]), dims.cl.shape),
        cd=reshape(@view(vec[dims.cd.index]), dims.cd.shape),
        beta1=reshape(@view(vec[dims.beta1.index]), dims.beta1.shape),
        alpha=reshape(@view(vec[dims.alpha.index]), dims.alpha.shape),
        reynolds=reshape(@view(vec[dims.reynolds.index]), dims.reynolds.shape),
        mach=reshape(@view(vec[dims.mach.index]), dims.mach.shape),

        # Wake Velocities
        vz_wake=reshape(@view(vec[dims.vz_wake.index]), dims.vz_wake.shape),
        vr_wake=reshape(@view(vec[dims.vr_wake.index]), dims.vr_wake.shape),
        Cm_avg=reshape(@view(vec[dims.Cm_avg.index]), dims.Cm_avg.shape),
        Gamma_tilde=reshape(@view(vec[dims.Gamma_tilde.index]), dims.Gamma_tilde.shape),
        H_tilde=reshape(@view(vec[dims.H_tilde.index]), dims.H_tilde.shape),
        deltaGamma2=reshape(@view(vec[dims.deltaGamma2.index]), dims.deltaGamma2.shape),
        deltaH=reshape(@view(vec[dims.deltaH.index]), dims.deltaH.shape),

        # State estimates
        vz_est=reshape(@view(vec[dims.vz_est.index]), dims.vz_est.shape),
        vtheta_est=reshape(@view(vec[dims.vtheta_est.index]), dims.vtheta_est.shape),
        Cm_est=reshape(@view(vec[dims.Cm_est.index]), dims.Cm_est.shape),
    )
end

"""
"""
function withdraw_post_parameter_cache(vec, dims)
    idmaps = (;
        ductidsaftofrotors=(), wake_panel_ids_along_casing_wake_interface=(), wake_panel_ids_along_centerbody_wake_interface=()
    )
    return (;)
end

