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
function allocate_solve_parameter_cache(
    solve_type::NewtonSolve, paneling_constants; fd_chunk_size=12, levels=2
)

    # - Get problem dimensions - #
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

    # - Initial Guesses - #
    s = (nrotor * nbe,)
    l = lfs(s)
    vz_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    vtheta_rotor = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp,)
    l = lfs(s)
    Cm_wake = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Operating Point - #
    s = (1,)
    l = lfs(s)
    Vinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rhoinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    muinf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    asound = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor,)
    l = lfs(s)
    Omega = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Rotors - #
    s = (nrotor * nbe, nbn, 2)
    l = lfs(s)
    v_rb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nrotor * nws, 2)
    l = lfs(s)
    v_rr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nrotor * nbe, nwn, 2)
    l = lfs(s)
    v_rw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Induced Velocities on Wakes - #
    s = (nwp, nbn, 2)
    l = lfs(s)
    v_wb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nrotor * nws, 2)
    l = lfs(s)
    v_wr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nwp, nwn, 2)
    l = lfs(s)
    v_ww = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # # - Induced Velocities on Bodies - #
    # TODO: move this to the post-processing cache
    # s = (nbp, nbn, 2)
    # l = lfs(s)
    # v_bb = (; index=(total_length + 1):(total_length + l), shape=s)
    # total_length += l

    # s = (nbp, nrotor * nws, 2)
    # l = lfs(s)
    # v_br = (; index=(total_length + 1):(total_length + l), shape=s)
    # total_length += l

    # s = (nbp, nwn, 2)
    # l = lfs(s)
    # v_bw = (; index=(total_length + 1):(total_length + l), shape=s)
    # total_length += l

    # - Linear System - #

    s = (nbn + 2, nbn + 2)
    l = lfs(s)
    A_bb = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    A_bb_LU = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbn + 2,)
    l = lfs(s)
    b_bf = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nwn)
    l = lfs(s)
    A_bw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nwn)
    l = lfs(s)
    A_pw = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbp, nrotor * nws)
    l = lfs(s)
    A_br = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (2, nrotor * nws)
    l = lfs(s)
    A_pr = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (1,)
    l = lfs(s)
    lu_decomp_flag = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Blade Elements - #
    s = (nrotor,)
    l = lfs(s)
    B = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rtip = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    Rhub = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    fliplift = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    s = (nbe, nrotor)
    l = lfs(s)
    chords = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    twists = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    stagger = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    solidity = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    rotor_panel_centers = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    inner_fraction = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # - Wake constants - #
    s = (nwn,)
    l = lfs(s)
    wakeK = (; index=(total_length + 1):(total_length + l), shape=s)
    total_length += l

    # TODO: WHAT TO DO ABOUT AIRFOILS?
    # NOTE: how to know the size of airfoil objects up front?  if it's a DFDC type, this is pretty simple, but if it's a CCBlade type, then the size could be anything depending on how many data points are used
    # airfoils = allocate_airfoil_cache(aftype, nrotor, nbe)

    # # - Index Maps - #
    # # TODO: decide if you actually need these here. these actually won't/can't change at the optimization level, so you could have a different cache for these
    # # TODO: figure out what you renamed these to.  Also probably need to just call the function to get the sizes and lengths for some of these
    # wake_node_ids_along_casing_wake_interface
    # wake_node_ids_along_centerbody_wake_interface
    # rotorwakenodeid
    # wake_nodemap
    # wake_endnodeidxs
    # rotor_indices_in_wake
    # body_totnodes

    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_parameter_cache_dims=(;
            vz_rotor,
            vtheta_rotor,
            Cm_wake,
            operating_point=(; Vinf, rhoinf, muinf, asound, Omega),
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(; v_wb, v_wr, v_ww),
            # ivb=(; v_bb, v_br, v_bw),
            linsys=(; A_bb, A_bb_LU, b_bf, A_bw, A_pw, A_br, A_pr, lu_decomp_flag),
            blade_elements=(;
                B,
                fliplift,
                chords,
                twists,
                stagger,
                solidity,
                rotor_panel_centers,
                inner_fraction,
                Rtip,
                Rhub,
            ),
            wakeK,
            # idmaps=(;),
        ),
    )
end

"""
TODO: add another version of this that sets up cache for DFDC-like CSOR solve
"""
function allocate_solve_container_cache(
    solve_type::NewtonSolve, paneling_constants; fd_chunk_size=12, levels=1
)
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
    rhs = (; index=(total_length + 1):(total_length + l), shape=s)
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
        solve_container_cache=PreallocationTools.DiffCache(
            zeros(total_length), fd_chunk_size; levels=levels
        ),
        solve_container_cache_dims=(;
            gamb,
            rhs,
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

    # - Initial Guesses - #
    vz_rotor = reshape(@view(vec[dims.vz_rotor.index]), dims.vz_rotor.shape)
    vtheta_rotor = reshape(@view(vec[dims.vtheta_rotor.index]), dims.vtheta_rotor.shape)
    Cm_wake = reshape(@view(vec[dims.Cm_wake.index]), dims.Cm_wake.shape)

    # - Operating Point - #
    operating_point = (;
        Vinf=reshape(
            @view(vec[dims.operating_point.Vinf.index]), dims.operating_point.Vinf.shape
        ),
        rhoinf=reshape(
            @view(vec[dims.operating_point.rhoinf.index]), dims.operating_point.rhoinf.shape
        ),
        muinf=reshape(
            @view(vec[dims.operating_point.muinf.index]), dims.operating_point.muinf.shape
        ),
        asound=reshape(
            @view(vec[dims.operating_point.asound.index]), dims.operating_point.asound.shape
        ),
        Omega=reshape(
            @view(vec[dims.operating_point.Omega.index]), dims.operating_point.Omega.shape
        ),
    )

    # - induced velocities on rotor - #
    ivr = (;
        v_rb=reshape(@view(vec[dims.ivr.v_rb.index]), dims.ivr.v_rb.shape),
        v_rr=reshape(@view(vec[dims.ivr.v_rr.index]), dims.ivr.v_rr.shape),
        v_rw=reshape(@view(vec[dims.ivr.v_rw.index]), dims.ivr.v_rw.shape),
    )

    # - induced velocities on wake - #
    ivw = (;
        v_wb=reshape(@view(vec[dims.ivw.v_wb.index]), dims.ivw.v_wb.shape),
        v_wr=reshape(@view(vec[dims.ivw.v_wr.index]), dims.ivw.v_wr.shape),
        v_ww=reshape(@view(vec[dims.ivw.v_ww.index]), dims.ivw.v_ww.shape),
    )

    # # - induced velocities on body - #
    # TODO: move to post process withdraw
    # ivb = (;
    #     v_bb=reshape(@view(vec[dims.ivr.v_bb.index]), dims.ivr.v_bb.shape),
    #     v_br=reshape(@view(vec[dims.ivr.v_br.index]), dims.ivr.v_br.shape),
    #     v_bw=reshape(@view(vec[dims.ivr.v_bw.index]), dims.ivr.v_bw.shape),
    # )

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
        lu_decomp_flag=reshape(@view(vec[dims.linsys.lu_decomp_flag.index]), dims.linsys.lu_decomp_flag.shape),
    )

    # - blade element geometry - #
    blade_elements = (;
        B=reshape(@view(vec[dims.blade_elements.B.index]), dims.blade_elements.B.shape),
        Rtip=reshape(@view(vec[dims.blade_elements.Rtip.index]), dims.blade_elements.Rtip.shape),
        Rhub=reshape(@view(vec[dims.blade_elements.Rhub.index]), dims.blade_elements.Rhub.shape),
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
        rotor_panel_centers=reshape(
            @view(vec[dims.blade_elements.rotor_panel_centers.index]),
            dims.blade_elements.rotor_panel_centers.shape,
        ),
        inner_fraction=reshape(
            @view(vec[dims.blade_elements.inner_fraction.index]),
            dims.blade_elements.inner_fraction.shape,
        ),
    )

    wakeK = reshape(@view(vec[dims.wakeK.index]), dims.wakeK.shape)

    # # - index maps - #
    # idmaps = (;
    #     wake_node_ids_along_casing_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_casing_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_casing_wake_interface.shape,
    #     ),
    #     wake_node_ids_along_centerbody_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_centerbody_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_centerbody_wake_interface.shape,
    #     ),
    #     id_of_first_casing_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.shape,
    #     ),
    #     id_of_first_centerbody_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.shape,
    #     ),
    #     rotorwakenodeid=reshape(
    #         @view(vec[dims.idmaps.rotorwakenodeid.index]), dims.idmaps.rotorwakenodeid.shape
    #     ),
    #     wake_nodemap=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    #     wake_endnodeidxs=reshape(
    #         @view(vec[dims.idmaps.wake_endnodeidxs.index]),
    #         dims.idmaps.wake_endnodeidxs.shape,
    #     ),
    #     rotor_indices_in_wake=reshape(
    #         @view(vec[dims.idmaps.rotor_indices_in_wake.index]),
    #         dims.idmaps.rotor_indices_in_wake.shape,
    #     ),
    #     body_totnodes=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    # )

    return (;
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
    )
end

"""
"""
function withdraw_solve_container_cache(vec, dims)
    return (;
        # Strengths
        gamb=reshape(@view(vec[dims.gamb.index]), dims.gamb.shape),
        rhs=reshape(@view(vec[dims.rhs.index]), dims.rhs.shape),
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
        ductidsaftofrotors=(),
        wake_panel_ids_along_casing_wake_interface=(),
        wake_panel_ids_along_centerbody_wake_interface=(),
    )
    return (;)
end

