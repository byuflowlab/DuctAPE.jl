"""

length from size
"""
function lfs(shape)
    if length(shape) == 1
        return shape[1]
    else
        return *(shape...)
    end
end

"""
initialize all the containers you'll want cached for use in the solve, etc.

note that the default chunk size threshold in ForwardDiff is 12, and the automatically chosen chunk size (see pickchunksize function in ForwardDiff) will always pick 12 for the size of cache we are creating here, so we set that as our default.  Will want to do some benchmarking and change that as makes sense depending on machine capabilities.
"""
function initialize_caches(paneling_constants; cache_levels=1, chunksize=12)

    # - Extract Paneling Constants - #
    (; npanels, ncenterbody_inlet, nduct_inlet, wake_length, nwake_sheets, dte_minus_cbte) =
        paneling_constants

    # - Get Problem Dimensions - #

    # number of rotors is one less than the length of npanels if the duct and hub trailing edges line up, and is two less if they don't
    nr = iszero(dte_minus_cbte) ? length(npanels) - 1 : length(npanels) - 2

    # number of wake sheets (blade nodes)
    nws = nwake_sheets

    # number of blade elements (panels)
    nbe = nws - 1

    # number of body panels
    nbp = nduct_inlet * 2 + ncenterbody_inlet
    # add rest of panels mutual between centerbody and duct
    nbp += if iszero(dte_minus_cbte)
        sum(npanels[1:(end - 1)]) * 3
    else
        sum(npanels[1:(end - 2)]) * 3
    end
    # add additional duct or centerbody panels if one extends further back
    if dte_minus_cbte > 0
        nbp += npanels[end - 1] * 2
    elseif dte_minus_cbte < 0
        nbp += npanels[end - 1]
    end

    # number of body nodes is number of panels + number of bodies
    npn = nbp + 2

    # number of wake panels is the total number of npanels times the number of wake sheets
    nwp = sum(npanels) * nwake_sheets

    # number of wake nodes is one more than the number of panels for each wake sheet
    nwn = nwp + nwake_sheets

    # - Initialize Solve Parameter Cache - #
    # TODO: write this function
    # TODO: this cache will include the wake grid and could include stuff used for the state initialization.
    precomp_cache, precomp_cache_dims = generate_precomp_cache()

    # TODO: write this function
    solve_parameter_cache, solve_parameter_cache_dims = generate_solve_parameter_cache()

    # - Initialize Solve Container Cache - #
    # TODO: write this function
    solve_cache, solve_cache_dims = generate_solve_cache()

    # - Initialize Post-process Cache - #
    # TODO: write this function
    post_cache, post_cache_dims = generate_post_cache()

    ######################################################
    #TODO: move everything below here into the appropriate sub function
    ######################################################

    # # EXAMPLE
    # xshape = (3,)
    # xlength = lfs(xshape)
    # xd = (; index=(total_length + 1):(total_length + xlength), shape=xshape)
    # total_length += xlength

    # yshape = (2, 3)
    # ylength = lfs(yshape)
    # yd = (; index=(total_length + 1):(total_length + ylength), shape=yshape)
    # total_length += ylength

    # - initialize - #
    total_length = 0

    ##### ----- Compute Container Sizes ----- #####
    #TODO: go through every item and figure out the shape, and assign indices in a large vector based on the length and what has come before.
    #each item will be a namedtuple with the fields `index` and `shape`
    #then at the end, initialize a DiffCache based on a vector of zeros with the total length of everything.
    #TODO: need to decide what actually goes in what cache

    # - Blade Element Velocities - #
    beshape = (nbe, nr)
    belength = lfs(beshape)

    # absolute
    Wz_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    Wr_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    Wtheta_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    Wm_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    Wmag_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength

    # total relative
    vz_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vr_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vtheta_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength

    # relative components
    vzb_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vrb_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vzw_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vrw_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vzr_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength
    vrr_rotor = (; index=(total_length + 1):(total_length + belength), shape=beshape)
    total_length += belength

    # - influences on body - #
    gamb
    A_bb
    A_bb_LU
    b_bf
    A_bw
    A_pw
    A_br
    A_pr

    # - Assemble/Return NamedTuple - #
    return (;
        cache_vec=pat.DiffCache(ones(total_length), chunksize; levels=cache_levels),
        cache_dims=(;),
    )
end

"""
# precomp cache includes ivr, ivw, linsys, blade_elements, idmaps
"""
function withdraw_solve_parameter_cache(vec, dims)
    # - induced velocities on rotor - #
    ivr = (;
        vz_rb=reshape(vec[dims.ivr.vz_rb.index], dims.ivr.vz_rb.shape),
        vr_rb=reshape(vec[dims.ivr.vr_rb.index], dims.ivr.vr_rb.shape),
        vz_rr=reshape(vec[dims.ivr.vz_rr.index], dims.ivr.vz_rr.shape),
        vr_rr=reshape(vec[dims.ivr.vr_rr.index], dims.ivr.vr_rr.shape),
        vz_rw=reshape(vec[dims.ivr.vz_rw.index], dims.ivr.vz_rw.shape),
        vr_rw=reshape(vec[dims.ivr.vr_rw.index], dims.ivr.vr_rw.shape),
    )

    # - induced velocities on wake - #
    ivw = (;
        vz_wb=reshape(vec[dims.ivr.vz_wb.index], dims.ivr.vz_wb.shape),
        vr_wb=reshape(vec[dims.ivr.vr_wb.index], dims.ivr.vr_wb.shape),
        vz_wr=reshape(vec[dims.ivr.vz_wr.index], dims.ivr.vz_wr.shape),
        vr_wr=reshape(vec[dims.ivr.vr_wr.index], dims.ivr.vr_wr.shape),
        vz_ww=reshape(vec[dims.ivr.vz_ww.index], dims.ivr.vz_ww.shape),
        vr_ww=reshape(vec[dims.ivr.vr_ww.index], dims.ivr.vr_ww.shape),
    )

    # - induced velocities on body - #
    ivb = (;
        vz_bb=reshape(vec[dims.ivr.vz_bb.index], dims.ivr.vz_bb.shape),
        vr_bb=reshape(vec[dims.ivr.vr_bb.index], dims.ivr.vr_bb.shape),
        vz_br=reshape(vec[dims.ivr.vz_br.index], dims.ivr.vz_br.shape),
        vr_br=reshape(vec[dims.ivr.vr_br.index], dims.ivr.vr_br.shape),
        vz_bw=reshape(vec[dims.ivr.vz_bw.index], dims.ivr.vz_bw.shape),
        vr_bw=reshape(vec[dims.ivr.vr_bw.index], dims.ivr.vr_bw.shape),
    )

    # - linear system - #
    linsys = (;
        A_bb=reshape(vec[dims.linsys.A_bb.index], dims.linsys.A_bb.shape),
        A_bb_LU=la.LU(
            reshape(view(vec, dims.linsys.A_bb_LU.index), dims.linsys.A_bb_LU.shape),
            [i for i in 1:dims.linsys.A_bb_LU.shape[1]],
            0,
        ),
        b_bf=reshape(vec[dims.linsys.b_bf.index], dims.linsys.b_bf.shape),
        A_bw=reshape(vec[dims.linsys.A_bw.index], dims.linsys.A_bw.shape),
        A_pw=reshape(vec[dims.linsys.A_pw.index], dims.linsys.A_pw.shape),
        A_br=reshape(vec[dims.linsys.A_br.index], dims.linsys.A_br.shape),
        A_pr=reshape(vec[dims.linsys.A_pr.index], dims.linsys.A_pr.shape),
    )

    # - blade element geometry - #
    blade_elements(;
        B=reshape(vec[dims.blade_elements.B.index], dims.blade_elements.B.shape),
        Omega=reshape(
            vec[dims.blade_elements.Omega.index], dims.blade_elements.Omega.shape
        ),
        fliplift=reshape(
            vec[dims.blade_elements.fliplift.index], dims.blade_elements.fliplift.shape
        ),
        chords=reshape(
            vec[dims.blade_elements.chords.index], dims.blade_elements.chords.shape
        ),
        twists=reshape(
            vec[dims.blade_elements.twists.index], dims.blade_elements.twists.shape
        ),
        stagger=reshape(
            vec[dims.blade_elements.stagger.index], dims.blade_elements.stagger.shape
        ),
        solidity=reshape(
            vec[dims.blade_elements.solidity.index], dims.blade_elements.solidity.shape
        ),
        rotor_panel_center=reshape(
            vec[dims.blade_elements.rotor_panel_center.index],
            dims.blade_elements.rotor_panel_center.shape,
        ),
        inner_fraction=reshape(
            vec[dims.blade_elements.inner_fraction.index],
            dims.blade_elements.inner_fraction.shape,
        ),
    )

    # - index maps - #
    idmaps = (;
        ductwakeinterfacenodeid=reshape(
            vec[dims.idmaps.ductwakeinterfacenodeid.index],
            dims.idmaps.ductwakeinterfacenodeid.shape,
        ),
        hubwakeinterfacenodeid=reshape(
            vec[dims.idmaps.hubwakeinterfacenodeid.index],
            dims.idmaps.hubwakeinterfacenodeid.shape,
        ),
        #TODO: probably don't need this one?
        #TODO: or it needs to be elsewhere, it's only used for the aero initializaiton
        # rotorwakepanelid=reshape(
        #     vec[dims.idmaps.rotorwakepanelid.index], dims.idmaps.rotorwakepanelid.shape
        # ),
        rotorwakeid=reshape(
            vec[dims.idmaps.rotorwakeid.index], dims.idmaps.rotorwakeid.shape
        ),
    )

    return (; ivr, ivw, linsys, blade_elements, idmaps)
end

"""
"""
function withdraw_post_parameter_cache(vec, dims)
    idmaps = (;
        ductidsaftofrotors=(), wakeductinterfacepanelid=(), wakehubinterfacepanelid=()
    )
    return (;)
end

"""
"""
function withdraw_solve_container_cache(vec, dims)
    return (;
        # Strengths
        gamb=reshape(vec[dims.gamb.index], dims.gamb.shape),
        Gamr=reshape(vec[dims.Gamr.index], dims.Gamr.shape),
        sigr=reshape(vec[dims.sigr.index], dims.sigr.shape),
        gamw=reshape(vec[dims.gamw.index], dims.gamw.shape),

        # Blade Element Values
        beta1=reshape(vec[dims.beta1.index], dims.beta1.shape),
        Cmag_rotor=reshape(vec[dims.Cmag_rotor.index], dims.Cmag_rotor.shape),
        cl=reshape(vec[dims.cl.index], dims.cl.shape),
        cd=reshape(vec[dims.cd.index], dims.cd.shape),

        # Wake Velocities
        vz_wake=reshape(vec[dims.vz_wake.index], dims.vz_wake.shape),
        vr_wake=reshape(vec[dims.vr_wake.index], dims.vr_wake.shape),

        # State estimates
        Vz_est=reshape(vec[dims.Vz_est.index], dims.Vz_est.shape),
        Vtheta_est=reshape(vec[dims.Vtheta_est.index], dims.Vtheta_est.shape),
        Cm_est=reshape(vec[dims.Cm_est.index], dims.Cm_est.shape),
    )
end

"""
zero out all containers used inside nlsolve residual (velocities, gamb, blade element loads and inflow angles, etc.)
"""
function reset_solve_containers!(c)

    # - Set all the relevant containers to zero - #

    # Strengths
    c.gamb .= 0
    c.Gamr .= 0
    c.sigr .= 0
    c.gamw .= 0

    # Blade Element Values
    c.beta1 .= 0
    c.Cmag_rotor .= 0
    c.cl .= 0
    c.cd .= 0

    # Wake Velocities
    c.vz_wake .= 0
    c.vr_wake .= 0

    # State estimates
    c.Vz_est .= 0
    c.Vtheta_est .= 0
    c.Cm_est .= 0

    return c
end
