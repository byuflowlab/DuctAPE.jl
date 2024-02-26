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
function initialize_cache(panels; cache_levels=1, chunksize=12)

    #TODO: will need some more inputs to know the sizes of everything required.

    # - Extract Propulsor Parameters - #
    (; body_vortex_panels, rotor_source_panels, wake_vortex_panels) = panels

    # rename for brevity
    b = body_vortex_panels
    r = rotor_source_panels
    w = wake_vortex_panels

    # - Get Problem Dimensions - #
    nr = length(r)       # number of rotors
    nbe = r[1].totpanel  # number of blade elements (panels)
    nws = r[1].totnode   # number of wake sheets (blade nodes)
    nbn = b.totnode      # number of body nodes
    nbp = b.totpanel     # number of body panels
    nwn = w.totnode      # number of wake nodes
    nwp = w.totpanel     # number of wake panels

    # - initialize - #
    total_length = 0

    # # EXAMPLE
    # xshape = (3,)
    # xlength = lfs(xshape)
    # xd = (; index=(total_length + 1):(total_length + xlength), shape=xshape)
    # total_length += xlength

    # yshape = (2, 3)
    # ylength = lfs(yshape)
    # yd = (; index=(total_length + 1):(total_length + ylength), shape=yshape)
    # total_length += ylength

    ##### ----- Compute Container Sizes ----- #####
    #TODO: go through every item and figure out the shape, and assign indices in a large vector based on the length and what has come before.
    #each item will be a namedtuple with the fields `index` and `shape`
    #then at the end, initialize a DiffCache based on a vector of zeros with the total length of everything.

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
    b_bf
    A_bw
    A_pw
    A_br
    A_pr
    LHS

    # - Assemble/Return NamedTuple - #
    return (;
        cache_vec=pat.DiffCache(ones(total_length), chunksize; levels=cache_levels),
        cache_dims=(;),
    )
end

"""
"""
function withdraw_cache(vec, dims)
    # TODO: need to get all the same names as contained in dims then do something like:
    # varname = reshape(vec[dims.varname.index], dims.varname.shape)
    # then return a namedtuple of all the variables
    # return (;varname, )
end

"""
Calculate all the AIC matrices and various velocity stuff
"""
function populate_cache!(cache, propulsor)
    # - Extract Fan Parameters - #

    # - Compute "Everything" (in place) - #

    return cache
end

"""
zero out all containers used inside nlsolve residual (velocities, gamb, blade element loads and inflow angles, etc.)
"""
function reset_containers!(c)

    # - Set all the relevant containers to zero - #

    # Strengths
    c.gamb .= 0
    c.Gamr .= 0
    c.sigr .= 0
    c.gamw .= 0

    # Blade Element Values
    c.beta1 .= 0
    c.Cmag_rotor .= 0
    c.Re .= 0
    c.Ma .= 0
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
