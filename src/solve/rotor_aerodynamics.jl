"""
    calculate_induced_velocities_on_rotors(blade_elements, Gamr, vx_rb, vr_rb, gamb,
        vx_rw, vr_rw, gamw)

Calculate axial, radial, and tangential induced velocity on the rotors.

`Gamr::Matrix{Float}` : rotor circulation values (must be a matrix, even if only 1 rotor is being used). Size is number of blade elements by number of rotors.
`vx_rw::Matrix{Matrix{Float}}` : unit axial induced velocities of wakes on rotors. Must be a matrix of size number of rotors by number of wakes made up of matrices of size number of blade elements by number of wake "panels."
`gamw::Matrix{Float}` : wake circulation values (must be a matrix). Size is number of wake sheets by number of "panels" per sheet.
"""
function calculate_induced_velocities_on_rotors(
    blade_elements,
    Gamr,
    vx_rw,
    vr_rw,
    gamw,
    vx_rr,
    vr_rr,
    sigr,
    vx_rb=nothing,
    vr_rb=nothing,
    gamb=nothing,
    vx_rbte=nothing,
    vr_rbte=nothing,
    TEidxs=nothing;
    debug=false,
)

    # problem dimensions
    _, nrotor = size(Gamr) # number of rotors

    # initialize outputs
    vx_rotor = similar(Gamr) .= 0 # axial induced velocity
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
        if gamb != nothing
            @views vx_rotor[:, irotor] .+= vx_rb[irotor] * gamb
            @views vr_rotor[:, irotor] .+= vr_rb[irotor] * gamb

            if debug
                @views vxb_rotor[:, irotor] .+= vx_rb[irotor] * gamb
                @views vrb_rotor[:, irotor] .+= vr_rb[irotor] * gamb
            end

            #note: v?_rbte[irotor] should have been defined as zeros for endpoints that are not coincident.
            @views vx_rotor .+= vx_rbte[irotor] * gamb[TEidxs]
            @views vr_rotor .+= vr_rbte[irotor] * gamb[TEidxs]

            if debug
                @views vxb_rotor .+= vx_rbte[irotor] * gamb[TEidxs]
                @views vrb_rotor .+= vr_rbte[irotor] * gamb[TEidxs]
            end
        end

        # add wake induced velocities
        @views vx_rotor[:, irotor] .+= vx_rw[irotor] * gamw
        @views vr_rotor[:, irotor] .+= vr_rw[irotor] * gamw

        if debug
            @views vxw_rotor[:, irotor] .+= vx_rw[irotor] * gamw[:]
            @views vrw_rotor[:, irotor] .+= vr_rw[irotor] * gamw[:]
        end

        # add rotor induced velocities
        for jrotor in 1:nrotor
            @views vx_rotor[:, irotor] .+= vx_rr[irotor, jrotor] * sigr[:, jrotor]
            @views vr_rotor[:, irotor] .+= vr_rr[irotor, jrotor] * sigr[:, jrotor]

            if debug
                @views vxr_rotor[:, irotor] .+= vx_rr[irotor, jrotor] * sigr[:, jrotor]
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

    # return raw induced velocities
    if debug
        return vx_rotor,
        vr_rotor, vtheta_rotor, vxb_rotor, vrb_rotor, vxw_rotor, vrw_rotor, vxr_rotor,
        vrr_rotor
    else
        return vx_rotor, vr_rotor, vtheta_rotor
    end
end

"""
"""
function reframe_rotor_velocities(
    vx_rotor, vr_rotor, vtheta_rotor, Vinf, Omega, rotor_panel_centers
)

    # the axial component also includes the freestream velocity ( see eqn 1.87 in dissertation)
    Wx_rotor = vx_rotor .+ Vinf

    # the tangential also includes the negative of the rotation rate (see eqn 1.87 in dissertation)
    Wtheta_rotor = similar(vtheta_rotor) .= vtheta_rotor

    for i in 1:length(Omega)
        Wtheta_rotor[:, i] .-= Omega[i] .* rotor_panel_centers[:, i]
    end

    # meridional component
    Wm_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2)

    # Get the inflow magnitude at the rotor as the combination of all the components
    Wmag_rotor = sqrt.(Wx_rotor .^ 2 .+ vr_rotor .^ 2 .+ Wtheta_rotor .^ 2)

    return Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor
end

"""
"""
function calculate_rotor_velocities(Gamr, gamw, sigr, gamb, inputs)

    # - Get induced velocities on rotor planes - #
    vx_rotor, vr_rotor, vtheta_rotor = calculate_induced_velocities_on_rotors(
        inputs.blade_elements,
        Gamr,
        inputs.vx_rw,
        inputs.vr_rw,
        gamw,
        inputs.vx_rr,
        inputs.vr_rr,
        sigr,
        inputs.vx_rb,
        inputs.vr_rb,
        gamb,
        inputs.vx_rbte,
        inputs.vr_rbte,
        (p -> p.idx).(inputs.body_doublet_panels.TEnodes),
    )

    # - Reframe rotor velocities into blade element frames
    Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor = reframe_rotor_velocities(
        vx_rotor,
        vr_rotor,
        vtheta_rotor,
        inputs.Vinf,
        inputs.blade_elements.Omega,
        inputs.rotor_panel_centers,
    )

    return vx_rotor, vr_rotor, vtheta_rotor, Wx_rotor, Wtheta_rotor, Wm_rotor, Wmag_rotor
end


"""
"""
function gamma_from_coeffs!(Gamr, Wmag_rotor, blade_elements, cl)

    for (irotor, be) in enumerate(blade_elements)
    # calculate vortex strength
    Gamr[:,irotor] = @. 1.0 / 2.0 * Wmag_rotor[:,irotor] * be.chords * cl[:,irotor]
end

    return Gamr
end

# function sigma_from_coeffs!(sigr, Wmag_rotor, blade_elements, rho)
function sigma_from_coeffs!(sigr, Wmag_rotor, blade_elements,cd, rho)

    # calculate source strength
    # TODO: need to figure out if there should be a division by r here or not.
    #NOTE: theory sigr doesn't make any sense. consider using actual drag equation until you have good reason not to
    #= it doesn't look like DFDC has the division, even though it's in the theory doc
    #it's possible (but I haven't been able to find it, that they divide by r later for some reason
    =#
    # this one is in the theory doc
    # things converge twice as fast with this one...
    # sigr[:] .= B / (4.0 * pi * r) * Wmag_rotor * c * cd
    # this one is in the DFDC source code
    # things converge slower, likely because the values are so small...
    # sigr[:] .= B / (4.0 * pi) * Wmag_rotor * c * cd

    for (irotor,be) in enumerate(blade_elements)

        B = be.B

        for inode in 1:size(sigr,1)
            if inode ==1
                W=Wmag_rotor[inode, irotor]
                c = be.chords[inode]
                # cd = be.cd[inode]
                d = cd[inode]

            elseif inode==size(sigr,1)
                W=Wmag_rotor[inode-1, irotor]
                c = be.chords[inode-1]
                # cd = be.cd[inode-1]
                d = cd[inode-1]

            else
                W = Wmag_rotor[inode-1:inode, irotor]
                c = be.chords[inode-1:inode]
                # cd = be.cd[inode-1:inode]
                d = cd[inode-1:inode]
            end

            sigr[inode,irotor] = B / (4 * pi) * sum(W .* c .* d) / length(d)
            # sigr[:] .= B / (4 * pi) * sum(W .* c .* cd ./ r) / length(cd)
            # sigr[:] .= B / (4 * pi) * rho * sum(W^2 .* c .* cd ./ r) / length(cd)
        end
    end

    return sigr
end

function calc_reynolds(chord, Wmag_rotor, rho, mu)
    return chord .* abs.(Wmag_rotor) * rho / mu
end

"""
    calculate_gamma_sigma(blade_elements, Vm, vtheta_rotor)

Calculate rotor circulation and source strengths using blade element data and inflow velocities.

**Arguments:**
`Vm::Matrix{Float}` : Meridional velocity (including freestream and induced velocity) at rotor planes. Must be a matrix of size number of blade elements by number of rotors.

**Returns:**
`Gamr::Matrix{Float}` : Rotor circulations [num blade_elements x num rotors]
`sigr::Matrix{Float}` : Rotor panel source strengths [num blade_elements x num rotors]
"""
function calculate_gamma_sigma(
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
    verbose=false,
)

    # get floating point type
    TF = promote_type(
        eltype(blade_elements[1].chords),
        eltype(blade_elements[1].twists),
        eltype(Wm_rotor),
        eltype(Wtheta_rotor),
    )

    # get problem dimensions (number of radial stations x number of rotors)
    nr, nrotor = size(Wm_rotor)

    # initialize outputs
    Gamr = zeros(TF, nr, nrotor)
    sigr = zeros(TF, nr + 1, nrotor)

    # call in-place function
    return calculate_gamma_sigma!(
        Gamr,
        sigr,
        blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        debug=debug,
        verbose=verbose,
    )
end

"""
    calculate_gamma_sigma!(Gamr, sigr, blade_elements, Vm, vtheta_rotor)

In-place version of [`calculate_gamma_sigma`](@ref)

Note that circulations and source strengths must be matrices of size number of blade elements by number of rotors. (same dimensions as Vm and vtheta_rotor)
"""
function calculate_gamma_sigma!(
    Gamr,
    sigr,
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
    verbose=false
    )

        #update coeffs
        if debug
        cl, cd, phidb, alphadb, cldb, cddb=update_coeffs!(
        blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        debug=debug,
        verbose=verbose
    )
else
    cl, cd  = update_coeffs!(
        blade_elements,
        Wm_rotor,
        Wtheta_rotor,
        Wmag_rotor,
        freestream;
        debug=debug,
        verbose=verbose
    )
end

    # - get circulation strength - #
    gamma_from_coeffs!(
        Gamr,
        Wmag_rotor,
        blade_elements,
        cl
    )

    # - get source strength - #
    sigma_from_coeffs!(
        sigr,
        Wmag_rotor,
        blade_elements,
        cd,
        freestream.rhoinf,
    )

    return Gamr, sigr
end


function update_coeffs!(
    blade_elements,
    Wm_rotor,
    Wtheta_rotor,
    Wmag_rotor,
    freestream;
    debug=false,
    verbose=false,
)

    # problem dimensions
    nr, nrotor = size(Wmag_rotor) #num radial stations, num rotors

    cl = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
    cd = zeros(eltype(Wmag_rotor), size(Wmag_rotor))

    if debug
        # initialize extra outputs
        phidb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        alphadb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        cldb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
        cddb = zeros(eltype(Wmag_rotor), size(Wmag_rotor))
    end

    # loop through rotors
    for (irotor, be) in enumerate(blade_elements)

        # loop through radial stations
        for ir in 1:nr

            # extract blade element properties
            B = be.B # number of blades
            omega = be.Omega # rotation rate
            fliplift = be.fliplift

            c = be.chords[ir] # chord length
            twist = be.twists[ir] # twist
            stagger = be.stagger[ir]
            r = be.rbe[ir] # radius
            solidity = be.solidity[ir]

            # calculate angle of attack
            phi, alpha = calculate_inflow_angles(
                Wm_rotor[ir, irotor], Wtheta_rotor[ir, irotor], stagger
            )

            # get local Reynolds number
            reynolds = calc_reynolds(
                c, Wmag_rotor[ir, irotor], freestream.rhoinf, freestream.muinf
            )
            # c * abs(Wmag_rotor[ir, irotor]) * freestream.rhoinf / freestream.muinf

            # get local Mach number
            mach = Wmag_rotor[ir, irotor] / freestream.asound

            # - Update cl and cd values - #
            # be.cl[ir], be.cd[ir] = lookup_clcd(
            cl[ir,irotor], cd[ir,irotor] = lookup_clcd(
                blade_elements[irotor].inner_airfoil[ir],
                blade_elements[irotor].outer_airfoil[ir],
                blade_elements[irotor].inner_fraction[ir],
                Wmag_rotor[ir, irotor],
                solidity,
                stagger,
                alpha,
                phi,
                reynolds,
                mach,
                freestream.asound;
                verbose=verbose,
                fliplift=fliplift,
            )


            if debug
                phidb[ir, irotor] = phi
                alphadb[ir, irotor] = alpha
                cldb[ir, irotor] = cl
                cddb[ir, irotor] = cd
            end
        end
    end

    if debug
        return cl, cd, phidb, alphadb, cldb, cddb
        # return phidb, alphadb, cldb, cddb
    else
        # return nothing
        return cl, cd
    end
end

function lookup_clcd(
    inner_airfoil,
    outer_airfoil,
    inner_fraction,
    Wmag,
    solidity,
    stagger,
    alpha,
    inflow,
    reynolds,
    mach,
    asound;
    verbose=false,
    fliplift=false,
)

    # look up lift and drag data for the nearest two input sections
    # TODO: this breaks rotor aero tests... need to update those.
    if typeof(inner_airfoil) <: DFDCairfoil
        # - DFDC Airfoil Parameter - #

        # get inner values
        clin, cdin, _ = dfdc_clcdcm(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            inner_airfoil,
            asound;
            verbose=verbose,
            fliplift=fliplift,
        )
        # get outer values
        clout, cdout, _ = dfdc_clcdcm(
            Wmag,
            reynolds,
            solidity,
            stagger,
            alpha,
            outer_airfoil,
            asound;
            verbose=verbose,
            fliplift=fliplift,
        )

    elseif typeof(inner_airfoil) <: DTCascade
        # - Cascade Lookups - #
        # get inner values
        clin, cdin = caseval(inner_airfoil, stagger, inflow, reynolds, mach, solidity)
        # get outer values
        clout, cdout = caseval(outer_airfoil, stagger, inflow, reynolds, mach, solidity)
    elseif typeof(inner_airfoil) <: ccb.AFType
        # - Airfoil Lookups - #
        # get inner values
        clin, cdin = search_polars(inner_airfoil, alpha)
        # get outer values
        clout, cdout = search_polars(outer_airfoil, alpha)
    else
        @error "No blade element datatype: $(typeof(inner_airfoil)) defined."
    end

    # interpolate inner and outer values
    cl = fm.linear([0.0; 1.0], [clin, clout], inner_fraction)
    cd = fm.linear([0.0; 1.0], [cdin, cdout], inner_fraction)

    return cl, cd
end

"""
"""
function calculate_inflow_angles(Wm_rotor, Wtheta_rotor, stagger)
    #inflow angle
    # phi = atan.(Wm_rotor, -Wtheta_rotor)
    phi = atan.(-Wtheta_rotor, Wm_rotor)

    #angle of attack
    alpha = phi - stagger
    # alpha = twist .- phi

    return phi, alpha
end

"""
    search_polars(airfoil, alpha, re=0.0, ma=0.0)

Look up lift and drag data for an airfoil using CCBlade
TODO: add in cascade database search at some point.
"""
search_polars(airfoil, alpha, re=0.0, ma=0.0) = ccb.afeval(airfoil, alpha, re, ma)

"""
DFDC-like polar function. copied from dfdc and adjusted for julia
TODO: add dfdc polar parameters object compatibility for airfoil field in blade elements objects
TODO: this is only for a single airfoil definition.  need to rememember to interpolate between sections if there are more than one (this happens near where this function is called, rather than in this function itself)
"""
function dfdc_clcdcm(
    inflow_magnitude,
    local_reynolds,
    local_solidity,
    local_stagger,
    alpha,
    afparams,
    asound;
    verbose=false,
    fliplift=false,
)

    #all these come from user defined inputs.
    (;
        alpha0,
        clmax,
        clmin,
        dclda,
        dclda_stall,
        dcl_stall,
        cdmin,
        clcdmin,
        dcdcl2,
        cmcon,
        Re_ref,
        Re_exp,
        mcrit,
    ) = afparams

    # factors for compressibility drag model, hhy 10/23/00
    # mcrit is set by user
    # effective mcrit is mcrit_eff = mcrit - clmfactor*(cl-clcdmin) - dmdd
    # dmdd is the delta mach to get cd=cdmdd (usually 0.0020)
    # compressible drag is cdc = cdmfactor*(mach-mcrit_eff)^mexp
    # cdmstall is the drag at which compressible stall begins

    cdmfactor = 10.0
    clmfactor = 0.25
    mexp = 3.0
    cdmdd = 0.0020
    cdmstall = 0.1000

    # prandtl-glauert compressibility factor
    msq = (inflow_magnitude / asound)^2
    msq_w = 2.0 * inflow_magnitude / asound^2
    if msq >= 1.0
        ma = sqrt(msq)
        if typeof(ma) != Float64
            maprint = ma.value
        else
            maprint = ma
        end
        if verbose
            @warn "clfactor: local mach number limited to 0.99, was $maprint"
        end
        msq = 0.99
        msq_w = 0.0
    end

    pgrt = 1.0 / sqrt(1.0 - msq)
    pgrt_w = 0.5 * msq_w * pgrt^3

    # mach number and dependence on velocity
    mach = sqrt(msq)
    mach_w = 0.0
    if mach != 0.0
        mach_w = 0.5 * msq_w / mach
    end

    # generate clfactor for cascade effects from section solidity
    clfactor = 1.0
    if local_solidity > 0.0
        clfactor = getclfactor(local_solidity, local_stagger)
    end

    # generate cl from dcl/dalpha and prandtl-glauert scaling
    cla = dclda * pgrt * (alpha - alpha0) * clfactor
    cla_alf = dclda * pgrt * clfactor
    cla_w = dclda * pgrt_w * (alpha - alpha0) * clfactor

    # effective clmax is limited by mach effects
    # reduces clmax to match the cl of onset of serious compressible drag
    clmx = clmax
    clmn = clmin
    dmstall = (cdmstall / cdmfactor)^(1.0 / mexp)
    clmaxm = max(0.0, (mcrit + dmstall - mach) / clmfactor) + clcdmin
    clmax = min(clmax, clmaxm)
    clminm = min(0.0, -(mcrit + dmstall - mach) / clmfactor) + clcdmin
    clmin = max(clmin, clminm)

    # cl limiter function (turns on after +-stall
    ecmax = exp(min(200.0, (cla - clmax) / dcl_stall))
    ecmin = exp(min(200.0, (clmin - cla) / dcl_stall))
    cllim = dcl_stall * log((1.0 + ecmax) / (1.0 + ecmin))
    cllim_cla = ecmax / (1.0 + ecmax) + ecmin / (1.0 + ecmin)

    # subtract off a (nearly unity) fraction of the limited cl function
    # this sets the dcl/dalpha in the stalled regions to 1-fstall of that
    # in the linear lift range
    fstall = dclda_stall / dclda
    clift = cla - (1.0 - fstall) * cllim
    cl_alf = cla_alf - (1.0 - fstall) * cllim_cla * cla_alf
    cl_w = cla_w - (1.0 - fstall) * cllim_cla * cla_w

    stallf = false
    if clift > clmax || clift < clmin
        stallf == true
    end

    # cm from cmcon and prandtl-glauert scaling
    cmom = pgrt * cmcon
    cm_al = 0.0
    cm_w = pgrt_w * cmcon

    # cd from profile drag, stall drag and compressibility drag
    # reynolds number scaling factor
    if (local_reynolds <= 0.0)
        rcorr = 1.0
        rcorr_rey = 0.0
    else
        rcorr = (local_reynolds / Re_ref)^Re_exp
        rcorr_rey = Re_exp / local_reynolds
    end

    # include quadratic lift drag terms from airfoil and annulus

    # cdcl2 = dcdcl2 + dcdcl2_stall
    cdcl2 = dcdcl2  # no chance of getting messed up...

    # in the basic linear lift range drag is a function of lift
    # cd = cd0 (constant) + quadratic with cl)
    cdrag = (cdmin + cdcl2 * (clift - clcdmin)^2) * rcorr
    cd_alf = (2.0 * cdcl2 * (clift - clcdmin) * cl_alf) * rcorr
    cd_w = (2.0 * cdcl2 * (clift - clcdmin) * cl_w) * rcorr
    cd_rey = cdrag * rcorr_rey

    # post-stall drag added
    fstall = dclda_stall / dclda
    dcdx = (1.0 - fstall) * cllim / (pgrt * dclda)
    dcd = 2.0 * dcdx^2
    dcd_alf = 4.0 * dcdx * (1.0 - fstall) * cllim_cla * cla_alf / (pgrt * dclda)
    dcd_w =
        4.0 *
        dcdx *
        ((1.0 - fstall) * cllim_cla * cla_w / (pgrt * dclda) - dcd / pgrt * pgrt_w)

    # compressibility drag (accounts for drag rise above mcrit with cl effects
    # cdc is a function of a scaling factor*(m-mcrit(cl))^mexp
    # dmdd is the mach difference corresponding to cd rise of cdmdd at mcrit
    dmdd = (cdmdd / cdmfactor)^(1.0 / mexp)
    critmach = mcrit - clmfactor * abs(clift - clcdmin) - dmdd
    critmach_alf = -clmfactor * abs(cl_alf)
    critmach_w = -clmfactor * abs(cl_w)
    if (mach < critmach)
        cdc = 0.0
        cdc_alf = 0.0
        cdc_w = 0.0
    else
        cdc = cdmfactor * (mach - critmach)^mexp
        cdc_w = mexp * mach_w * cdc / mach - mexp * critmach_w * cdc / critmach
        cdc_alf = -mexp * critmach_alf * cdc / critmach
    end

    fac = 1.0
    fac_w = 0.0
    # although test data does not show profile drag increases due to mach #
    # you could use something like this to add increase drag by prandtl-glauert
    # (or any function you choose)
    #   fac   = pg
    #    fac_w = pg_w
    # total drag terms
    cdrag = fac * cdrag + dcd + cdc
    cd_alf = fac * cd_alf + dcd_alf + cdc_alf
    cd_w = fac * cd_w + fac_w * cdrag + dcd_w + cdc_alf
    cd_rey = fac * cd_rey

    #jm: if flip lift is true, return negative of clift (for stators)
    return fliplift ? -clift : clift, cdrag, cmom
end

"""
dfdc function copied and adjusted for julia
calculates multi-plane cascade effects on lift slope as a function of solidity and stagger angle
solidty: b*c/(2*pi*r)
stagger angle is from axis (not plane of rotation), in radians
originally from a table-drive quadratic fit to a figure 6-29 in wallis, axial flow fans and ducts.
"""
function getclfactor(solidity, stagger)

    # data from table (these are factors for a quadratic fit
    x = [0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.1; 1.2; 1.3; 1.4; 1.5]
    a0 = [
        0.4755
        0.5255
        0.5722
        0.6142
        0.6647
        0.7016
        0.7643
        0.8302
        0.8932
        0.9366
        0.9814
    ]
    a1 = [
        -0.367495
        -0.341941
        -0.300058
        -0.255883
        -0.200593
        -0.114993
        -0.118602
        -0.130921
        -0.133442
        -0.077980
        -0.123071
    ]
    a2 = [
        0.489466
        0.477648
        0.453027
        0.430048
        0.381462
        0.310028
        0.298309
        0.285309
        0.263084
        0.184165
        0.251594
    ]

    sigi = 1.0 / solidity

    aa0 = quadspline(x, a0, sigi)
    aa1 = quadspline(x, a1, sigi)
    aa2 = quadspline(x, a2, sigi)

    # only valid for stagger 20deg to 90deg,
    # limit low stagger to 20deg value to give constant lift ratio below that
    dtr = pi / 180.0 #degrees to radians
    stagr = stagger
    if stagr < 20.0 * dtr
        stagr = 20.0 * dtr
    end
    if stagr > 90.0 * dtr
        stagr = 90.0 * dtr
    end

    # quadratic fit for clfactor at this sigma as function of stagger
    clfactor = aa0 + aa1 * stagr + aa2 * stagr * stagr

    # maximum value of lift ratio should be limited to 1.0
    clfactor = min(1.0, clfactor)

    return clfactor
end

# """
# """
# function quadspline(xdata, ydata, xpoint)
#     n = length(xdata)

#     if n == 1
#         return xdata[1]
#     end

#     ilow = 1
#     i = n

#     while (i - ilow > 1)
#         imid = round(Int, (i + ilow) / 2)
#         if (xpoint < xdata[imid])
#             i = imid
#         else
#             ilow = imid
#         end
#     end

#     ds = xdata[i] - xdata[i - 1]
#     t = (xpoint - xdata[i - 1]) / ds
#     ypoint = t * ydata[i] + (1.0 - t) * ydata[i - 1]
#     # xxs =  (ydata(i) - ydata(i-1))/ds

#     return ypoint
# end

######################################################################
#                                                                    #
#                           STUFF FOR WAKE                           #
#                                                                    #
######################################################################

"""
only used in post-process for cp.
expression not in dfdc theory, comes from source code.
todo: derive in theory doc
"""
function calculate_entropy_jumps(sigr, Vm)
    return sigr .* Vm
end

"""
    calculate_enthalpy_jumps(Gamr, Omega, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Gamr, Omega, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Gamr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            omega = Omega[irotor]
            B = num_blades[irotor]

            G = Gamr[ir, irotor]
            H_tilde[ir, irotor] += omega * B * G / (2.0 * pi)

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                H_tilde[ir, irotor] += H_tilde[ir, jrotor]
            end
        end
    end

    return H_tilde
end

"""
    calculate_net_circulation(Gamr, num_blades)

Calculate net circulation from upstream rotors.
"""
function calculate_net_circulation(Gamr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Gamr)

    # calculate net circulations
    Gam_tilde = similar(Gamr) .= 0
    for ir in 1:nr
        for irotor in 1:nrotor
            B = num_blades[irotor]

            G = Gamr[ir, irotor]
            Gam_tilde[ir, irotor] += B * G

            # add the upstream contributions
            for jrotor in 1:(irotor - 1)
                Gam_tilde[ir, irotor] += Gam_tilde[ir, jrotor]
            end
        end
    end

    return Gam_tilde
end