"""
# Fields:
- `alpha0::Float` : zero lift angle of attack
- `clmax::Float` : maximum cl
- `clmin::Float` : minimum cl
- `dclda::Float` : lift curve slope (1/radians)
- `dclda_stall::Float` :  lift curve slope post-stall (1/radians)
- `dcl_stall::Float` : cl increment from initial to total stall.
- `cdmin::Float` : minimum cd
- `cldmin::Float` : cl at cdmin
- `dcddcl2::Float` : quadratic curve factor for cl vs cd curve \$\\left(\\frac{d(c_d)}{d(c_l^2)}\\right)\$
- `cmcon::Float` : pitching moment constant (unused right now)
- `Re_ref::Float` : reference Reynolds number at which cd values apply
- `Re_exp::Float` : Reynolds number exponent scaling \$\\left( c_d = c_d(Re/Re_{ref})^{Re_{exp}}\\right)\$ should be 0.2 for fully laminar and 0.5 for fully turbulent
- `mcrit::Float` : critical Mach number
- `correct_for_mach::Bool` : flag to add Prandtl-Glauert correction
- `correct_for_cascade::Bool` : flag to add cascade corrections
- `correct_for_reynolds::Bool` : flag to add reynolds drag correction
- `correct_for_transonic::Bool` : flag to add drag correction above critical mach number
"""
@kwdef struct DFDCairfoil{TF,TB}
    alpha0::TF = 0.0
    clmax::TF = 1.5
    clmin::TF = -0.5
    dclda::TF = 2.0 * pi
    dclda_stall::TF = 0.1
    dcl_stall::TF = 0.1
    cdmin::TF = 0.01
    clcdmin::TF = 0.5
    dcddcl2::TF = 0.005
    cmcon::TF = 0.0
    Re_ref::TF = 1e6
    Re_exp::TF = 0.35
    mcrit::TF = 0.7
    correct_for_mach::TB = true
    correct_for_cascade::TB = true
    correct_for_reynolds::TB = true
    correct_for_transonic::TB = true
end

"""
    dfdceval(
        inflow_magnitude,
        local_reynolds,
        local_solidity,
        local_stagger,
        alpha,
        afparams,
        asound;
        verbose=false,
        is_stator=0,
    )

DFDC-like polar function.

# Arguments:
- `inflow_magnitude::Float` : Velocity magnitude of inflow
- `local_reynolds::Float` : local Reynolds number
- `local_solidity::Float` : local Solidity
- `local_stagger::Float` : local Stagger angle (radians)
- `alpha::Float` : local angle of attack (radians)
- `afparams::C4Blade.DFDCairfoil` : DFDCairfoil object
- `asound::Float` : speed of sound

# Keyword Arguments:
- `verbose::Bool=false::` : print verbose statements
- `is_stator::Int=0` : flag to flip lift values (e.g. for stators)

# Returns:
- `cl::Float` : lift coefficient corrected for compressibility, transonic regime, solidity and stagger as required.
- `cd::Float` : lift coefficient corrected for compressibility, transonic regime, solidity and stagger as required.
- `cmom::Float` : returned from input
"""
function dfdceval(
    inflow_magnitude,
    local_reynolds,
    local_solidity,
    local_stagger,
    alpha,
    afparams,
    asound;
    verbose=false,
    is_stator=0,
)
    if !iszero(is_stator)
        # if stator and you want to use a normal airfoil, flip lift will flip the angle of attack to effectively flip the airfoil camber direction.
        alpha *= -1.0
    end

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
        dcddcl2,
        cmcon,
        Re_ref,
        Re_exp,
        mcrit,
    ) = afparams

    (; correct_for_mach, correct_for_reynolds, correct_for_cascade, correct_for_transonic) =
        afparams

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

    mach = 0.0
    mach_w = 0.0
    pgrt = 1.0
    pgrt_w = 0.0
    if correct_for_mach
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
    end

    # generate clfactor for cascade effects from section solidity
    clfactor = 1.0
    if correct_for_cascade
        if local_solidity > 0.0
            clfactor = getclfactor(local_solidity, local_stagger)
        end
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
    if (local_reynolds <= 0.0) || !correct_for_reynolds
        rcorr = 1.0
        rcorr_rey = 0.0
    else
        rcorr = (local_reynolds / Re_ref)^Re_exp
        rcorr_rey = Re_exp / local_reynolds
    end

    # include quadratic lift drag terms from airfoil and annulus

    # cdcl2 = dcddcl2 + dcddcl2_stall
    cdcl2 = dcddcl2  # no chance of getting messed up...

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
    if (mach < critmach) || !correct_for_transonic
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

    #jm: if flip lift is true, return negative of clift (for stators) to put the force in the correct direction.
    # return !iszero(is_stator) ? -clift : clift, cdrag, cmom
    return clift, cdrag, cmom
end

"""
    getclfactor(solidity, stagger)

Calculates multi-plane cascade effects on lift slope as a function of solidity and stagger angle.

Comes from a table-driven quadratic fit to a figure 6-29 in Wallis, "Axial Flow Fans and Ducts."

# Arguments:
- `solidty::Float` : local solidity: \$\\frac{bc}{2*\\pi*r}\$
- `stagger::Float` : local stagger angle is from axis (not plane of rotation), in radians
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
