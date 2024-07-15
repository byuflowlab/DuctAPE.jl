######################################################################
#                                                                    #
#                       Robustness Adjustments                       #
#                                                                    #
######################################################################

"""
    stall_limiters(
        aoa,
        cl,
        cd;
        clminid=nothing,
        clmaxid=nothing,
        cutoff_slope=0.1,
        N=20,
        blend_hardness=50
    )

Cuts off coefficient vs alpha curve at min and max coefficient and places rest of curve from -pi to min coeff and max coeff to pi according to user defined cutoff_slope (default 0.1)

# Arguments:
- `aoa::AbstractVector{Float}` : input angles of attack, in radians
- `cl::AbstractVector{Float}` : input lift coefficients
- `cd::AbstractVector{Float}` : input drag coefficients

# Keyword Arguments:
- `clminid::Float=nothing` : manually set index for minimum cl
- `clmaxid::Float=nothing` : manually set index for maximum cl
- `cutoff_slope::Float=0.1` : "post-stall" slope
- `blend_hardness::Float=50` : hardenss of blend between nominal polar and post-stall modifications.

# Returns:
- `aoa_ext::AbstractVector{Float}` : angles of attack for modified polar, in radians
- `cl_ext::AbstractVector{Float}` : modified lift coefficients
- `cd_ext::AbstractVector{Float}` : modified drag coefficients
"""
function stall_limiters(
    aoa, cl, cd; clminid=nothing, clmaxid=nothing, cutoff_slope=0.1, N=20, blend_hardness=50
)

    # find cl min, associated index, and angle of attack
    if isnothing(clminid)
        clmin, clminid = findmin(cl)
    else
        clmin = cl[clminid]
    end
    aoamin = aoa[clminid]
    cdmin = cd[clminid]

    # find cl max, associated index, and angle of attack
    if isnothing(clmaxid)
        clmax, clmaxid = findmax(cl)
    else
        clmax = cl[clmaxid]
    end
    aoamax = aoa[clmaxid]
    cdmax = cd[clmaxid]

    # get full span of aoa's
    aoaext = [
        range(-pi, aoamin, N)
        aoa[(clminid + 1):(clmaxid - 1)]
        range(aoamax, pi, N)
    ]

    #get function for positive stall region.
    cl_ps = @. cutoff_slope * (aoaext - aoamax) + clmax
    cd_ps = @. cutoff_slope * (aoaext - aoamax) + cdmax

    #get function for negative stall region flip slope sign if applied to drag curve
    cl_ns = @. cutoff_slope * (aoaext - aoamin) + clmin
    cd_ns = @. -1.0 * cutoff_slope * (aoaext - aoamin) + cdmin

    # fill nominal cls
    fillcl = [
        # cl_ns[1:(N - 1)]
        cl_ns[1:N]
        cl[(clminid + 1):(clmaxid - 1)]
        cl_ps[(end - N + 1):end]
        # cl_ps[(end - N + 2):end]
    ]
    # fill nominal cds
    fillcd = [
        # cd_ns[1:(N - 1)]
        cd_ns[1:N]
        cd[(clminid + 1):(clmaxid - 1)]
        cd_ps[(end - N + 1):end]
        # cd_ps[(end - N + 2):end]
    ]

    clblend1 = FLOWMath.sigmoid_blend.(cl_ns, fillcl, aoaext, aoamin, blend_hardness)
    clblend2 = FLOWMath.sigmoid_blend.(clblend1, cl_ps, aoaext, aoamax, blend_hardness)

    cdblend1 = FLOWMath.sigmoid_blend.(cd_ns, fillcd, aoaext, aoamin, blend_hardness)
    cdblend2 = FLOWMath.sigmoid_blend.(cdblend1, cd_ps, aoaext, aoamax, blend_hardness)

    aoa_ret = collect(range(-pi, pi, 361))
    cl_ret = FLOWMath.akima(aoaext, clblend2, range(-pi, pi, 361))
    cd_ret = FLOWMath.akima(aoaext, cdblend2, range(-pi, pi, 361))

    nanid = findfirst(x -> !isnan(x), cl_ret)
    if nanid > 1
        cl_ret[1:nanid] .= cl_ret[nanid]
    end

    nanid = findfirst(x -> !isnan(x), cd_ret)
    if nanid > 1
        cd_ret[1:nanid] .= cd_ret[nanid]
    end

    nanid = findlast(x -> !isnan(x), cl_ret)
    if nanid < length(cl_ret)
        cl_ret[nanid:end] .= cl_ret[nanid]
    end

    nanid = findlast(x -> !isnan(x), cd_ret)
    if nanid < length(cd_ret)
        cd_ret[nanid:end] .= cd_ret[nanid]
    end

    return aoa_ret, cl_ret, cd_ret
end

######################################################################
#                                                                    #
#                 Solidity and Stagger Corrections                   #
#                                                                    #
######################################################################

"""
    solidity_and_stagger_factor_smooth(solidity, stagger; blend_hardness=100)

A smoothed version of `solidity_and_stagger_factor`.

# Arguments:
- `solidity::Float` : local solidity
- `stagger::Float` : local stagger (in radians)

# Keyword Arguments:
- `blend_hardness::Float=100` : hardness for smoothing blends
"""
function solidity_and_stagger_factor_smooth(solidity, stagger; blend_hardness=100)

    # data from table (these are factors for a quadratic fit)
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

    # spline the data
    aa0 = quadspline(x, a0, sigi)
    aa1 = quadspline(x, a1, sigi)
    aa2 = quadspline(x, a2, sigi)

    # quadratic fit for clfactor at this sigma as function of stagger
    ssfun(s) = aa0 + aa1 * s + aa2 * s^2

    # if stagger <= pi / 4.0
    stagblendlow = FLOWMath.sigmoid_blend(
        ssfun(pi / 9.0), ssfun(stagger), stagger, pi / 9.0, blend_hardness
    )
    # else
    stagblendhi = FLOWMath.sigmoid_blend(
        stagblendlow, ssfun(pi / 2.0), stagger, pi / 2.0, blend_hardness
    )
    # end

    ssfactor =
        stagblendhi -
        FLOWMath.sigmoid_blend.(0.0, stagblendhi - 1.0, stagblendhi, 0.99, blend_hardness)

    return ssfactor
end

"""
    solidity_and_stagger_factor(solidity, stagger; blend_hardness=100)

Correction for airfoil data used in a high-solidity cascade application.
Correction is used in DFDC airfoils nominally and come from quadratic fits to curves in fig 6-29 "Axial Flow Fans and Ducts" by Wallis (1983).
Note that the corrections are really only meant for Wallis' custom airfoil design and specific conditions mentioned in the book.

# Arguments:
- `solidity::Float` : local solidity
- `stagger::Float` : local stagger (in radians)

# Keyword Arguments:
- `blend_hardness::Float=100` : hardness for smoothing blends
"""
function solidity_and_stagger_factor(solidity, stagger;)
    # data from table (these are factors for a quadratic fit)
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

    # spline the data
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

"""
    quadspline(xdata, ydata, xpoint)

Sample data in quadratic spline at give point.
"""
function quadspline(xdata, ydata, xpoint)
    n = length(xdata)

    if n == 1
        return xdata[1]
    end

    ilow = 1
    i = n

    while (i - ilow > 1)
        imid = round(Int, (i + ilow) / 2)
        if (xpoint < xdata[imid])
            i = imid
        else
            ilow = imid
        end
    end

    ds = xdata[i] - xdata[i - 1]
    t = (xpoint - xdata[i - 1]) / ds
    ypoint = t * ydata[i] + (1.0 - t) * ydata[i - 1]

    return ypoint
end

"""
    solidity_and_stagger(cl, solidity, stagger; blend_hardness=100)

Apply smoothed Wallis' cascade correction (see `solidity_and_stagger_factor_smooth`) to local lift.

# Arguments:
- `cl::AbstractVector{Float}` : input lift coefficients
- `solidity::Float` : local solidity
- `stagger::Float` : local stagger (in radians)

# Keyword Arguments:
- `blend_hardness::Float=100` : hardness of smoothing blends

# Returns:
- `cl_corr::AbstractVector{Float}` : corrected lift coefficients.
"""
function solidity_and_stagger(cl, solidity, stagger; blend_hardness=100)
    return cl .* solidity_and_stagger_factor_smooth(
        solidity, stagger; blend_hardness=blend_hardness
    )
end

"""
    solidity_and_stagger!(cl, solidity, stagger; blend_hardness=100)

In-place version of `solidity_and_stagger`.
"""
function solidity_and_stagger!(cl, solidity, stagger; blend_hardness=100)
    cl .*= solidity_and_stagger_factor_smooth(
        solidity, stagger; blend_hardness=blend_hardness
    )
    return cl
end

######################################################################
#                                                                    #
#                    Compressibility Corrections                     #
#                                                                    #
######################################################################

#---------------------------------#
#         Prandtl-Glauert         #
#---------------------------------#

"""
    prandtl_glauert_factor(mach; verbose=false, blend_range=0.02)

Smoothed Prandtl-Glauert Mach correction factor

# Arguments:
- `mach::Float` : Mach number

# Keyword Arguments:
- `blend_range::Float=0.02` : range for blending factor and max cutoff (allowing Mach >= 1.0 for continuity)
"""
function prandtl_glauert_factor(mach; blend_range=0.02)
    pgf(M) = 1.0 / sqrt(1.0 - min(M, 0.999)^2)

    return FLOWMath.quintic_blend(pgf(mach), pgf(0.99), mach, 0.975, blend_range)
end

# """
# generates Akima spline of Prandtl-Glauert correction factors from Mach = 0.0 to Mach = 0.99
# """
# function genpgspline(; N=100)
#     return FLOWMath.Akima(
#         range(0.0, 0.99; length=100), prandtl_glauert_factor.(range(0.0, 0.99; length=100))
#     )
# end

"""
    prandtl_glauert(cl, ma)

Applies Prandtl-Glauert correction

# Arguments:
- `cl::Float` : local lift coefficient

# Returns
- `cl_corr::AbstractVector{Float}` : corrected lift coefficients
"""
function prandtl_glauert(cl, ma)
    return cl .* prandtl_glauert_factor(ma)
end

"""
    prandtl_glauert!(cl, ma)

In place version of `pradtl_glauert`.
"""
function prandtl_glauert!(cl, ma)
    cl .*= prandtl_glauert_factor(ma; verbose=verbose)
    return cl
end

#---------------------------------#
#        Reynolds Effects         #
#---------------------------------#

"""
    re_drag(cd, re, re_ref; re_exp=0.5)

# Arguments:
- `cd::AbstractVector{Float}` : input drag coefficients
- `re::Float` : Current Reynolds number
- `re_ref::Float` : Reference Reynolds number (at which the cd's were generated)

# Keyword Arguments:
- `re_exp::Float=0.5` : should be 0.2 for laminar and 0.5 for turbulent flow

# Returns:
- `cd_corr::AbstractVector{Float}` : Reynolds corrected drag coefficients
"""
function re_drag(cd, re, re_ref; re_exp=0.5)
    return cd .* (re_ref / re)^re_exp
end

"""
    re_drag!(cd, re, re_ref; re_exp=0.5)

In-place version of `re_drag`.
"""
function re_drag!(cd, re, re_ref; re_exp=0.5)
    cd .*= (re_ref / re)^re_exp
    return cd
end

#---------------------------------#
#        Transonic Effects        #
#---------------------------------#

"""
    transonic_lift_limiter_smooth!(
        cl,
        mach,
        clcdmin,
        clmax,
        clmin,
        dclda;
        mcrit=0.7,
        dcl_stall=0.1,
        dclda_stall=0.1,
        cdmfactor=10.0,
        clmfactor=0.25,
        mexp=3.0,
        cdmstall=0.1000,
        blend_hardness=50,
    )

Smoothed, vectorized, in-place version of `transonic_lift_limiter`.

# Different Arguments:
- `cl::AbstractVector{Float}` : vector of lift coefficients

# Additional Keyword Argument:
- `blend_hardness::Float=50` : hardenss of smoothing blends
"""
function transonic_lift_limiter_smooth!(
    cl,
    mach,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1000,
    # misc
    blend_hardness=50,
)
    if dclda <= 0.0
        @error "lift curve slope MUST be positive"
    end

    dmstall = (cdmstall / cdmfactor)^(1.0 / mexp)

    # get new cl_max
    clmaxm = (mcrit + dmstall - mach) / clmfactor
    clmaxmsmooth =
        clmaxm - FLOWMath.sigmoid_blend(clmaxm, 0.0, clmaxm, 0.05, blend_hardness) + clcdmin
    clmaxsmooth =
        clmaxmsmooth - FLOWMath.sigmoid_blend(
            clmaxmsmooth - clmax, 0.0, clmax, clmaxmsmooth, blend_hardness
        )

    # get new cl_min
    clminm = -(mcrit + dmstall - mach) / clmfactor
    clminmsmooth =
        clminm - FLOWMath.sigmoid_blend(0.0, clminm, clminm, -0.05, blend_hardness) +
        clcdmin
    clminsmooth =
        clminmsmooth - FLOWMath.sigmoid_blend(
            0.0, clminmsmooth - clmin, clmin, clminmsmooth, blend_hardness
        )
    # println("clminsmooth: ", clminsmooth)
    #     return clminsmooth

    # calculate exponents for convenience
    cmax = (cl .- clmaxsmooth) / dcl_stall
    cmaxsmooth = @. cmax -
        FLOWMath.sigmoid_blend(0.0, cmax - 100.0, cmax, 100.0, blend_hardness)
    ecmaxsmooth = exp.(cmaxsmooth)

    cmin = (clminsmooth .- cl) / dcl_stall
    cminsmooth = @. cmin -
        FLOWMath.sigmoid_blend(0.0, cmin - 100.0, cmin, 100.0, blend_hardness)
    ecminsmooth = exp.(cminsmooth)

    # get limiter factor
    cllim = @. dcl_stall * log((1.0 + ecmaxsmooth) / (1.0 + ecminsmooth))

    # subtract off a (nearly unity) fraction of the limited cl function
    # this sets the dcl/dalpha in the stalled regions to 1-fstall of that
    # in the linear lift range
    if dclda == 0.0
        fstall = 0.0
    else
        fstall = dclda_stall / dclda
    end

    # Final Lift
    # clift = cl - (1.0 - fstall) * cllim
    @. cl -= (1.0 - fstall) * cllim

    return cl
end

"""
    transonic_lift_limiter(
        cl,
        mach,
        clcdmin,
        clmax,
        clmin,
        dclda;
        mcrit=0.7,
        dcl_stall=0.1,
        dclda_stall=0.1,
        cdmfactor=10.0,
        clmfactor=0.25,
        mexp=3.0,
        cdmstall=0.1000,
    )

Airfoil polar corrections due to transonic effects as found in XROTOR and DFDC.
Note that this correction is done nominally in the DFDC airfoil evaluation.

# Arguments:
- `cl::Float` : input lift coefficient
- `mach::Float` : Mach number
- `clcdmin::Float` : lift coefficient at minimum drag coefficient.
- `clmax::Float` : maximum lift coefficient
- `clmin::Float` : minimum lift coefficient
- `dclda::Float` : lift-curve slope
- `mcrit::Float=0.7` : critical Mach number
- `dcl_stall::Float=0.1` : cl increment from initial to total stall
- `dclda_stall::Float=0.1` :  lift curve slope post-stall (1/radians)
- `cdmfactor::Float=10.0` : factor hard coded in XROTOR and DFDC
- `clmfactor::Float=0.25` : factor hard coded in XROTOR and DFDC
- `mexp::Float=3.0` : factor hard coded in XROTOR and DFDC
- `cdmstall::Float=0.1000` : factor hard coded in XROTOR and DFDC

# Returns:
- `cl_corr:Float` : corrected lift coefficient
"""
function transonic_lift_limiter(
    cl,
    mach,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1000,
)
    dmstall = (cdmstall / cdmfactor)^(1.0 / mexp)

    clmaxm = max(0.0, (mcrit + dmstall - mach) / clmfactor) + clcdmin
    clmax = min(clmax, clmaxm)
    clminm = min(0.0, -(mcrit + dmstall - mach) / clmfactor) + clcdmin
    clmin = max(clmin, clminm)

    # calculate exponents for convenience
    ecmax = exp(min(200.0, (cl - clmax) / dcl_stall))
    ecmin = exp(min(200.0, (clmin - cl) / dcl_stall))

    # get limiter factor
    cllim = dcl_stall * log((1.0 + ecmax) / (1.0 + ecmin))

    # subtract off a (nearly unity) fraction of the limited cl function
    # this sets the dcl/dalpha in the stalled regions to 1-fstall of that
    # in the linear lift range
    fstall = dclda_stall / dclda

    # Final Lift
    clift = cl - (1.0 - fstall) * cllim

    return clift
end

"""
    transonic_drag_addition!(
        cd,
        cl,
        clcdmin,
        mach;
        mcrit=0.7,
        cdmfactor=10.0,
        clmfactor=0.25,
        mexp=3.0,
        cdmdd=0.0020,
        cdmstall=0.1000,
        absdx=0.0625,
        blend_hardness=50,
    )

Smoothed, vecotrized, in-place version of `transonic_drag_addition`.

# Different Arguments:
- `cd::AbstractVector{Float}` : vector of drag coefficients
- `cl::AbstractVector{Float}` : vector of lift coefficients

# Additional Keyword Argument:
- `blend_hardness::Float=50` : hardenss of smoothing blends
"""
function transonic_drag_addition!(
    cd,
    cl,
    clcdmin,
    mach;
    mcrit=0.7,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmdd=0.0020,
    cdmstall=0.1000,
    # smooth absolute value dx
    absdx=0.0625,
    blend_hardness=50,
)

    # if mach < mcrit
    #     return cd
    # else
    dmdd = (cdmdd / cdmfactor)^(1.0 / mexp)
    critmach = @. mcrit - clmfactor * FLOWMath.abs_smooth((cl - clcdmin), absdx) - dmdd
    cdc = @. cdmfactor * (mach - critmach)^mexp

    cdcsmooth = @. FLOWMath.sigmoid_blend(0.0, cdc, mach, mcrit, blend_hardness)
    cd .+= cdcsmooth

    return cd
end

"""
    transonic_drag_addition(
        cd,
        cl,
        clcdmin,
        mach;
        mcrit=0.7,
        cdmfactor=10.0,
        clmfactor=0.25,
        mexp=3.0,
        cdmdd=0.0020,
        cdmstall=0.1000,
        absdx=0.0625,
    )

Drag augmentation due to transonic effects as found in XROTOR and DFDC.
Note this is nominally applied to DFDC airfoil evaluation.

# Arguments:
- `cd::Float` : input drag coefficient
- `cl::Float` : input lift coefficient
- `clcdmin::Float` : lift coefficient at minimum drag coefficient.
- `mach::Float` : Mach number

# Keyword Arguments
- `mcrit::Float=0.7` : critical Mach number
- `cdmfactor::Float=10.0` : factor hard coded in XROTOR and DFDC
- `clmfactor::Float=0.25` : factor hard coded in XROTOR and DFDC
- `mexp::Float=3.0` : factor hard coded in XROTOR and DFDC
- `cdmstall::Float=0.1000` : factor hard coded in XROTOR and DFDC
- `absdx::Float=0.0625` : smoothing factor for smooth absolute value function

# Returns:
- `cl_corr:Float` : corrected lift coefficient
"""
function transonic_drag_addition(
    cd,
    cl,
    clcdmin,
    mach;
    mcrit=0.7,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmdd=0.0020,
    cdmstall=0.1000,
    # smooth absolute value dx
    absdx=0.0625,
)
    if mach < mcrit
        return cd
    else
        dmdd = (cdmdd / cdmfactor)^(1.0 / mexp)
        critmach = @. mcrit - clmfactor * FLOWMath.abs_smooth((cl - clcdmin), absdx) - dmdd
        cdc = @. cdmfactor * (mach - critmach)^mexp

        return cd .+ cdc
    end
end

######################################################################
#                                                                    #
#                     On-the-fly Correction Functions                #
#                                                                    #
######################################################################

"""
    corrected_clcd(af::AlphaReAF, alpha, Re, Mach, solidity, stagger; kwargs...)

Evaluates and applies on-the-fly corrections for airfoil lift and drag.
On-the-fly airfoil polar corrections include solidity/stagger corrections, Prandtl-Glauert compressibility corrections, and transonic lift limits and drag additions.

    corrected_clcd!(cl, cd, af::AlphaReAF, Re, alpha, Mach, solidity, stagger; kwargs...)

Evaluates and applies on-the-fly corrections for airfoil lift and drag in place.

    corrected_clcd!(cl, cd, Mach, solidity, stagger; kwargs...)

Applies on-the-fly corrections for airfoil lift and drag in place.

    corrected_clcd!(cl, cd, af::AlphaAF, alpha, Re, Mach, solidity, stagger; kwargs...)

Evaluates and applies on-the-fly corrections, including Reynolds corrections, for airfoil lift and drag in place


    corrected_clcd(cas::InReStSoMaCAS, inflow, Re, Mach, solidity, stagger)

Evaluates cascade lift and drag.


# Arguments:
**Coefficients**
- `cl::Float` : local lift coefficient
- `cd::Float` : local drag coefficient

**Airfoil Object**
- `af::AlphaReAF` : airfoil object of CCBlade type dependent on angle of attack and Reynolds number
or
- `af::AlphaAF` : airfoil object of CCBlade type dependent on angle of attack only
or
- `cas::InReStSoMaCAS` : cascade object depentent on inflow angle, Reynolds number, stagger, solidity, and Mach number.

**Flow Angle**
- `alpha::Float` : angle of attack, radians.  Used with airfoil types
or
- `inflow::Float` : inflow angle, radians.  Used with cascade types

**Flow Conditions**
- `Re::Float` : Reynolds number
- `Mach::Float` : Mach number

**Geometry**
- `solidity::Float` : Local solidity
- `stagger::Float` : Stagger angle, radians

# Keyword Arguments:
- `mcrit::Float`=0.7 : Critical Mach number

**rotorzloc airfoil type parameters for post-stall behavior**
- `dcl_stall::Float`=0.1 : change in cl from incipient to total stall, used in transonic lift limiter correction
- `dclda_stall::Float`=0.1 : Post-stall lift curve slope

**Correction factors that were hard coded in rotorzloc and DFDC**
- `cdmfactor::Float`=10.0 :
- `clmfactor::Float`=0.25 :
- `mexp::Float`=3.0 :
- `cdmstall::Float`=0.1 :
- `cdmdd::Float`=0.0020 :

**Smoothing Paramters**
- `ssblend_hardness::Float`=100.0 : sigmoid blending hardness for solidity/stagger corrections
- `transblendhardness::Float`=75.0 : sigmoid blending hardness for transonic corrections
- `absdx::Float`=0.0625 : smooth absolute value Δα (radians) for transonic drag addition

**Miscellaneous**
- `verbose::Bool`=false : Boolean of whether to print warnings, etc.
"""
function corrected_clcd(
    af::AlphaReAF,
    alpha,
    Re,
    Mach,
    solidity,
    stagger,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1,
    cdmdd=0.0020,
    # smooth absolute value dx
    absdx=0.0625,
    # smoothing hardness for sigmoid blends
    ssblend_hardness=100,
    transblendhardness=50,
    # misc
    verbose=false,
)
    # - Get cl, cd - #
    cl, cd = afeval(af, alpha, Re, Mach) #note: Mach is not actually used since the type is AlphaReAF and not AlphaReMachAF

    # - Apply Corrections - #
    corrected_clcd!(
        cl,
        cd,
        Re,
        Mach,
        solidity,
        stagger,
        clcdmin,
        clmax,
        clmin,
        dclda;
        mcrit=0.7,
        # rotorzloc airfoil type parameters for post-stall behavior
        dcl_stall=0.1,
        dclda_stall=0.1,
        # factors hard coded in rotorzloc and dfdc
        cdmfactor=10.0,
        clmfactor=0.25,
        mexp=3.0,
        cdmstall=0.1,
        cdmdd=0.0020,
        # smooth absolute value dx
        absdx=0.0625,
        # smoothing hardness for sigmoid blends
        ssblend_hardness=100,
        transblendhardness=50,
        # misc
        verbose=false,
    )

    return cl, cd
end

function corrected_clcd!(
    cl,
    cd,
    af::AlphaReAF,
    alpha,
    Re,
    Mach,
    solidity,
    stagger,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1,
    cdmdd=0.0020,
    # smooth absolute value dx
    absdx=0.0625,
    # smoothing hardness for sigmoid blends
    ssblend_hardness=100,
    transblendhardness=50,
    # misc
    verbose=false,
)
    # - Get cl, cd - #
    cl, cd .= afeval(af, alpha, Re, Mach) #note: Mach is not actually used since the type is AlphaReAF and not AlphaReMachAF

    # - Apply Corrections - #
    corrected_clcd!(cl, cd, Mach, solidity, stagger)

    return cl, cd
end

function corrected_clcd!(
    cl,
    cd,
    Mach,
    solidity,
    stagger,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1,
    cdmdd=0.0020,
    # smooth absolute value dx
    absdx=0.0625,
    # smoothing hardness for sigmoid blends
    ssblend_hardness=100,
    transblendhardness=50,
    # misc
    verbose=false,
)

    # - Apply Solidity Correction to Lift - #
    solidityandstagger!(cl, solidity, stagger; blend_hardness=ssblend_hardness)

    # - Apply Prandtl-Glauert Correction to Lift - #
    prandtl_glauert!(cl, Mach; verbose=verbose)

    # - Apply Transonic Corrections - #
    # To Lift
    transonic_lift_limiter_smooth!(
        cl,
        Mach,
        clcdmin,
        clmax,
        clmin,
        dclda;
        mcrit=mcrit,
        # rotorzloc airfoil type parameters for post-stall behavior
        dcl_stall=dcl_stall,
        dclda_stall=dclda_stall,
        # factors hard coded in rotorzloc and dfdc
        cdmfactor=cdmfactor,
        clmfactor=clmfactor,
        mexp=mexp,
        cdmstall=cdmstall,
        # misc
        blend_hardness=transblendhardness,
        verbose=verbose,
    )

    # To Drag
    transonic_drag_addition!(
        cd,
        cl,
        clcdmin,
        Mach;
        mcrit=mcrit,
        # factors hard coded in rotorzloc and dfdc
        cdmfactor=cdmfactor,
        clmfactor=clmfactor,
        mexp=mexp,
        cdmdd=cdmdd,
        cdmstall=cdmstall,
        # smooth absolute value dx
        absdx=absdx,
        # sigmoid blending hardness
        blend_hardness=transblendhardness,
        verbose=verbose,
    )

    return cl, cd
end

function corrected_clcd!(
    cl,
    cd,
    af::AlphaAF,
    alpha,
    Re,
    Mach,
    solidity,
    stagger,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # Reynolds correction numbers
    re_ref=2e6,
    re_exp=0.2,
    # rotorzloc airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in rotorzloc and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1,
    cdmdd=0.0020,
    # smooth absolute value dx
    absdx=0.0625,
    # smoothing hardness for sigmoid blends
    ssblend_hardness=100,
    transblendhardness=50,
    # misc
    verbose=false,
)

    # - Get cl, cd - #
    cl, cd .= afeval(af, alpha, Re, Mach)

    # - Apply Solidity Correction to Lift - #
    solidityandstagger!(cl, solidity, stagger; blend_hardness=ssblend_hardness)

    # - Apply Prandtl-Glauert Correction to Lift - #
    prandtl_glauert!(cl, Mach; verbose=verbose)

    # - APPLY REYNOLDS CORRECTION TOO - #
    re_drag!(cd, Re, re_ref; re_exp=re_exp)

    # - Apply Transonic Corrections - #
    # To Lift
    transonic_lift_limiter_smooth!(
        cl,
        mach,
        clcdmin,
        clmax,
        clmin,
        dclda;
        mcrit=mcrit,
        # rotorzloc airfoil type parameters for post-stall behavior
        dcl_stall=dcl_stall,
        dclda_stall=dclda_stall,
        # factors hard coded in rotorzloc and dfdc
        cdmfactor=cdmfactor,
        clmfactor=clmfactor,
        mexp=mexp,
        cdmstall=cdmstall,
        # misc
        blend_hardness=transblendhardness,
        verbose=verbose,
    )

    # To Drag
    transonic_drag_addition!(
        cd,
        cl,
        clcdmin,
        Mach;
        mcrit=mcrit,
        # factors hard coded in rotorzloc and dfdc
        cdmfactor=cdmfactor,
        clmfactor=clmfactor,
        mexp=mexp,
        cdmdd=cdmdd,
        cdmstall=cdmstall,
        # smooth absolute value dx
        absdx=absdx,
        # sigmoid blending hardness
        blend_hardness=transblendhardness,
        verbose=verbose,
    )

    return cl, cd
end

function corrected_clcd(cas::InReStSoMaCAS, inflow, Re, Mach, solidity, stagger)
    return caseval(cas, inflow, Re, stagger, solidity, Mach)
end
