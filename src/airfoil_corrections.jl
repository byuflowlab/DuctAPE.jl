"""
Prandtl-Glauert Correction
"""
function prandtlglauert(cl, ma)
    return cl ./ (1.0 - ma^2)
end

"""
Prandtl-Glauert Correction

In place: overwrites `cl`.
"""
function prandtlglauert!(cl, ma)
    cl ./= (1.0 - ma^2)
    return cl
end

"""
"""
function transonicliftlimiter(
    cl,
    mach,
    clcdmin,
    clmax,
    clmin,
    dclda;
    mcrit=0.7,
    # xrotor airfoil type parameters for post-stall behavior
    dcl_stall=0.1,
    dclda_stall=0.1,
    # factors hard coded in xrotor and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmstall=0.1000,
)
    dmstall = (cdmstall / cdmfactor)^(1.0 / mexp)

    # get new cl_max
    clmaxm = max(0.0, (mcrit + dmstall - mach) / clmfactor) + clcdmin
    clmax = min(clmax, clmaxm)

    # get new cl_min
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
reexp should be 0.2 for laminar and 0.5 for turbulent
"""
function redrag(cd, re, reref; reexp=0.5)
    return cd .* (reref / re)^reexp
end

"""
reexp should be 0.2 for laminar and 0.5 for turbulent
"""
function redrag!(cd, re, reref; reexp=0.5)
    cd .*= (reref / re)^reexp
    return cd
end

"""
"""
function transonicdragaddition(
    cd,
    cl,
    clcdmin,
    mach;
    mcrit=0.7,
    # factors hard coded in xrotor and dfdc
    cdmfactor=10.0,
    clmfactor=0.25,
    mexp=3.0,
    cdmdd=0.0020,
    cdmstall=0.1000,
    # smooth absolute value dx
    absdx=0.0625,
)
    dmdd = (cdmdd / cdmfactor)^(1.0 / mexp)
    critmach = @. mcrit - clmfactor * FLOWMath.abs_smooth((cl - clcdmin), absdx) - dmdd
    cdc = @. cdmfactor * (mach - critmach)^mexp

    return cd .+ cdc
end

"""
cuts off coefficient vs alpha curve at min and max coefficient and places rest of curve from -pi to min coeff and max coeff to pi according to user defined cutoff_slope (default 0.1)
"""
function stalllimiters(aoa, coeff; drag=false, cutoff_slope=0.1, N=20, blend_hardness=50)

    # find cl min, associated index, and angle of attack
    coeffmin, coeffminid = drag ? (coff[1], 1) : findmin(coeff)
    aoamin = aoa[coeffminid]

    # find coeff max, associated index, and angle of attack
    coeffmax, coeffmaxid = drag ? (coeff[end], length(coeff)) : findmax(coeff)
    aoamax = aoa[coeffmaxid]

    # get full span of aoa's
    aoaext = [
        range(-pi, aoamin, N)
        aoa[(coeffminid + 1):(coeffmaxid - 1)]
        range(aoamax, pi, N)
    ]

    #get function for positive stall region.
    coeff_ps = @. cutoff_slope * (aoaext - aoamax) + coeffmax

    #get function for negative stall region flip slope sign if applied to drag curve
    coeff_ns = @. (drag ? -1.0 : 1.0) * cutoff_slope * (aoaext - aoamin) + coeffmin

    # fill nominal coeffs
    fillcoeff = [
        coeff_ns[1:(N - 1)]
        coeff[coeffminid:coeffmaxid]
        coeff_ps[(end - N + 2):end]
    ]

    blend1 = FLOWMath.sigmoid_blend.(coeff_ns, fillcoeff, aoaext, aoamin, blend_hardness)
    blend2 = FLOWMath.sigmoid_blend.(blend1, coeff_ps, aoaext, aoamax, blend_hardness)

    return collect(range(-pi, pi, 361)), FLOWMath.akima(aoaext, blend2, range(-pi, pi, 361))
end
