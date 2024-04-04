"""
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
    # xxs =  (ydata(i) - ydata(i-1))/ds

    return ypoint
end

"""
"""
function vectorize_airfoils(airfoils::AbstractArray{T}) where {T<:DFDCairfoil}

    # get type of airfoils
    TF = promote_type(eltype.(getproperty.(airfoils, :dclda))...)
    lenaf = length(airfoils)
    lenfields = length(fieldnames(T))

    afvec = zeros(TF, lenaf * lenfields)

    for (i, af) in enumerate(airfoils)
        afvec[(i - 1) * lenfields + 1] = af.alpha0
        afvec[(i - 1) * lenfields + 1 + 1] = af.clmax
        afvec[(i - 1) * lenfields + 1 + 2] = af.clmin
        afvec[(i - 1) * lenfields + 1 + 3] = af.dclda
        afvec[(i - 1) * lenfields + 1 + 4] = af.dclda_stall
        afvec[(i - 1) * lenfields + 1 + 5] = af.dcl_stall
        afvec[(i - 1) * lenfields + 1 + 6] = af.cdmin
        afvec[(i - 1) * lenfields + 1 + 7] = af.clcdmin
        afvec[(i - 1) * lenfields + 1 + 8] = af.dcddcl2
        afvec[(i - 1) * lenfields + 1 + 9] = af.cmcon
        afvec[(i - 1) * lenfields + 1 + 10] = af.Re_ref
        afvec[(i - 1) * lenfields + 1 + 11] = af.Re_exp
        afvec[(i - 1) * lenfields + 1 + 12] = af.mcrit
    end

    return afvec
end

function vectorize_airfoils(airfoils::AbstractArray{T}) where {T<:AFType}
    @error "no dispatch for vectorizing CCBlade airfoil types yet"

    return nothing
end

function vectorize_airfoils(airfoils::AbstractArray{T}) where {T<:DTCascade}
    @error "no dispatch for vectorizing C4Blade cascade types yet"

    return nothing
end
