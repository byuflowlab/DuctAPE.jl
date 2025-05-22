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