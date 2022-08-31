function ismonotonic(A, cmp=<)
    current = A[1]
    for i in 2:length(A)
        newval = A[i]
        cmp(newval, current) && return false
        current = newval
    end
    return true
end

function lintran(rb1, rbend, ra1, raend, ra)
    return rb1 .+ (rbend - rb1) / (raend - ra1) .* (ra .- ra1)
end

function sinespace(n)
    theta = range(0, pi / 2; length=n)
    return sin.(theta)
end

function cosinespace(n)
    theta = range(0, pi / 2; length=n)
    return cos.(theta)
end

function get_omega(rpm)
    return rpm * pi / 30.0
end
