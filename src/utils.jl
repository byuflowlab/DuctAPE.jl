function ismonotonic(A, cmp=<)
    current = A[1]
    for i in 2:length(A)
        newval = A[i]
        cmp(newval, current) && return false
        current = newval
    end
    return true
end
