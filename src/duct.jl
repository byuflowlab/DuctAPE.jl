#=
Duct Setup Functions
=#

function split_wall(x,z)

    _, xminidx = findmin(x)

    return x[1:xminidx], x[xminidx:end], z[1:xminidx], z[xminidx:end]
end
