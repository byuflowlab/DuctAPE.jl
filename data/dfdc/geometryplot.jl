include("dstestr2_case.jl")

using Plots
include("../../plots_default.jl")

plot(ductx, ductr; markershape=:circle, markersize=1.5)
plot!(hubx, hubr; markershape=:circle, markersize=1.5)

plot!(dragx, dragr; markershape=:circle, markersize=1.5)

plot!(xdisk1 * ones(length(r1)), r1; color=4, markershape=:circle, markersize=1.5)
chord1rot = chord1 .* sind.(beta1)
plot!(xdisk1 .- 0.25 * chord1rot, r1; color=4)
plot!(xdisk1 .+ 0.75 * chord1rot, r1; color=4)

plot!(xdisk2 * ones(length(r2)), r2; color=5, markershape=:circle, markersize=1.5)
chord2rot = chord2 .* sind.(beta2)
plot!(xdisk2 .- 0.25 * chord2rot, r2; color=5)
plot!(xdisk2 .+ 0.75 * chord2rot, r2; color=5)
