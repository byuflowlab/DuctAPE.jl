include("gen_logo.jl")

# axis of rotation
plot!(
    [nz[1], wg[1, end, 1]],
    -0.001 * ones(2);
    color=plotsgray,
    label="",
    lw=lw - 1,
)
plot!(nz[13] * ones(2), [-0.075, 0.075]; color=plotsgray, label="", lw=lw - 1)
plot!(
    nz[13] * ones(2) .+ 0.075,
    [-0.075, 0.075];
    color=plotsgray,
    label="",
    lw=lw - 1,
)
plot!(
    wg[1, end - 5, 1] * ones(2),
    [-0.075, 0.075];
    color=plotsgray,
    label="",
    lw=lw - 1,
)
plot!(
    wg[1, end - 5, 1] * ones(2) .- 0.075,
    [-0.075, 0.075];
    color=plotsgray,
    label="",
    lw=lw - 1,
)

##### ----- SAVE ----- #####
plot!(; grid=false, background_color=nothing)
savefig("logo-dark.svg")