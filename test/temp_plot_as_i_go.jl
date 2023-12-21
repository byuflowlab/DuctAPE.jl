include("../visualize/plots_default.jl")
include("data/basic_two_rotor_for_test.jl")

pg = plot(; aspectratio=1, size=(1000, 1000))

plot!(
    pg,
    duct_coordinates[:, 1],
    duct_coordinates[:, 2];
    marker=true,
    markersize=2,
    color=myblue,
    label="",
)
plot!(
    pg,
    hub_coordinates[:, 1],
    hub_coordinates[:, 2];
    marker=true,
    markersize=2,
    color=myblue,
    label="",
)

for rsp in rotorstator_parameters
    plot!(
        pg,
        ones(length(rsp.r)) * rsp.rotorzloc,
        rsp.r .* rsp.Rtip;
        color=myred,
        marker=true,
        markersize=2,
        label="",
    )
end

savefig(pg, "geomplot1.pdf")

pg = plot(; aspectratio=1, size=(1000, 1000))

plot!(
    pg,
    rp_duct_coordinates[1, :],
    rp_duct_coordinates[2, :];
    markersize=3,
    markershape=:rect,
    color=myblue,
    label="",
)
plot!(
    pg,
    rp_hub_coordinates[1, :],
    rp_hub_coordinates[2, :];
    markersize=3,
    markershape=:rect,
    color=myblue,
    label="",
)

for rsp in rotorstator_parameters
    plot!(
        pg,
        ones(length(rsp.r)) * rsp.rotorzloc,
        rsp.r .* rsp.Rtip;
        color=myred,
        marker=true,
        markersize=2,
        label="",
    )
end

for (z, r) in zip(eachcol(grid[1, :, :]), eachcol(grid[2, :, :]))
    plot!(pg, z, r; color=mygray, marker=true, markersize=2)
end

savefig(pg, "geomplot.pdf")
