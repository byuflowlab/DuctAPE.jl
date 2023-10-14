#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/no_rotor/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

# - load experimental data - #
include(project_dir * "/test/data/naca_662-015.jl")
# coordinates = lewis_duct_coordinates

# - load geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
coordinates = reverse(duct_coordinates; dims=1)

npans = [21, 31, 41, 51, 61, 71, 81, 91, 101, 151, 201, 301, 401, 501, 601, 701, 801]#, 901, 1001, 2001, 3001, 4001, 5001]
cpsums = zeros(length(npans))
for (i, npan) in enumerate(npans)
    println("N Panels = ", npan - 1)
    repanel = dt.repanel_airfoil(coordinates; N=npan, normalize=false)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    # panels = dt.generate_panels([repanel];body=true)
    panels = dt.generate_panels([repanel])

    # ##### ----- Visualize to Check ----- #####
    # visualize_paneling(;
    #     body_panels=panels,
    #     coordinates=[coordinates],
    #     controlpoints=true,
    #     nodes=true,
    #     normals=true,
    #     savepath=savepath,
    #     filename=["duct-geometry.pdf"],
    # )

    xn = panels.node[:, 1]
    xcp = panels.controlpoint[:, 1]

    #---------------------------------#
    #       Operating Conditions      #
    #---------------------------------#

    # Define freestream on panels
    Vinf = 1.0 #magnitude doesn't matter yet.
    Vs = Vinf * [1.0 0.0] # axisymmetric, so no radial component
    Vsmat = repeat(Vs, size(panels.controlpoint, 1)) # need velocity on each panel

    #---------------------------------#
    #        Induced Velocities       #
    #---------------------------------#

    # - Initial System Matrices - #
    AICn, AICt = dt.vortex_panel_influence_matrices(
        panels.controlpoint,
        panels.normal,
        panels.tangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    kids = [
        size(AICn)[1]+1 1
        size(AICn)[1]+1 size(AICn)[2]
    ]

    LHS = zeros(size(AICn)[1] + 1, size(AICn)[2])

    dt.add_kutta!(LHS, AICn, kids)

    RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
    push!(RHS, 0.0)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#
    gamb = LHS \ RHS

    pg = plot(xn, gamb; xlabel="x", ylabel="node strengths", label="")
    savefig(pg, savepath * "duct-gammas.pdf")

    #---------------------------------#
    #         Post-Processing         #
    #---------------------------------#

    ## --- Velocity Contributions --- ###

    # get tangent
    Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

    # add in body induced tangent velocity
    Vtan .+= AICt * gamb

    # add in jump term
    jump = (gamb[1:(end - 1)] + gamb[2:end]) / 2
    Vtan .-= jump / 2.0

    ### --- Steady Surface Pressure --- ###
    cp = 1.0 .- (Vtan / Vinf) .^ 2

    #---------------------------------#
    #             PLOTTING            #
    #---------------------------------#
    pp = plot(;
        xlabel="x",
        ylabel=L"c_p",
        yflip=true,
        extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
    )
    sp = pp[1]
    plot!(
        pp,
        pressurexupper,
        pressureupper;
        seriestype=:scatter,
        color=byublue,
        markershape=:utriangle,
        label="exp outer",
    )
    plot!(
        pp,
        pressurexlower,
        pressurelower;
        seriestype=:scatter,
        color=byublue,
        markershape=:dtriangle,
        label="exp inner",
    )
    plot!(pp, xcp, cp; label="DuctAPE", color=1)

    savefig(savepath * "duct-pressure-comp-$(npan)-panels.pdf")
    savefig(savepath * "duct-pressure-comp-$(npan)-panels.tikz")
    cpsums[i] = sum(cp .* panels.influence_length)
end

id150 = findfirst(x -> x == 151, npans)

relerr = (cpsums[end] - cpsums[id150]) / cpsums[end] * 100

print("relative err from 100 to 800 panels: ", relerr, "%")

pconv = plot(;
    xlabel="Number of Panels",
    xscale=:log10,
    ylabel=L"\sum_{i=1}^N \left[c_{p_i} \Delta s_i\right]",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
plot!(
    npans .- 1,
    cpsums;
    linestyle=:dot,
    linecolor=byublue,
    markercolor=myblue,
    markershape=:rect,
    markersize=3,
    label="",
)
plot!(
    [npans[id150]] .- 1, [cpsums[id150]]; seriestype=:scatter, markershape=:rect, label=""
)
annotate!(
    npans[id150] + 75,
    cpsums[id150] + 0.001,
    text("150 panels", 10, :left, :bottom; color=myred),
)

savefig(savepath * "duct-grid-refinement.pdf")
savefig(savepath * "duct-grid-refinement.tikz")

