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

# - load geometry - #
# read data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
cut = 2
coordinates = [x_hub[1:(end - cut)] r_hub[1:(end - cut)]]

npans = [21, 31, 41, 51, 61, 71, 81, 91, 101, 151, 201, 301, 401, 501]# 601, 701, 801]#, 901, 1001, 2001, 3001, 4001, 5001]
cpsums = zeros(length(npans))
for (i, npan) in enumerate(npans)
    println("N Panels = ", npan - 1)
    repanel = dt.repanel_revolution(coordinates; N=npan, normalize=false)

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
    #     filename=["hub-geometry.pdf"],
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
    ]

    LHS = zeros(size(AICn)[1] + 1, size(AICn)[2])

    dt.add_kutta!(LHS, AICn, kids)

    RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
    push!(RHS, 0.0)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#
    gamb = LHS \ RHS

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
    pg = plot(; xlabel="x", ylabel=L"\mathrm{Panel~strengths~}(\gamma^B)")
    plot!(pg, xn, gamb; label="")
    savefig(savepath * "hub-gammas.pdf")

    pp = plot(;
        xlabel="x",
        ylabel=L"\frac{V_s}{V_\infty}",
        extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
    )
    plot!(
        pp,
        Vs_over_Vinf_x,
        Vs_over_Vinf_vs;
        seriestype=:scatter,
        color=byublue,
        markershape=:utriangle,
        label="experimental",
    )
    plot!(pp, xcp, Vtan ./ Vinf; color=myblue, label="DuctAPE")

    savefig(savepath * "hub-velocity-comp-$(npan)-panels.pdf")
    savefig(savepath * "hub-velocity-comp-$(npan)-panels.tikz")
    cpsums[i] = sum(cp .* panels.influence_length)
end

id80 = findfirst(x -> x == 81, npans)

relerr = (cpsums[end] - cpsums[id80]) / cpsums[end] * 100

print("relative err: ", relerr, "%")

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
plot!([npans[id80]] .- 1, [cpsums[id80]]; seriestype=:scatter, markershape=:rect, label="")
annotate!(npans[id80] + 40, cpsums[id80], text("80 panels", 10, :left, :top; color=myred))

savefig(savepath * "hub-grid-refinement.pdf")
savefig(savepath * "hub-grid-refinement.tikz")

