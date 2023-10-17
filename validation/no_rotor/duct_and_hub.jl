#---------------------------------#
#              SETUP              #
#---------------------------------#

# - Get Project Directory - #
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

# create save path
savepath = project_dir * "/validation/no_rotor/figs/"
dispath =
    project_dir * "/../../Writing/dissertation/src/ductsolvercontents/ductsolverfigures/"

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

# - load experimental data - #
include(project_dir * "/test/data/naca_662-015.jl")
# coordinates = lewis_duct_coordinates

# - load duct geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
duct_coordinates = reverse(duct_coordinates; dims=1)

# - load hub geometry - #
# read data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
cut = 2
hub_coordinates = [x_hub[1:(end - cut)] r_hub[1:(end - cut)]]

npansduct = [41, 51, 61, 71, 81, 91, 101, 161, 201, 301, 401, 501, 601, 701, 801]
npanshub = ceil.(Int, npansduct ./ 2)
cpsums = zeros(length(npansduct))
for (i, (npanduct, npanhub)) in enumerate(zip(npansduct, npanshub))
    println("N Duct Panels = ", npanduct - 1)
    println("N Hub Panels = ", npanhub - 1)

    repanel_duct = dt.repanel_airfoil(duct_coordinates; N=npanduct, normalize=false)
    repanel_hub = dt.repanel_revolution(hub_coordinates; N=npanhub, normalize=false)

    # f = open(savepath * "duct-coordinates-$(npanduct-1)-panels.dat", "w")
    # for (x, r) in zip(repanel_duct[:, 1], repanel_duct[:, 2])
    #     write(f, "$x $r\n")
    # end
    # close(f)

    # f = open(savepath * "hub-coordinates-$(npanhub-1)-panels.dat", "w")
    # for (x, r) in zip(repanel_hub[:, 1], repanel_hub[:, 2])
    #     write(f, "$x $r\n")
    # end
    # close(f)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    # panels = dt.generate_panels([repanel];body=true)
    panels = dt.generate_panels([repanel_duct, repanel_hub])

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
        size(AICn)[1]+1 npanduct
        size(AICn)[1]+2 npanduct+1
    ]

    LHS = zeros(size(AICn)[1] + 2, size(AICn)[2])

    dt.add_kutta!(LHS, AICn, kids)

    RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
    push!(RHS, 0.0)
    push!(RHS, 0.0)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#
    gamb = LHS \ RHS

    # pg = plot(xn, gamb; xlabel="x", ylabel="node strengths", label="")
    # savefig(pg, savepath * "duct-gammas.pdf")

    #---------------------------------#
    #         Post-Processing         #
    #---------------------------------#

    ## --- Velocity Contributions --- ###

    # get tangent
    Vtan = [dt.dot(v, t) for (v, t) in zip(eachrow(Vsmat), eachrow(panels.tangent))]

    # add in body induced tangent velocity
    Vtan .+= AICt * gamb

    # add in jump term
    jumpduct = (gamb[1:(npanduct - 1)] + gamb[2:npanduct]) / 2
    jumphub = (gamb[(npanduct + 1):(end - 1)] + gamb[(npanduct + 2):end]) / 2
    jump = [jumpduct; jumphub]
    Vtan .-= jump / 2.0

    ### --- Steady Surface Pressure --- ###
    cp = 1.0 .- (Vtan / Vinf) .^ 2

    #---------------------------------#
    #             PLOTTING            #
    #---------------------------------#
    pcp = plot(;
        xlabel="x",
        ylabel=L"c_p",
        yflip=true,
        extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
    )
    plot!(
        pcp,
        pressurexupper,
        pressureupper;
        seriestype=:scatter,
        color=byublue,
        markershape=:utriangle,
        label="Experimental Nacelle",
    )
    plot!(
        pcp,
        pressurexlower,
        pressurelower;
        seriestype=:scatter,
        color=byublue,
        markershape=:dtriangle,
        label="Experimental Casing",
    )
    include(savepath * "duct-xvcp-$(npanduct-1)-panels.jl")
    plot!(pcp, ductxcp[:, 1], ductxcp[:, 2]; label="DuctAPE Isolated Duct", color=1)
    plot!(
        pcp,
        xcp[1:(npanduct - 1)],
        cp[1:(npanduct - 1)];
        label="DuctAPE Duct with Center Body",
        color=2,
    )

    savefig(
        pcp,
        savepath *
        "system-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
    )
    if npanduct == 161
        savefig(
            pcp,
            dispath *
            "system-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
        )
    end

    pvs = plot(;
        xlabel="x",
        ylabel=L"\frac{V_s}{V_\infty}",
        extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
    )
    plot!(
        pvs,
        Vs_over_Vinf_x,
        Vs_over_Vinf_vs;
        seriestype=:scatter,
        color=byublue,
        markershape=:utriangle,
        label="Experimental Center Body",
    )
    include(savepath * "hub-xvvs-$(npanhub-1)-panels.jl")
    plot!(
        pvs,
        hubxvvs[:, 1],
        hubxvvs[:, 2];
        color=myblue,
        label="DuctAPE Isolated Center Body",
    )
    plot!(
        pvs,
        xcp[(npanduct):end],
        Vtan[(npanduct):end] ./ Vinf;
        color=myred,
        label="DuctAPE Center Body with Duct",
    )

    savefig(
        pvs,
        savepath *
        "system-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
    )
    if npanduct == 161
        savefig(
            pvs,
            dispath *
            "system-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
        )
    end
    cpsums[i] = sum(cp .* panels.influence_length)
end

id160 = findfirst(x -> x == 161, npansduct)

relerr = (cpsums[end] - cpsums[id160]) / cpsums[end] * 100

print("relative err: ", relerr, "%")

pconv = plot(;
    xlabel="Total Number of Panels",
    xscale=:log10,
    ylabel=L"\sum_{i=1}^N \left[c_{p_i} \Delta s_i\right]",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)

plot!(
    npansduct .+ npanshub .- 2,
    cpsums;
    linestyle=:dot,
    linecolor=byublue,
    markercolor=myblue,
    markershape=:rect,
    markersize=3,
    label="",
)

plot!(
    [npansduct[id160] + npanshub[id160]] .- 2,
    [cpsums[id160]];
    seriestype=:scatter,
    markershape=:rect,
    label="",
)

annotate!(
    npansduct[id160] + npanshub[id160] + 300,
    cpsums[id160] - 0.0005,
    text("240 total panels", 10, :left, :bottom; color=myred),
)

savefig(savepath * "duct-and-hub-grid-refinement.pdf")
savefig(dispath * "duct-and-hub-grid-refinement.tikz")
