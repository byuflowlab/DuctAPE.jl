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

npanduct = 161
npanhub = ceil.(Int, npanduct ./ 2)
repanel_duct = dt.repanel_airfoil(duct_coordinates; N=npanduct, normalize=false)
repanel_hub = dt.repanel_revolution(hub_coordinates; N=npanhub, normalize=false)
ductrs = collect(range(0.5, 2.0; length=5))

rrange = range(myred.r, darkred.r, length(ductrs))
grange = range(myred.g, darkred.g, length(ductrs))
brange = range(myred.b, darkred.b, length(ductrs))

pvs = plot(;
    xlabel="x",
    ylabel=L"\frac{V_s}{V_\infty}",
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
# plot!(
#     pvs,
#     Vs_over_Vinf_x,
#     Vs_over_Vinf_vs;
#     seriestype=:scatter,
#     color=byublue,
#     markershape=:utriangle,
#     label="Experimental Center Body",
# )
include(savepath * "hub-xvvs-$(npanhub-1)-panels.jl")
g = open(dispath * "tikzmovecolor.col", "w")

for (i, r) in enumerate(ductrs)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    # panels = dt.generate_panels([repanel];body=true)
    duct_coords = copy(repanel_duct)
    duct_coords[:, 2] .-= duct_coords[1,2]
    duct_coords[:, 2] .+= r
    panels = dt.generate_panels([duct_coords, repanel_hub])

    f = open(savepath * "moveduct-coordinates-r-$(r).dat", "w")
    for (x, r) in zip(duct_coords[:, 1], duct_coords[:, 2])
        write(f, "$x $r\n")
    end
    close(f)

    write(g, "color={rgb,1:red,$(rrange[i]); green,$(grange[i]); blue,$(brange[i])}\n")

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

    plot!(
        pvs,
        xcp[(npanduct):end],
        Vtan[(npanduct):end] ./ Vinf;
        color=RGB(rrange[i], grange[i], brange[i]),
        # label=i == 1 ? "Duct close and moving away" : "",
        label="",
    )
end
plot!(pvs, hubxvvs[:, 1], hubxvvs[:, 2]; color=myblue,label="")#, label="Without Duct")
close(g)

savefig(
    pvs,
    savepath *
    "ductmove-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
)
savefig(
    pvs,
    dispath *
    "ductmove-velocity-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
)
