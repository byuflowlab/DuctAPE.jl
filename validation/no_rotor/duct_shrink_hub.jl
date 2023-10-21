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

hubscale = range(3.0, 0.1; length=5)

pcp = plot(;
    xlabel="x",
    ylabel=L"c_p",
    yflip=true,
    extra_kwargs=Dict(:subplot => Dict("ylabel style" => "{rotate=-90}")),
)
include(savepath * "duct-xvcp-$(npanduct-1)-panels.jl")

rrange = range(myred.r, darkred.r, length(ductrs))
grange = range(myred.g, darkred.g, length(ductrs))
brange = range(myred.b, darkred.b, length(ductrs))

for (i, hs) in enumerate(hubscale)
    shrink_hub = copy(repanel_hub)
    shrink_hub[:, 2] .*= hs
    f = open(dispath * "shrinkhub-coordinates-scale-$(hs).dat", "w")
    for (x, r) in zip(shrink_hub[:, 1], shrink_hub[:, 2])
        write(f, "$x $r\n")
    end
    close(f)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    # panels = dt.generate_panels([repanel];body=true)
    panels = dt.generate_panels([repanel_duct, shrink_hub])

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
    AICn, AICt = dt.vortex_aic_boundary_on_boundary(
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
        pcp,
        xcp[1:(npanduct - 1)],
        cp[1:(npanduct - 1)];
        label="",
        color=RGB(rrange[i], grange[i], brange[i]),
    )
end

plot!(pcp, ductxcp[:, 1], ductxcp[:, 2]; label="", color=1)
savefig(
    pcp,
    savepath *
    "shrinkhub-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.pdf",
)
savefig(
    pcp,
    dispath *
    "shrinkhub-pressure-comp-$(npanduct-1)-duct-panels-$(npanhub-1)-hub-panels.tikz",
)
