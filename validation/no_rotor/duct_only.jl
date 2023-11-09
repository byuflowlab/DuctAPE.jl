#=
Isolated Duct Validation using data from Lewis annular airofils section.
Geometry is a NACA 66-015, with coordinates generated from OpenVSP, with the repeated LE point manually removed.
In addition, geometry was interpolated using a cosine spaced scheme
=#

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
# convenience path for saving directly to dissertation directory #TODO: remove before public release
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

# - load geometry - #
# read data file
include(project_dir * "/test/data/naca_662-015_smooth.jl")
# put coordinates together
coordinates = reverse(duct_coordinates; dims=1)

# Run a grid refinement loop
# Note that numerical integration has difficulty above 700 panels when using cosine spacing on this geometry
npans = [21, 31, 41, 51, 61, 71, 81, 91, 101, 151, 161, 201, 301, 401, 501, 601, 701]
# initialize refinement metric
cpsums = zeros(length(npans))

# loop through refinement levels
for (i, npan) in enumerate(npans)
    # print number of panels
    println("N Panels = ", npan - 1)

    # interpolate geometry with cosine spacing and given number of panels
    repanel = dt.repanel_airfoil(coordinates; N=npan, normalize=false)

    # save specific example refinements
    if npan == 161
        f = open(dispath * "duct-coordinates-$(npan-1)-panels.dat", "w")
    else
        f = open(savepath * "duct-coordinates-$(npan-1)-panels.dat", "w")
    end
    for (x, r) in zip(repanel[:, 1], repanel[:, 2])
        write(f, "$x $r\n")
    end
    close(f)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    panels = dt.generate_panels([repanel])

    # rename axial coordinates for convenience in later plotting
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


    # - Boundary on internal psuedo control point influence coefficients - #
    AICpcp, _ = dt.vortex_aic_boundary_on_field(
        panels.itcontrolpoint,
        panels.itnormal,
        panels.ittangent,
        panels.node,
        panels.nodemap,
        panels.influence_length,
    )

    ## -- Manually Assemble Linear System -- ##
    # need to manually assemble since we don't have a way of automating the isolated case right now
    # initialize LHS Matrix
    LHS = zeros(size(AICn)[2] + 1, size(AICn)[2] + 1)
    # fill LHS matrix with standard AIC terms
    LHS[1:size(AICn, 1), 1:size(AICn, 2)] .= AICn
    # add on pseudo control point influence terms
    LHS[size(AICn, 2), 1:size(AICn, 2)] .= AICpcp'
    # add in dummy variable influence terms
    LHS[1:size(AICn, 1), size(AICn, 2)+1] .= 1.0
    # add kutta condition
    LHS[size(AICn, 2) + 1, 1] = LHS[size(AICn, 2) + 1, size(AICn, 2)] = 1.0

    # - assemble RHS - #
    RHS = dt.freestream_influence_vector(panels.normal, Vsmat)
    push!(RHS, -1.0)
    push!(RHS, 0.0)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#
    # standard linear solve
    gamb = LHS \ RHS
    # don't include dummy variable associated with internal panel
    gamb = gamb[1:(end - 1)]

    # plot raw strengths
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

    # save outputs for comparision
    f = open(savepath * "duct-xvcp-$(npan-1)-panels.jl", "w")
    write(f, "ductxcp = [\n")
    for (x, c) in zip(panels.controlpoint[:, 1], cp)
        write(f, "$x $c\n")
    end
    write(f, "]")
    close(f)

    # save outputs for comparision
    f = open(savepath * "duct-xvvs-$(npan-1)-panels.jl", "w")
    write(f, "ductxvvs = [\n")
    for (x, c) in zip(panels.controlpoint[:, 1], Vtan)
        write(f, "$x $c\n")
    end
    write(f, "]")
    close(f)

    #---------------------------------#
    #             PLOTTING            #
    #---------------------------------#
    #plot outputs
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

    savefig(savepath * "duct-pressure-comp-$(npan-1)-panels.pdf")
    if npan == 161
        savefig(dispath * "duct-pressure-comp-$(npan-1)-panels.tikz")
    end

    # save refinement metric
    cpsums[i] = sum(cp .* panels.influence_length)
end

# look at relative error for example case and most refined case
id160 = findfirst(x -> x == 161, npans)

relerr = (cpsums[end] - cpsums[id160]) / cpsums[end] * 100

print("relative err from 160 to 800 panels: ", relerr, "%")

# - plot convergence - #
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
    [npans[id160]] .- 1, [cpsums[id160]]; seriestype=:scatter, markershape=:rect, label=""
)
annotate!(
    npans[id160] + 75,
    cpsums[id160] + 0.001,
    text("160 panels", 10, :left, :bottom; color=myred),
)

savefig(savepath * "duct-grid-refinement.pdf")
savefig(dispath * "duct-grid-refinement.tikz")
