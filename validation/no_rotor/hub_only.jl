#=
Isolated Hub Validation using data from Lewis bodies of revolution section.
Geometry is based on that provided by lewis, but manually generated based on the lead edge circle radius, overall lengths, etc.
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

# - load DuctAPE - #
using DuctAPE
const dt = DuctAPE

# - load plotting defaults - #
include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default.jl")

# - load geometry - #
# read data file
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# hub final r-coordinate needs to be set to zero so that it's not negative
r_hub[end] = 0.0
# put coordinates together
coordinates = [x_hub r_hub]

#TODO: check to see if all these cases still work
# use the same number of panels as the isolated duct case
npansduct = [41, 51, 61, 71, 81, 91, 101, 161, 201, 301, 401, 501, 601, 701]
# except divide by 2 since it's just a half
npans = ceil.(Int, npansduct ./ 2)
# initialize refinement metric
cpsums = zeros(length(npans))

# loop through refinement levels
for (i, npan) in enumerate(npans)
    # print number of panels
    println("N Panels = ", npan - 1)

    # interpolate geometry
    repanel = dt.repanel_revolution(coordinates; N=npan, normalize=false)

    #---------------------------------#
    #             Paneling            #
    #---------------------------------#
    ##### ----- Generate Panels ----- #####
    panels = dt.generate_panels([repanel])

    # rename axial coordinates for convenience in later plotting
    xn = panels.node[:, 1]
    xcp = panels.controlpoint[:, 1]

    # save specific example refinements
    if npan == 81
        f = open(dispath * "hub-coordinates-$(npan-1)-panels.dat", "w")
    else
        f = open(savepath * "hub-coordinates-$(npan-1)-panels.dat", "w")
    end
    for (x, r) in zip(repanel[:, 1], repanel[:, 2])
        write(f, "$x $r\n")
    end
    close(f)

    plot(repanel[:, 1], repanel[:, 2]; xlabel="x", ylabel="r", label="", aspectratio=1)
    plot(repanel[:, 1], repanel[:, 2]; xlabel="x", ylabel="r", label="", aspectratio=1)
    savefig(savepath * "hub-geometry.pdf")

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

    @time "AIC Matrix" begin
        # - Boundary on boundary influence coefficients - #
        AICn, AICt = dt.vortex_aic_boundary_on_boundary(
            panels.controlpoint,
            panels.normal,
            panels.tangent,
            panels.node,
            panels.nodemap,
            panels.influence_length,
        )
    end

    ## -- Manually Assemble Linear System -- ##
    # need to manually assemble since we don't have a way of automating the isolated case right now
    # initialize LHS matrix
    LHS = zeros(size(AICn, 2) + 1, size(AICn, 2) + 1)
    # fill in nominal AIC elements
    LHS[1:size(AICn, 1), 1:size(AICn, 2)] .= AICn
    # first node is on axis and is prescribed
    LHS[size(AICn, 1) + 1, 1] = 1.0
    # last node is on axis and is prescribed
    LHS[size(AICn, 1) + 2, size(AICn, 2)] = 1.0
    # need to add dummy variable due to second prescribed node
    LHS[1:size(AICn, 1), end] .= 1.0

    # - Freetream influence for RHS vector - #
    RHS = dt.freestream_influence_vector(panels.normal, repeat(Vs, panels.totpanel))
    push!(RHS, 0.0)
    push!(RHS, 0.0)

    #---------------------------------#
    #             Solving             #
    #---------------------------------#
    # standard linear solve
    gamb = LHS \ RHS
    # don't include dummy variable associated with internal panel
    gamb = gamb[1:(end - 1)]

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

    # save outputs for comparision
    f = open(savepath * "hub-xvvs-$(npan-1)-panels.jl", "w")
    write(f, "hubxvvs = [\n")
    for (x, v) in zip(panels.controlpoint[:, 1], Vtan)
        write(f, "$x $v\n")
    end
    write(f, "]")
    close(f)

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

    savefig(savepath * "hub-velocity-comp-$(npan-1)-panels.pdf")
    if npan == 81
        savefig(dispath * "hub-velocity-comp-$(npan-1)-panels.tikz")
    end

    # save refinement metric
    cpsums[i] = sum(cp .* panels.influence_length)
end

# look at relative error for example case and most refined case
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
annotate!(
    npans[id80] + 40, cpsums[id80], text("80 panels", 10, :left, :bottom; color=myred)
)

savefig(savepath * "hub-grid-refinement.pdf")
savefig(dispath * "hub-grid-refinement.tikz")

