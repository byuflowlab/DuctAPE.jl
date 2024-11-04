#---------------------------------#
#         DISPATCH TYPES          #
#---------------------------------#

# Geometry
struct plotGeometry end
struct plotDuctGeometry end
struct plotBodyGeometry end
struct underlayGeometry end

# Surface Distributions
struct plotCP end
struct plotVtan end

# Boundary Layer
struct plotStagnation end
struct plotMomentum end

# Streamlines
struct plotStreamlines end

#---------------------------------#
#        PLOT Duct GEOMETRY        #
#---------------------------------#
@recipe function plot_duct_geometry(
    ::plotDuctGeometry,
    bvp;
    plot_panels=false,
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Default Labels
    xguide --> L"z"
    yguide --> L"r"

    # Aspect Ratio
    aspect_ratio --> 1
    xlim --> (
        minimum(bvp.node[1, 1:Int(bvp.npanel[1])]) - maximum(bvp.node[1, :]) * 0.1,
        maximum(bvp.node[1, :]) * 1.1,
    )

    # - Plot Body Geometry - #
    @series begin
        label --> false
        seriescolor --> 1
        if plot_panels
            marker --> true
            markersize --> 0.75
            markerstrokecolor --> 1
        end
        return bvp.node[1, Int(bvp.endnodeidxs[1, 1]):Int(bvp.endnodeidxs[2, 1])],
        bvp.node[2, Int(bvp.endnodeidxs[1, 1]):Int(bvp.endnodeidxs[2, 1])]
    end

    return nothing
end
# --------------------------------#
#       PLOT BODY GEOMETRY        #
#---------------------------------#
@recipe function plot_body_geometry(
    ::plotBodyGeometry,
    bvp;
    plot_panels=false,
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Default Labels
    xguide --> L"z"
    yguide --> L"r"

    # Aspect Ratio
    aspect_ratio --> 1

    # - Plot Body Geometry - #
    for b in 1:Int(bvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> 1
            if plot_panels
                marker --> true
                markersize --> 0.75
                markerstrokecolor --> 1
            end
            return bvp.node[1, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])],
            bvp.node[2, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])]
        end
    end

    return nothing
end

#---------------------------------#
#        PLOT ALL GEOMETRY        #
#---------------------------------#
@recipe function plot_geometry(
    ::plotGeometry,
    bvp,
    rsp,
    wvp;
    plot_panels=false,
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Default Labels
    xguide --> L"z"
    yguide --> L"r"

    # Aspect Ratio
    aspect_ratio --> 1
    ylim --> (0.0, maximum(bvp.node[2, :]) * 1.05)

    # - Plot Rotor Geometry - #
    for r in 1:Int(rsp.nbodies[])
        @series begin
            label --> false
            seriescolor --> 2
            if plot_panels
                marker --> true
                markersize --> 1.25
                markerstrokecolor --> 2
            end
            return rsp.node[1, Int(rsp.endnodeidxs[1, r]):Int(rsp.endnodeidxs[2, r])],
            rsp.node[2, Int(rsp.endnodeidxs[1, r]):Int(rsp.endnodeidxs[2, r])]
        end
    end

    # - Plot Body Geometry - #
    for b in 1:Int(bvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> 1
            if plot_panels
                marker --> true
                markersize --> 0.75
                markerstrokecolor --> 1
            end
            return bvp.node[1, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])],
            bvp.node[2, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])]
        end
    end

    # - Plot Wake Geometry - #
    for w in 1:Int(wvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> 3
            linewidth --> 0.5
            if plot_panels
                marker --> true
                markersize --> 0.5
                markerstrokecolor --> 3
            end
            return wvp.node[1, Int(wvp.endnodeidxs[1, w]):Int(wvp.endnodeidxs[2, w])],
            wvp.node[2, Int(wvp.endnodeidxs[1, w]):Int(wvp.endnodeidxs[2, w])]
        end
    end

    return nothing
end

#---------------------------------#
#       PLOT CP DISTRIBUTINS      #
#---------------------------------#
@recipe function plot_cp(
    ::plotCP,
    bvp,
    bout,
    rsp=nothing;
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
    cp_ylim =nothing
)


    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Flip y axis
    yflip --> true
    # Label Axes
    yguide --> L"C_P"
    xguide --> L"z"

    if !isnothing(rsp)
        xticks --> determine_geometry_xlabels(bvp, rsp)
        xgrid --> true
        ygrid --> false
    else
        grid --> false
    end

    if !isnothing(cp_ylim)
        ylim --> cp_ylim
    end

    # - Plot Body Surface Pressure - #
    for b in 1:Int(bvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> b
            linewidth = 2

            return bvp.controlpoint[
                1, Int(bvp.endpanelidxs[1, b]):Int(bvp.endpanelidxs[2, b])
            ],
            bout.cp_out[Int(bvp.endpanelidxs[1, b]):Int(bvp.endpanelidxs[2, b])]
        end
    end
end

#---------------------------------#
#      PLOT Vtan DISTRIBUTINS     #
#---------------------------------#
@recipe function plot_vtan(
    ::plotVtan,
    bvp,
    bout,
    Vref,
    rsp=nothing;
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
    labels=["Duct"; "Center Body"],
    vtan_ylim=nothing
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Label Axes
    yguide --> L"\frac{V_\mathrm{surf}}{V_\infty}"
    xguide --> L"z"
    if !isnothing(vtan_ylim)
        ylim --> vtan_ylim
    end

    if !isnothing(rsp)
        xticks --> determine_geometry_xlabels(bvp, rsp)
        xgrid --> true
        ygrid --> false
    else
        grid --> false
    end

    # - Plot Body Surface Pressure - #
    for b in 1:Int(bvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> b
            linewidth = 2
            return bvp.controlpoint[
                1, Int(bvp.endpanelidxs[1, b]):Int(bvp.endpanelidxs[2, b])
            ],
            bout.Vtan_out[Int(bvp.endpanelidxs[1, b]):Int(bvp.endpanelidxs[2, b])] ./ Vref
        end
    end
end

#---------------------------------#
#   SET XTICKS AT LE, TE, ETC.    #
#---------------------------------#
function determine_geometry_xlabels(bvp, rsp; tol=1e-2)
    xt = []
    xl = []

    # Rotor
    for r in 1:Int(rsp.nbodies[])
        push!(xt, rsp.node[1, Int(rsp.endnodeidxs[r, r])])
        push!(xl, round(rsp.node[1, Int(rsp.endnodeidxs[r, r])]; digits=2))
    end

    # Centerbody
    for i in 1:2
        tmp = bvp.node[1, Int(bvp.endnodeidxs[i, 2])]
        if !any(abs.(tmp .- xt) .< tol)
            push!(xl, round(tmp; digits=2))
        end
        if !any(abs.(tmp .- xt) .< tol / 10)
            push!(xt, tmp)
        end
    end

    # Duct
    for f in [maximum, minimum]
        tmp = f(bvp.node[1, Int(bvp.endnodeidxs[1, 1]):Int(bvp.endnodeidxs[2, 1])])
        if !any(abs.(tmp .- xt) .< tol)
            push!(xl, round(tmp; digits=2))
        end
        if !any(abs.(tmp .- xt) .< tol / 10)
            push!(xt, tmp)
        end
    end

    return (xt, xl)
end

#---------------------------------#
#PLOT GEOMETRY UNDER DISTRIBUTIONS#
#---------------------------------#
@recipe function plot_overlayed_geometry(
    ::underlayGeometry,
    bvp,
    rsp;
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
    discrete_labels=true,
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Plot specific values
    axis --> false
    ticks --> false
    aspect_ratio --> 1

    # - Plot Body Geometry - #
    for b in 1:Int(bvp.nbodies[])
        @series begin
            label --> false
            seriescolor --> b
            linealpha --> 1 / 2
            linewidth --> 0.5
            return bvp.node[1, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])],
            bvp.node[2, Int(bvp.endnodeidxs[1, b]):Int(bvp.endnodeidxs[2, b])]
        end
    end

    # - Close Center Body - #
    @series begin
        label --> false
        seriescolor --> Int(bvp.nbodies[])
        linealpha --> 1 / 2
        linewidth --> 0.5
        return ones(2) * bvp.node[1, Int(bvp.endnodeidxs[2, 2])],
        [0.0; bvp.node[2, Int(bvp.endnodeidxs[2, 2])]]
    end

    # - Plot Rotor Geometry - #
    for r in 1:Int(rsp.nbodies[])
        @series begin
            z_order --> 1
            label --> false
            seriescolor --> 3
            linewidth --> 1.5
            return rsp.node[1, Int(rsp.endnodeidxs[1, r]):Int(rsp.endnodeidxs[2, r])],
            rsp.node[2, Int(rsp.endnodeidxs[1, r]):Int(rsp.endnodeidxs[2, r])]
        end
    end
end

#---------------------------------#
#         PLOT STAGNATION         #
#---------------------------------#
@recipe function plot_stagnation(
    ::plotStagnation,
    blo,
    bvp,
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
)
    @series begin
        seriestype --> :scatter
        color --> 2
        markerstrokealpha --> 0
        label --> ""

        stagz = sum(
            bvp.controlpoint[1, blo.stagnation_indices] .*
            [1.0 - blo.split_ratio; blo.split_ratio],
        )
        stagr = sum(
            bvp.controlpoint[2, blo.stagnation_indices] .*
            [1.0 - blo.split_ratio; blo.split_ratio],
        )
        return [stagz], [stagr]
    end

    return nothing
end

#---------------------------------#
#              PLOT δ₂            #
#---------------------------------#
@recipe function plot_momentum_thickness(
    ::plotMomentum,
    blo,
    bvp;
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
    scale_thickness=50.0,
    bl_ylim=nothing
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Plot specific values
    aspect_ratio --> 1

    # stagnation point
    stagz = sum(
        bvp.controlpoint[1, blo.stagnation_indices] .*
        [1.0 - blo.split_ratio; blo.split_ratio],
    )
    stagr = sum(
        bvp.controlpoint[2, blo.stagnation_indices] .*
        [1.0 - blo.split_ratio; blo.split_ratio],
    )

    # stagnation point normal
    stagnz = sum(
        bvp.normal[1, blo.stagnation_indices] .* [1.0 - blo.split_ratio; blo.split_ratio]
    )
    stagnr = sum(
        bvp.normal[2, blo.stagnation_indices] .* [1.0 - blo.split_ratio; blo.split_ratio]
    )

    if !isnothing(bl_ylim)
        ylim --> bl_ylim
    end

    @series begin
        color --> 1
        label --> false
        # linewidth --> 0.5
        if isapprox(blo.separation_point_ratio_upper, 1.0)
            arrow --> :closed
        end

        # spline the states vs steps
        d2_upper = blo.upper_solved_states[1, :]

        # spline the control points vs surface length
        cpzu = akima_smooth(
            blo.surface_length_upper,
            [stagz; bvp.controlpoint[1, blo.stagnation_indices[2]:Int(bvp.npanel[1])]],
            blo.upper_solved_steps,
        )
        cpru = akima_smooth(
            blo.surface_length_upper,
            [stagr; bvp.controlpoint[2, blo.stagnation_indices[2]:Int(bvp.npanel[1])]],
            blo.upper_solved_steps,
        )

        nhatzu = akima_smooth(
            blo.surface_length_upper,
            [stagnz; bvp.normal[1, blo.stagnation_indices[2]:Int(bvp.npanel[1])]],
            blo.upper_solved_steps,
        )
        nhatru = akima_smooth(
            blo.surface_length_upper,
            [stagnr; bvp.normal[2, blo.stagnation_indices[2]:Int(bvp.npanel[1])]],
            blo.upper_solved_steps,
        )

        # thickness point is thickness normal to control points
        d2zu = cpzu .+ (d2_upper .* nhatzu) .* scale_thickness
        d2ru = cpru .+ (d2_upper .* nhatru) .* scale_thickness

        return d2zu, d2ru
    end

    @series begin
        color --> 2
        label --> false
        # linewidth --> 0.5
        if isapprox(blo.separation_point_ratio_lower, 1.0)
            arrow --> :closed
        end

        # spline the states vs steps
        d2_lower = blo.lower_solved_states[1, :]

        # spline the control points vs surface length
        cpzl = akima_smooth(
            blo.surface_length_lower,
            [bvp.controlpoint[1, blo.stagnation_indices[2]:-1:1]; stagz],
            blo.lower_solved_steps,
        )
        cprl = akima_smooth(
            blo.surface_length_lower,
            [bvp.controlpoint[2, blo.stagnation_indices[2]:-1:1]; stagr],
            blo.lower_solved_steps,
        )

        nhatzl = akima_smooth(
            blo.surface_length_lower,
            [bvp.normal[1, blo.stagnation_indices[2]:-1:1]; stagnz],
            blo.lower_solved_steps,
        )
        nhatrl = akima_smooth(
            blo.surface_length_lower,
            [bvp.normal[2, blo.stagnation_indices[2]:-1:1]; stagnr],
            blo.lower_solved_steps,
        )

        # thickness point is thickness normal to control points
        d2zl = cpzl .+ (d2_lower .* nhatzl) .* scale_thickness
        d2rl = cprl .+ (d2_lower .* nhatrl) .* scale_thickness

        return d2zl, d2rl
    end

    return nothing
end

#---------------------------------#
#         PLOT STREAMLINES        #
#---------------------------------#
@recipe function plot_streamlines(
    ::plotStreamlines,
    bvp,
    bvs,
    wvp,
    wvs,
    rsp,
    rss,
    vinf;
    starting_radial_points=range(0.001, 1.0; length=20),
    axial_range=[-0.2, 1.0],
    nominal_step_size=1e-2,
    step_limit=Int(1e2),
    dot_tol=0.999,
    stag_tol=0.4,
    integration_options=IntegrationOptions(),
    default_colors=(;
        primary=RGB(1 / 255, 149 / 255, 226 / 255), #blue
        secondary=RGB(189 / 255, 10 / 255, 53 / 255), #red
        tertiary=RGB(76 / 255, 173 / 255, 59 / 255), #green
        quaternary=RGB(238 / 255, 167 / 255, 46 / 255), #orange
        quinary=RGB(155 / 255, 82 / 255, 162 / 255), #purple
        plotsgray=RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ),
    scale_thickness=50.0,
)
    color_palette --> [
        default_colors.primary,
        default_colors.secondary,
        default_colors.tertiary,
        default_colors.quaternary,
        default_colors.quinary,
        default_colors.plotsgray,
    ]

    # Plot specific values
    aspect_ratio --> 1

    p, _ = DuctAPE.calculate_streamlines(
        bvp,
        bvs,
        wvp,
        wvs,
        rsp,
        rss,
        vinf;
        starting_radial_points=starting_radial_points,
        axial_range=axial_range,
        step_limit=step_limit,
        nominal_step_size=nominal_step_size,
        integration_options=integration_options,
        stag_tol=stag_tol,
    )

    for (z, r) in zip(eachcol(p[1, :, :]), eachcol(p[2, :, :]))
        @series begin
            color --> default_colors.plotsgray
            label --> ""
            linewidth --> 0.5
            return z, r
        end
    end

    return nothing
end
