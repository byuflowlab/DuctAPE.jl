#=
Default Plots Settings for creating tikz native figures for dissertation
=#

using Plots;
pgfplotsx();
# pythonplot()
using LaTeXStrings
using Measures

default()
default(;
    #     #:Plot
    # background_color=RGBA(1, 1, 1, 0),
    # background_color = nothing,
    # background_color_outside = nothing,
    #     display_type,
    #     dpi,
    #     extra_kwargs,
    #     extra_plot_kwargs,
    fontfamily="Monaco",
    #     foreground_color,
    #     html_output_format,
    #     inset_subplots,
    #     layout,
    #     link,
    #     overwrite_figure,
    #     plot_title,
    #     plot_title_location,
    #     plot_titlefontcolor,
    #     plot_titlefontfamily,
    #     plot_titlefonthalign,
    #     plot_titlefontrotation,
    #     plot_titlefontsize,
    #     plot_titlefontvalign,
    #     pos,
    #     show,
    # size=(400, 300), #it appears that 100 ≈ 1inch in LaTeX
    size=(300, 225), #it appears that 100 ≈ 1inch in LaTeX
    # size=(600, 450), #it appears that 100 ≈ 1inch in LaTeX
    # size=(800, 600), #it appears that 100 ≈ 1inch in LaTeX
    #     tex_output_standalone,
    #     thickness_scaling,
    #     warn_on_unsupported,
    #     window_title,

    #######################
    #      :Series
    #######################
    #     arrow,
    #     bar_edges,
    #     bar_position,
    #     bar_width,
    #     bins,
    #     colorbar_entry,
    #     connections,
    #     contour_labels,
    #     contours,
    #     extra_kwargs,
    #     fill_z,
    fillalpha=0.125,
    fillcolor=RGB(128 / 255, 128 / 255, 128 / 255),
    #     fillrange,
    #     group,
    #     hover,
    #     label,
    #     levels,
    #     line_z,
    #     linealpha,
    #     linecolor,
    #     linestyle,
    linewidth=1.0,
    #     marker_z,
    #     markeralpha,
    #     markercolor,
    #     markershape,
    #     markersize,
    markerstrokealpha=0,
    #     markerstrokecolor,
    #     markerstrokestyle,
    #     markerstrokewidth,
    #     normalize,
    #     orientation,
    #     primary,
    #     quiver,
    #     ribbon,
    #     series_annotations,
    #     seriesalpha,
    #     seriescolor,
    #     seriestype,
    #     show_empty_bins,
    #     smooth,
    #     stride,
    #     subplot,
    #     weights,
    #     x,
    #     xerror,
    #     y,
    #     yerror,
    #     z,
    #     zerror

    #######################
    #      :Subplot
    #######################
    #     annotationcolor,
    annotationfontfamily="Monaco",
    annotationfontsize=10,
    annotationhalign=:left,
    #     annotationrotation,
    #     annotations,
    annotationvalign=:bottom,
    #     aspect_ratio,
    background_color_inside=nothing,
    background_color_legend=nothing,
    background_color_subplot=nothing,
    #     bottom_margin,
    #     camera,
    #     clims,
    color_palette=[
        RGB(0.0 / 255, 92.0 / 255, 171.0 / 255) # royal blue
        RGB(192.0 / 255, 83.0 / 255, 103.0 / 255) # royal red
        RGB(143.0 / 255, 166.0 / 255, 81.0 / 255) # royal green
        RGB(130.0 / 255, 130.0 / 255, 130.0 / 255) # royal gray
    ],
    #     colorbar,
    #     colorbar_continuous_values,
    #     colorbar_discrete_values,
    #     colorbar_fontfamily,
    #     colorbar_formatter,
    #     colorbar_scale,
    #     colorbar_tickfontcolor,
    #     colorbar_tickfontfamily,
    #     colorbar_tickfonthalign,
    #     colorbar_tickfontrotation,
    #     colorbar_tickfontsize,
    #     colorbar_tickfontvalign,
    #     colorbar_ticks,
    #     colorbar_title,
    #     colorbar_title_location,
    #     colorbar_titlefontcolor,
    #     colorbar_titlefontfamily,
    #     colorbar_titlefonthalign,
    #     colorbar_titlefontrotation,
    #     colorbar_titlefontsize,
    #     colorbar_titlefontvalign,
    #     extra_kwargs,
    #     fontfamily_subplot,
    foreground_color_legend=nothing,
    #     foreground_color_subplot,
    #     foreground_color_title,
    #     framestyle = :zerolines,
    #     left_margin,
    # legend=false, # include legend true/false
    #     legendfontcolor,
    #     legendfontfamily,
    #     legendfonthalign,
    #     legendfontrotation,
    #     legendfontsize,
    #     legendfontvalign,
    #     legendtitle,
    #     legendtitlefontcolor,
    #     legendtitlefontfamily,
    #     legendtitlefonthalign,
    #     legendtitlefontrotation,
    #     legendtitlefontsize,
    #     legendtitlefontvalign,
    margin=0mm,
    #     projection,
    #     right_margin,
    #     subplot_index,
    #     title,
    #     titlefontcolor,
    #     titlefontfamily,
    #     titlefonthalign,
    #     titlefontrotation,
    #     titlefontsize,
    #     titlefontvalign,
    #     titlelocation,
    #     top_margin

    #####################
    #       :Axis
    #####################
    #     discrete_values,
    #     draw_arrow,
    #     flip,
    #     foreground_color_axis,
    #     foreground_color_border,
    #     foreground_color_grid,
    #     foreground_color_guide,
    #     foreground_color_minor_grid,
    #     foreground_color_text,
    #     formatter,
    grid=false, # background grid true/false
    #     gridalpha,
    # gridlinewidth=0.5,
    #     gridstyle,
    #     guide,
    #     guide_position,
    #     guidefontcolor,
    #     guidefontfamily,
    #     guidefonthalign,
    #     guidefontrotation,
    #     guidefontsize,
    #     guidefontvalign,
    # ylims=(0, 3),
    # xlims=(0, 2),
    #     link,
    #     minorgrid,
    #     minorgridalpha,
    #     minorgridlinewidth,
    #     minorgridstyle,
    #     minorticks,
    #     mirror,
    #     rotation,
    # scale,
    # showaxis = false, #turns off spines and tick labels, but not ticks
    #     tick_direction,
    #     tickfontcolor,
    #     tickfontfamily,
    #     tickfonthalign,
    #     tickfontrotation,
    #     tickfontsize,
    #     tickfontvalign,
    # ticks=false, #turns off tick marks
    #     widen,
)

byublue = RGB(0.0, 46.0 / 255, 93.0 / 255) #BYU Blue
darkblue = RGB(0 / 255, 25 / 255, 50 / 255)
byured = RGB(155.0 / 255, 0, 0) #"BYU" Red
darkred = RGB(50 / 255, 0 / 255, 25 / 255)
middlegray = RGB(128.0 / 255, 128.0 / 255, 128.0 / 255) #Middle Gray
myblue = RGB(0.0 / 255, 92.0 / 255, 171.0 / 255) # royal blue
myred = RGB(192.0 / 255, 83.0 / 255, 103.0 / 255) # royal red
mygreen = RGB(143.0 / 255, 166.0 / 255, 81.0 / 255) # royal green
mygrey = RGB(130.0 / 255, 130.0 / 255, 130.0 / 255) # royal gray

"""
"""
function tightplot(p; spidx=1)

    # Get figsize
    figsize = p.attr[:size]

    #get width and height of plot area based on min and max
    # xmin, xmax = Plots.ignorenan_extrema(p)
    # ymin = minimum([Plots.ignorenan_minimum(series[:y]) for series in p.series_list])
    # ymax = maximum([Plots.ignorenan_maximum(series[:y]) for series in p.series_list])
    xmin, xmax = xlims(p)
    ymin, ymax = ylims(p)

    width = xmax - xmin
    height = ymax - ymin
    cropratio = height / width

    if width > height
        newsize = (figsize[1], figsize[1] * cropratio)
    elseif width < height
        newsize = (figsize[2] / cropratio, figsize[2])
    else
        newsize = (figsize[1], figsize[1])
    end

    plot!(p; size=newsize)

    return nothing
end

"""
"""
function savetightplot(p, filename; spidx=1)
    tightplot(p; spidx=spidx)
    savefig(filename)
    return nothing
end
