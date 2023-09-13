#=
Default Plots Settings for creating tikz native figures for dissertation
=#

using Plots;
pgfplotsx();
# pyplot()
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
    fontfamily="Computer Modern",
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
    size=(400, 300) .* 5.0 ./ 4.0, #it appears that 100 ≈ 1inch in LaTeX
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
    markersize=1.5,
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
    annotationfontfamily="Computer Modern",
    annotationfontsize=10,
    #     annotationhalign,
    #     annotationrotation,
    #     annotations,
    #     annotationvalign,
    #     aspect_ratio,
    background_color_inside=nothing,
    background_color_legend=nothing,
    background_color_subplot=nothing,
    #     bottom_margin,
    #     camera,
    #     clims,
    color_palette=[
        #polynesian blue
        RGB(0 / 255, 75 / 255, 150 / 255),
        #004B96
        #argentinian blue
        RGB(100 / 255, 175 / 255, 250 / 255),
        #64AFFA
        #imperial red
        RGB(250 / 255, 75 / 255, 75 / 255),
        #FA4B4B
        #melon
        RGB(250 / 255, 175 / 255, 150 / 255),
        #FAAF96
        #lime green
        RGB(75 / 255, 200 / 255, 75 / 255),
        #4BC84B
        #light green
        RGB(150 / 255, 250 / 255, 150 / 255),
        #96FA96
    ],
    #color_palette=[
    #    RGB(0.0, 46.0 / 255.0, 93.0 / 255.0), #BYU Blue
    #    RGB(155.0 / 255.0, 0.0, 0.0), #"BYU" Red
    #    RGB(128.0 / 255.0, 128.0 / 255.0, 128.0 / 255.0), #Middle Gray
    #    RGB(162.0 / 255.0, 227.0 / 255.0, 162.0 / 255.0), #Light Green
    #    RGB(243.0 / 255.0, 209.0 / 255.0, 243.0 / 255.0), #Pink
    #    RGB(205.0 / 255.0, 179.0 / 255.0, 0.0), #Yellow
    #    RGB(161.0 / 255.0, 161.0 / 255.0, 226.0 / 255.0), #Purple
    #],
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
    margin=5mm,
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

myblue = [
    #oxford blue
    RGB(0 / 255, 25 / 255, 50 / 255),
    #01932
    #polynesian blue
    RGB(0 / 255, 75 / 255, 150 / 255),
    #004B96
    #argentinian blue
    RGB(100 / 255, 175 / 255, 250 / 255),
    #64AFFA
]
myred = [
    #dark purple (very dark red)
    RGB(50 / 255, 0 / 255, 25 / 255),
    #320019
    #imperial red
    RGB(250 / 255, 75 / 255, 75 / 255),
    #FA4B4B
    #melon
    RGB(250 / 255, 175 / 255, 150 / 255),
    #FAAF96
]
mygreen = [
    #dark green
    RGB(0 / 255, 50 / 255, 25 / 255),
    #003219
    #lime green
    RGB(75 / 255, 200 / 255, 75 / 255),
    #4BC84B
    #light green
    RGB(150 / 255, 250 / 255, 150 / 255),
    #96FA96
]
mygray = [
    #battleship gray
    RGB(150 / 255, 150 / 255, 150 / 255),
    #969696
    #davy's gray
    RGB(75 / 255, 75 / 255, 75 / 255),
    #4B4B4B
]

#mycolors = [
#    RGB(0, 46 / 255, 93 / 255), #BYU Blue
#    RGB(155 / 255, 0, 0), #"BYU" Red
#    RGB(128 / 255, 128 / 255, 128 / 255), #Middle Gray
#    RGB(162 / 255, 277 / 255, 162 / 255), #Light Green
#    RGB(243 / 255, 2090 / 255, 243 / 255), #Pink
#    RGB(205 / 255, 179 / 255, 0.0), #Yellow
#    RGB(161 / 255, 161 / 255, 226 / 255), #Purple
#]
