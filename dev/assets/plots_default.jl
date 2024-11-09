#=
Default Plots Settings for creating tikz native figures for dissertation
=#

using Plots
gr()
using LaTeXStrings
using Measures
using Printf

# - COLORS - #
byublue = RGB(0.0, 46.0 / 255, 93.0 / 255) #BYU Navy Blue
darkblue = RGB(0 / 255, 25 / 255, 50 / 255)
byured = RGB(155.0 / 255, 0, 0) #"BYU" Red
darkred = RGB(50 / 255, 0 / 255, 25 / 255)

primary = RGB(0, 92 / 255, 171 / 255) #blue
secondary = RGB(190 / 255, 76 / 255, 77 / 255) #red
tertiary = RGB(105 / 255, 174 / 255, 95 / 255) #green
quternary = RGB(167 / 255, 84 / 255, 164 / 255) #purple
quinary = RGB(190 / 255, 147 / 255, 61 / 255) #yellow
plotsgray = RGB(128 / 255, 128 / 255, 128 / 255) #gray

function plots_default(;
    background_color=nothing,
    fontfamily="cmunrm",
    foreground_color=RGB(128 / 255, 128 / 255, 128 / 255), #gray, plotsgray,
    size=(400, 300),
    annotationfontfamily="cmunrm",
    background_color_inside=nothing,
    background_color_legend=nothing,
    background_color_subplot=nothing,
    color_palette=[
        RGB(0, 92 / 255, 171 / 255), #blue
        RGB(190 / 255, 76 / 255, 77 / 255),
        RGB(105 / 255, 174 / 255, 95 / 255),
        RGB(167 / 255, 84 / 255, 164 / 255),
        RGB(190 / 255, 147 / 255, 61 / 255),
        RGB(128 / 255, 128 / 255, 128 / 255), #gray
    ],
    foreground_color_legend=nothing,
    margin=5mm,
    grid=false, # background grid true/false
    yguidefontrotation=-90,
)
    Plots.default()
    Plots.default(;
        background_color=background_color,
        fontfamily=fontfamily,
        foreground_color=foreground_color,
        size=size,
        annotationfontfamily=annotationfontfamily,
        background_color_inside=background_color_inside,
        background_color_legend=background_color_legend,
        background_color_subplot=background_color_subplot,
        color_palette=color_palette,
        foreground_color_legend=foreground_color_legend,
        margin=margin,
        grid=grid,
        yguidefontrotation=yguidefontrotation,
    )

    return (; color_palette)
end

custom_defaults = plots_default()
