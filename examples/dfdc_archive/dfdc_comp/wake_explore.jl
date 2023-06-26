
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

savepath = project_dir * "/examples/dfdc_comp/"

include(project_dir * "/plots_default.jl")

using DuctTAPE
const dt = DuctTAPE

include(project_dir * "/examples/dfdc_comp/DFDC_WAKE_STRENGTHS.jl")
# column 1: control point index (where QC is)
# column 2: x location for QC
# column 3: r location for QC
# column 4: Vm at QC point
# column 5: wake panel index (where gamma is)
# column 6: x location for gamma
# column 7: r location for gamma
# column 8: Vmavg for gamma
# column 9: gamma value

#---------------------------------#
#          DFDC Geometry          #
#---------------------------------#
ms = 1.0
lw = 0.5

## -- Duct and Hub -- ##
# - "control points" - #
plot(; xlabel="x", ylabel="r", aspectratio=1)
plot!(
    WAKE2[:, 2],
    WAKE2[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[1],
    label="Duct Control Points",
)
plot!(
    WAKE1[:, 2],
    WAKE1[:, 3];
    seriestype=:scatter,
    color=mycolors[2],
    markersize=ms,
    markershape=:circ,
    label="Hub Control Points",
)

# - "panel nodes" - #
plot!(
    WAKE2[:, 6],
    WAKE2[:, 7];
    linewidth=lw,
    markershape=:rect,
    markersize=ms,
    color=mycolors[4],
    label="Duct Panel Nodes",
)
plot!(
    WAKE1[:, 6],
    WAKE1[:, 7];
    linewidth=lw,
    color=mycolors[7],
    markershape=:rect,
    markersize=ms,
    label="Hub Panel Nodes",
)

savefig(savepath * "dfdc_hubductgeom.pdf")

## -- Rotor -- ##
# - "control points" - #
plot!(
    WAKE3[:, 2],
    WAKE3[:, 3];
    seriestype=:scatter,
    color=mycolors[5],
    markersize=ms,
    markershape=:circ,
    label="Rotor Control Points",
)

# - "panel nodes" - #
plot!(
    WAKE3[:, 6],
    WAKE3[:, 7];
    linewidth=lw,
    markershape=:rect,
    markersize=ms,
    color=:black,
    label="Rotor Panel Nodes",
)
savefig(savepath * "dfdc_hubductrotorgeom.pdf")

## -- Duct and Hub Wakes -- ##

# - "control points" - #
plot!(
    WAKE14[:, 2],
    WAKE14[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[4],
    label="Duct Wake Control Points",
)
plot!(
    WAKE4[:, 2],
    WAKE4[:, 3];
    seriestype=:scatter,
    color=mycolors[7],
    markersize=ms,
    markershape=:circ,
    label="Hub Wake Control Points",
)

# - "panel nodes" - #
plot!(
    WAKE14[:, 6],
    WAKE14[:, 7];
    linewidth=lw,
    markershape=:rect,
    markersize=ms,
    color=mycolors[1],
    label="Duct Wake Panel Nodes",
)
plot!(
    WAKE4[:, 6],
    WAKE4[:, 7];
    linewidth=lw,
    color=mycolors[2],
    markershape=:rect,
    markersize=ms,
    label="Hub Wake Panel Nodes",
)

savefig(savepath * "dfdc_rotorbodybodywakegeom.pdf")

## -- Intermediate Wakes -- ##

# - "control points" - #
plot!(
    WAKE5[:, 2],
    WAKE5[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="Wake Control Points",
)
plot!(
    WAKE6[:, 2],
    WAKE6[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE7[:, 2],
    WAKE7[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE8[:, 2],
    WAKE8[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE9[:, 2],
    WAKE9[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE10[:, 2],
    WAKE10[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE11[:, 2],
    WAKE11[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE12[:, 2],
    WAKE12[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)
plot!(
    WAKE13[:, 2],
    WAKE13[:, 3];
    seriestype=:scatter,
    markersize=ms,
    markershape=:circ,
    color=mycolors[3],
    label="",
)

# - "panel nodes" - #
plot!(
    WAKE5[:, 6],
    WAKE5[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="Wake Panel Nodes",
)
plot!(
    WAKE6[:, 6],
    WAKE6[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE7[:, 6],
    WAKE7[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE8[:, 6],
    WAKE8[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE9[:, 6],
    WAKE9[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE10[:, 6],
    WAKE10[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE11[:, 6],
    WAKE11[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE12[:, 6],
    WAKE12[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)
plot!(
    WAKE13[:, 6],
    WAKE13[:, 7];
    linewidth=lw,
    color=mycolors[6],
    markershape=:rect,
    markersize=ms,
    label="",
)

savefig(savepath * "dfdc_fullgeom.pdf")

