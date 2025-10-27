using DuctTAPE
const dt = DuctTAPE
include("../../plots_default_new.jl")

x = [0.0; 1.0]
y = [0.25; 0.0]
xy = [x y]
panel = dt.generate_panels(xy)

annoteoffset = 0.05

p = plot(aspectratio=1, axis=false)

plot!(panel.nodes[1,:,1],panel.nodes[1,:,2], marker=true, arrow=9,markersize=4, color=byublue[1], label="")
plot!([panel.nodes[1,2,1]], [panel.nodes[1,2,2]], seriestype=:scatter, label="", color=byublue[1], markersize=4)
annotate!(panel.nodes[1,1,1], annoteoffset+panel.nodes[1,1,2], text(L"p_1", 12, :center, color=byublue[1]))
annotate!(panel.nodes[1,2,1], annoteoffset+panel.nodes[1,2,2], text(L"p_2", 12, :center, color=byublue[1]))

normal_scale=0.25
normal_x = panel.controlpoint[1,1].+[0.0; normal_scale*panel.normal[1,1]]
normal_y = panel.controlpoint[1,2].+[0.0; normal_scale*panel.normal[1,2]]

plot!(normal_x, normal_y, arrow=7, color=byured[1], label="")
annotate!(normal_x[2]+annoteoffset, normal_y[2]-annoteoffset, text(L"\hat{n}", 12, :center, color=byured[1]))

plot!([panel.controlpoint[1,1]], [panel.controlpoint[1,2]], seriestype=:scatter, markersize=5, markershape=:rect, color=byured[3], label="")
annotate!(normal_x[1], normal_y[1]-annoteoffset, text(L"\overline{p}", 12, :center, color=byured[3]))


function makecircle(R,h,v,n; start=pi/2, stop=5*pi/2)
    t = range(start,stop,n)
    x = R*cos.(t) .+ h
    y = R*sin.(t) .+ v
    return x,y
end

x,y = makecircle(0.1,panel.controlpoint[1,1], panel.controlpoint[1,2],120, start=5*pi/6, stop=11*pi/5)
plot!(x[1:end-2],y[1:end-2],linecolor=byublue[2],arrow=7, color=byublue[2], label="")
annotate!(panel.controlpoint[1,1], -2.5*annoteoffset+panel.controlpoint[1,2], text(L"\gamma", 12, :left, color=byublue[2]))

savefig(p,"constant_vortex_panel.tikz")
# savefig(p,"margin_doublet_panel.pdf")
