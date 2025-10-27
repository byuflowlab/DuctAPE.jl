using DuctTAPE
const dt = DuctTAPE
include("../../plots_default_new.jl")

x = [0.0; 1.0]
y = [0.25; 0.0]
xy = [x y]
panel = dt.generate_panels(xy)

annoteoffset = 0.05

p = plot(aspectratio=1, axis=false)

plot!(panel.nodes[1,:,1],panel.nodes[1,:,2], marker=true, arrow=true,markersize=2, color=blue[1], label="")
plot!([panel.nodes[1,2,1]], [panel.nodes[1,2,2]], seriestype=:scatter, label="", color=blue[1], markersize=2)
annotate!(panel.nodes[1,1,1], annoteoffset+panel.nodes[1,1,2], text(L"p_1", 12, :center, color=blue[1]))
annotate!(panel.nodes[1,2,1], annoteoffset+panel.nodes[1,2,2], text(L"p_2", 12, :center, color=blue[1]))

normal_scale=0.25
normal_x = panel.control_point[1,1].+[0.0; normal_scale*panel.normal[1,1]]
normal_y = panel.control_point[1,2].+[0.0; normal_scale*panel.normal[1,2]]

plot!(normal_x, normal_y, arrow=true, color=red[1], label="")
annotate!(normal_x[2]+annoteoffset, normal_y[2]-annoteoffset, text(L"\hat{n}", 12, :center, color=red[1]))

plot!([panel.control_point[1,1]], [panel.control_point[1,2]], seriestype=:scatter, markersize=3, markershape=:rect, color=red[2], label="")
annotate!(normal_x[1], normal_y[1]-annoteoffset, text(L"\bar{p}", 12, :center, color=red[2]))


function makecircle(R,h,v,n; start=pi/2, stop=5*pi/2)
    t = range(start,stop,n)
    x = R*cos.(t) .+ h
    y = R*sin.(t) .+ v
    return x,y
end

x,y = makecircle(0.05,panel.nodes[1,1,1], panel.nodes[1,1,2],120, start=5*pi/6, stop=11*pi/5)
plot!(x[1:end-2],y[1:end-2],linecolor=1,arrow=true, color=blue[2], label="")
annotate!(panel.nodes[1,1,1], -1.5*annoteoffset+panel.nodes[1,1,2], text(L"\gamma_1=+\mu", 12, :left, color=blue[2]))

x,y = makecircle(0.05,panel.nodes[1,2,1], panel.nodes[1,2,2],120, stop=4*pi/6, start=11*pi/5)
plot!(x[1:end-2],y[1:end-2],linecolor=1,arrow=true,color=blue[2], label="")
annotate!(panel.nodes[1,2,1], -1.5*annoteoffset+panel.nodes[1,2,2], text(L"\gamma_2=-\mu", 12, :right, color=blue[2]))

savefig(p,"margin_doublet_panel.tikz")
# savefig(p,"margin_doublet_panel.pdf")
