include("../../plots_default.jl")

plots_default(;
    size=(200, 200),
    linewidth=2.0,
    markerstrokealpha=1,
    markerstrokewidth=2,
    markersize=2,
    legend=false, # include legend true/false
    ylims=(0, 3),
    xlims=(0, 2),
    ticks=false, #turns off tick marks
)

savepath = joinpath("prelims", "figures/")

## -- Constant Source Distribution Plot
x = [0; 1.75]
y = [1.25; 1.25]

plot1 = plot(
    x,
    y;
    fillrange=[0; 0],
    fillcolor=1,
    fillalpha=0.125,
    xlabel="s",
    ylabel=L"$q(s)=$ const",
);

savefig(plot1, joinpath(savepath, "constantsourcedist.tikz"))

## -- Linear Source Distribution Plot

x = [0; 1.75]
y = [1.25; 1.25]

plot2 = plot(x, y; fillrange=[0; 0], fillcolor=1, xlabel="s", ylabel=L"$q(s)=$ linear");

y2 = [1.25; 2.25]
plot!(x, y2; fillrange=[1.25; 1.25], fillcolor=2);

savefig(plot2, joinpath(savepath, "linearsourcedist.tikz"))

## -- Uniform Flow

x = [0; 2]
y = [0.05; 0.25]

uplot = plot(
    x,
    y;
    linecolor=1,
    showaxis=false,
    arrow=arrow(7.0),
    markerstrokealpha=1,
    markerstrokecolor=1,
    markercolor=1,
)
for i in 1:7
    plot!(
        x,
        (i * 0.35) .+ y;
        linecolor=1,
        arrow=arrow(7.0),
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
end

savefig(uplot, joinpath(savepath, "uniformflow.tikz"))

## -- Source Flow
function makecircle(R, h, v, n)
    t = range(pi / 2; stop=5 * pi / 2, length=n)
    x = R * cos.(t) .+ h
    y = R * sin.(t) .+ v

    return x, y
end

x, y = makecircle(1.5, 1.5, 1.5, 12)

s1plot = plot(
    [1.5],
    [1.5];
    seriestype=:scatter,
    markershape=:circle,
    showaxis=false,
    aspect_ratio=1,
    xlims=(0, 3),
    markerstrokealpha=1,
    markerstrokecolor=1,
    markercolor=1,
);
for i in 1:12
    plot!(
        [1.5; x[i]],
        [1.5, y[i]];
        linecolor=1,
        arrow=arrow(7.0),
        aspect_ratio=1,
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
end

savefig(s1plot, joinpath(savepath, "sourceflow.tikz"))

## -- Sink Flow
x1, y1 = makecircle(1.5, 1.5, 1.5, 12)
x2, y2 = makecircle(0.35, 1.5, 1.5, 12)

s1plot = plot(
    [1.5],
    [1.5];
    seriestype=:scatter,
    markershape=:circle,
    showaxis=false,
    aspect_ratio=1,
    xlims=(0, 3),
    markerstrokealpha=1,
    markerstrokecolor=1,
    markercolor=1,
);
for i in 1:12
    plot!(
        [x1[i]; x2[i]],
        [y1[i], y2[i]];
        linecolor=1,
        arrow=arrow(7.0),
        aspect_ratio=1,
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
end

savefig(s1plot, joinpath(savepath, "sinkflow.tikz"))

## -- Vortex Flow

vplot = plot(
    [1.5],
    [1.5];
    seriestype=:scatter,
    markershape=:circle,
    showaxis=false,
    aspect_ratio=1,
    xlims=(0, 3),
    markerstrokealpha=1,
    markerstrokecolor=1,
    markercolor=1,
);
for i in 1:4
    x, y = makecircle(i * 0.35, 1.5, 1.5, 120)
    plot!(
        x[1:(end - 2)],
        y[1:(end - 2)];
        linecolor=1,
        arrow=arrow(7.0),
        aspect_ratio=1,
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
end

savefig(vplot, joinpath(savepath, "vortexflow.tikz"))

## -- doublet Flow

dplot = plot(
    [0],
    [0];
    seriestype=:scatter,
    markershape=:circle,
    showaxis=false,
    aspect_ratio=1,
    xlims=(-1.5, 1.5),
    ylims=(-1.5, 1.5),
    markercolor=1,
    markerstrokecolor=1,
);
for i in 1:4
    x, y = makecircle(i * 0.15, 0, i * 0.15, 36)
    plot!(
        -x[1:(end - 1)],
        y[1:(end - 1)];
        linecolor=1,
        arrow=arrow(7.0),
        aspect_ratio=1,
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
    plot!(
        -x[1:(end - 1)],
        -y[1:(end - 1)];
        linecolor=1,
        arrow=arrow(7.0),
        aspect_ratio=1,
        markerstrokealpha=1,
        markerstrokecolor=1,
        markercolor=1,
    )
end

savefig(dplot, joinpath(savepath, "doubletflow.tikz"))
