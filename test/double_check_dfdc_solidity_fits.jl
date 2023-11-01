
# - Get Project Directory - #
project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

# - load DuctTAPE - #
using DuctTAPE
const dt = DuctTAPE

# - load plotting defaults - #
include(project_dir * "/visualize/plots_default.jl")
savepath = project_dir * "/test/figures/"
include(project_dir * "/test/data/wallis_fig6-29_data.jl")

pss = plot(
    oneoversig05[:, 1],
    oneoversig05[:, 2];
    seriestype=:scatter,
    label="Wallis Model",
    xlabel="Stagger Angle",
    ylabel=L"\frac{c_\ell}{c_{\ell_i}}",
    color=myblue,
)

plot!(
    pss, oneoversig06[:, 1], oneoversig06[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig07[:, 1], oneoversig07[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig08[:, 1], oneoversig08[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig09[:, 1], oneoversig09[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig10[:, 1], oneoversig10[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig11[:, 1], oneoversig11[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig12[:, 1], oneoversig12[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig13[:, 1], oneoversig13[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig14[:, 1], oneoversig14[:, 2]; seriestype=:scatter, label="", color=myblue
)

plot!(
    pss, oneoversig15[:, 1], oneoversig15[:, 2]; seriestype=:scatter, label="", color=myblue
)

N = 100
solids = 1.0 ./ [0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 1.1; 1.2; 1.3; 1.4; 1.5]
stags = range(0.0, 100.0, N) * pi / 180.0
# solids = range(0.0,5.0,N)
# stags = range(0.0, 180.0, N) * pi / 180.0
dfdcvals = zeros(length(stags), length(solids))
smoothvals = zeros(length(stags), length(solids))
for (is, solidity) in enumerate(solids)
    for (i, stagger) in enumerate(stags)
        dfdcvals[i, is] = dt.solidityandstaggerfactor(solidity, stagger)
        smoothvals[i, is] = dt.solidityandstaggerfactorsmooth(solidity, stagger)
    end
    plot!(
        pss,
        stags * 180.0 / pi,
        dfdcvals[:, is];
        color=myred,
        label=is == 1 ? "DFDC Fits" : "",
    )
    plot!(
        pss,
        stags * 180.0 / pi,
        smoothvals[:, is];
        color=mygreen,
        label=is == 1 ? "Smoothed Fits" : "",
    )
end
plot!()
savefig(savepath * "solidity_correction_double_check.pdf")
