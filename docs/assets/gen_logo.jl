using Plots
using FLOWMath
using DuctAPE
const dt = DuctAPE

include("geometry.jl")

##### ----- COLORS ----- #####

julia_blue = RGB(0.251, 0.388, 0.847)
julia_green = RGB(0.22, 0.596, 0.149)
julia_purple = RGB(0.584, 0.345, 0.698)
julia_red = RGB(0.796, 0.235, 0.2)

##### ----- GEOMETRY ----- #####

Rtip = 0.75
Rhub = 0.25
duct_chord = 2.0
duct_le_radius = 0.025
te_camber_angle = 9.0
wedge_angle = 10.0

# - Plotting Options - #
plot(; axis=false)
lw = 3
fa = 1 / 3

# blue rotor
rotorzloc = 0.35 * duct_chord
r = range(Rhub + 0.01, Rtip - 0.025, 11)
c = range(0.25, 0.2, 11) .* Rtip
t = range(70.0, 30.0, 11)

lez = rotorzloc .- c .* 0.25 .* sind.(t)
tez = rotorzloc .+ c .* 0.75 .* sind.(t)

plot!(
    lez,
    r;
    label="",
    color=julia_blue,
    # linewidth=lw,
    fillrange=Rhub * ones(11),
    fillcolor=julia_blue,
    fillalpha=fa,
)
plot!(
    tez,
    r;
    label="",
    color=julia_blue,
    # linewidth=lw,
    fillrange=Rhub * ones(11),
    fillcolor=julia_blue,
    fillalpha=fa,
)
plot!(
    [lez[end]; tez[end]],
    [Rtip; Rtip] .- 0.02;
    label="",
    linewidth=0,
    fillrange=Rhub * ones(11),
    fillcolor=julia_blue,
    fillalpha=fa,
)

# green duct
nz, nr, cz, cr, _, _ = duct_geom(
    Rtip, duct_chord, duct_le_radius, te_camber_angle, wedge_angle; duct_alpha=2, N=60
)
nr .+= Rtip
cr .+= Rtip

plot!(
    cz,
    cr;
    aspectratio=1,
    label="",
    fillrange=nr,
    color=julia_green,
    fillcolor=julia_green,
    fillalpha=fa,
    linewidth=lw,
)
plot!(nz, nr; label="", color=julia_green, linewidth=lw)

# red center body
cbz, cbr, _, _, _ = centerbody_geom(
    Rhub,
    duct_chord;
    cb_nc_le=0.125,
    cb_nc_stop=0.35,
    cb_tc_start=0.5,
    cb_ncp_z=0.125,
    cb_tcp_z=0.5,
    cb_te_r=0.0,
    cb_te_z=0.9,
    N=60,
    fmspline=(x, y) -> FLOWMath.Akima(x, y),
    smooth=true,
)

plot!(
    cbz[2:(end - 2)],
    cbr[2:(end - 2)];
    label="",
    color=julia_red,
    linewidth=lw,
    fillrange=zero(cbr),
    fillcolor=julia_red,
    fillalpha=fa,
)
cbr[end] = 0.0

# purple wake
# TODO: need to load DuctAPE and run the wake generation functions

# assemble propulsor
include("define_propulsor.jl")

# get wake geometry
problem_dimensions, prepost_containers, _, _, _, _, _, _ = dt.setup_analysis(
    propulsor, dt.set_options(; finterp=FLOWMath.akima)
)
wg = prepost_containers.wake_grid
for i in 2:2:size(wg, 3)
    plot!(wg[1, :, i], wg[2, :, i]; color=julia_purple, label="", linewidth=lw / 2)
end

# finish rotor bits
plot!(
    rotorzloc .* ones(length(wg[2, 1, 2:2:end])),
    wg[2, 1, 2:2:end];
    label="",
    color=julia_blue,
    markerstrokecolor=julia_blue,
    # markershape=:hline,
    markersize=3,
    linewidth=3,
    seriestype=:scatter,
)

# axis of rotation
plot!([nz[1], wg[1, end, 1]], -0.001 * ones(2); color=:black, label="", lw=1)
plot!(nz[13] * ones(2), [-0.05, 0.05]; color=:black, label="")
plot!(nz[13] * ones(2) .+ 0.03, [-0.05, 0.05]; color=:black, label="")
plot!(wg[1, end - 5, 1] * ones(2), [-0.05, 0.05]; color=:black, label="")
plot!(wg[1, end - 5, 1] * ones(2) .- 0.03, [-0.05, 0.05]; color=:black, label="")

##### ----- SAVE ----- #####
plot!(;grid=false, background_color=nothing)
savefig("assets/logo.svg")
