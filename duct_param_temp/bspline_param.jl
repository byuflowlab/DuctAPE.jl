#=

Attempt to make a Bspline parameterization for duct cross section

=#

using Splines
using Plots

function generate_duct_section(
    Rin,
    Rd,
    Rout,
    Lin,
    Ltot,
    r_le,
    wedge_angle,
    te_camber_angle,
    max_thick;
    N=161,
    disk_ratio=0.125,
)
    boattailangle = pi * wedge_angle / 360 #boattail angle is half of wedge angle

    te_camber_angle *= pi / 180

    #--Initialize knot vector
    knots_inner = [
        0.0,
        0.0,
        0.0,
        Lin * (1.0 - disk_ratio) / Ltot,
        Lin * (1.0 - disk_ratio) / Ltot,
        Lin * (1.0 + disk_ratio) / Ltot,
        Lin * (1.0 + disk_ratio) / Ltot,
        1.0,
        1.0,
        1.0,
    ]

    knots_outer = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    #calculate unweighted controlpoints
    Pouter = [
        [0.0, Rin], [0.0, Rin + (Rin - Rd)], [Lin, Rd + 2.0 * max_thick], [Ltot, Rout]
    ]

    Pinner = [
        [0.0, Rin], #sqrt(2.0 * leradius) / 3.0],
        [0.0, Rd],
        [Lin * (1.0 - disk_ratio), Rd],
        [Lin * (1.0 - disk_ratio), Rd],
        [Lin * (1.0 + disk_ratio), Rd],
        # [Lin * 1.1, Rd],
        [(Ltot + Lin) / 2.0, Rd],
        [Ltot, Rout],
    ]

    inner_spline = Splines.NURBS(2, knots_inner, ones(length(Pinner)), Pinner)
    outer_spline = Splines.NURBS(3, knots_outer, ones(length(Pouter)), Pouter)

    u = range(0.0, 1.0; length=N + 1)
    n = length(inner_spline.ctrlpts) - 1
    Cw_inner = [[0.0; 0.0] for i in 1:length(u)]
    Cw_outer = [[0.0; 0.0] for i in 1:length(u)]

    #loop through parametric points to get curve points
    for i in 1:length(u)
        Cw_inner[i] = Splines.curvepoint(inner_spline, u[i])
        Cw_outer[i] = Splines.curvepoint(outer_spline, u[i])
    end

    return Pinner, Pouter, Cw_inner, Cw_outer
end

Rin = 2.75
Rd = 2.5
Rout = 2.35
Ltot = Rd * 2.0
Lin = 0.4 * Ltot
r_le = 0.25
wedge_angle = 20.0
te_camber_angle = 5.0
max_thick = 0.5
N = 161
disk_ratio = 0.125

pin, pout, xrin, xrout = generate_duct_section(
    Rin,
    Rd,
    Rout,
    Lin,
    Ltot,
    r_le,
    wedge_angle,
    te_camber_angle,
    max_thick;
    N=N,
    disk_ratio=disk_ratio,
)

plot(; aspectratio=:equal)
plot!(getindex.(pin, 1), getindex.(pin, 2); seriestype=:scatter)
plot!(getindex.(pout, 1), getindex.(pout, 2); seriestype=:scatter)
plot!(getindex.(xrout, 1), getindex.(xrout, 2))
plot!(getindex.(xrin, 1), getindex.(xrin, 2))

using FLOWMath
using FLOWFoil
t = Ltot .* FLOWFoil.cosine_spacing(161)
osp = akima(getindex.(xrout, 1), getindex.(xrout, 2), t)
isp = akima(getindex.(xrin, 1), getindex.(xrin, 2), t)

xx = [reverse(t); t[2:end]]
rr = [reverse(isp); osp[2:end]]
if rr[1] != rr[end]
    rr[1] = rr[end]
end

f = open("jenky.dat", "w")
for i in 1:length(xx)
    write(f, "$(xx[i]),  $(rr[i])\n")
end
close(f)
