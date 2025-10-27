using Splines
using FLOWMath
const fm = FLOWMath
using DuctAPE
const dt = DuctAPE
using Plots

# generate hub coordinates
cbchord = 4
Ncb = 101

# knot vector
knots3 = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
knotscb = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]

# control points
cb_cps = [
    [0.0, 0.0], [0.125, 0.5], [cbchord / 2, 0.75], [cbchord * 0.8, 0.25], [cbchord, 0.25]
]

# spline
cb_spline = Splines.NURBS(3, knotscb, ones(length(cb_cps)), cb_cps)

# - Get cb Points - #
u = range(0.0, 1.0; length=Ncb + 1)
n = length(cb_spline.ctrlpts) - 1
Cw_cb = [[0.0; 0.0] for i in 1:length(u)]

#loop through parametric points to get curve points
for i in 1:length(u)
    Cw_cb[i] = Splines.curvepoint(cb_spline, u[i])
end

# - Respline the whole geometry to make things smooth - #
# assemble array
cb_xs = getindex.(Cw_cb, 1)
cb_rs = getindex.(Cw_cb, 2)

# generate duct inner coordinates
duct_le_radius = 3.75
duct_te_radius = 3.0
chord = 2.5
N = round(Int, Ncb * chord / cbchord)

# - Define nacelle spline - #

# control points
casing_cps = [
    [0.0, duct_le_radius],
    [0.0, duct_le_radius - 0.25],
    [chord / 9.0, duct_te_radius],
    [chord, duct_te_radius],
]

# spline
casing_spline = Splines.NURBS(3, knots3, ones(length(casing_cps)), casing_cps)

# - Get casing Points - #
u = range(0.0, 1.0; length=N + 1)
n = length(casing_spline.ctrlpts) - 1
Cw_casing = [[0.0; 0.0] for i in 1:length(u)]

#loop through parametric points to get curve points
for i in 1:length(u)
    Cw_casing[i] = Splines.curvepoint(casing_spline, u[i])
end

# - Respline the whole geometry to make things smooth - #
# assemble array
casing_xs = getindex.(Cw_casing, 1)
casing_rs = getindex.(Cw_casing, 2)

# control points
nacelle_cps = [
    [0.0, duct_le_radius],
    [0.0, duct_le_radius + 0.35],
    [chord / 1.75, duct_te_radius + 0.25],
    [chord, duct_te_radius],
]

# spline
nacelle_spline = Splines.NURBS(3, knots3, ones(length(nacelle_cps)), nacelle_cps)

# - Get nacelle Points - #
u = range(0.0, 1.0; length=N + 1)
n = length(nacelle_spline.ctrlpts) - 1
Cw_nacelle = [[0.0; 0.0] for i in 1:length(u)]

#loop through parametric points to get curve points
for i in 1:length(u)
    Cw_nacelle[i] = Splines.curvepoint(nacelle_spline, u[i])
end

# - Respline the whole geometry to make things smooth - #
# assemble array
nacelle_xs = getindex.(Cw_nacelle, 1)
nacelle_rs = getindex.(Cw_nacelle, 2)

# extend coordinates
Nex = 50
lower_wall_x = [
    range(-1.5, 0.0, Nex)
    cb_xs[2:(end - 1)]
    range(cbchord, cbchord + 1.5, Nex)
]
lower_wall_r = [cb_rs[1] * ones(Nex); cb_rs[2:(end - 1)]; cb_rs[end] * ones(Nex)]
upper_wall_x = [
    range(-1.5, 0.0, Nex)
    casing_xs[2:(end - 1)]
    range(chord, cbchord + 1.5, Nex + Ncb - N)
]
upper_wall_r = [
    casing_rs[1] * ones(Nex)
    casing_rs[2:(end - 1)]
    casing_rs[end] * ones(Ncb - N + Nex)
]

# initialize grid
nr = 101
nx = 151
grid = zeros(2, nx, nr)

# interpolate coordinates
xglob = range(-1.5, cbchord + 1.5, nx)
lrs = fm.akima(lower_wall_x, lower_wall_r, xglob)
urs = fm.akima(upper_wall_x, upper_wall_r, xglob)
grid[1, :, :] .= xglob
grid[2, :, :] .= reduce(hcat, range(lrs, urs, nr))
# grid[1, :, :] .= reduce(hcat, range(lower_wall_x, upper_wall_x, nr))
# grid[2, :, :] .= reduce(hcat, range(lower_wall_r, upper_wall_r, nr))

# relax grid
dt.relax_grid!(grid; max_iterations=100, tol=1e-9, verbose=false)

# save duct and centerbody geometry
f = open("grid-coord-hub.dat", "w")
for (x, r) in zip(cb_xs, cb_rs)
    write(f, "$(x) $(r)\n")
end
write(f, "$(cb_xs[end]) 0.0")
close(f)

f = open("grid-coord-duct.dat", "w")
for (x, r) in
    zip([reverse(casing_xs); nacelle_xs[2:end]], [reverse(casing_rs); nacelle_rs[2:end]])
    write(f, "$(x) $(r)\n")
end
close(f)

# save grid vertical lines
xixcoarse = grid[1, 1:20:end, :]'
xircoarse = grid[2, 1:20:end, :]'

g = open("etaarrow.dat", "w")
for (i, (xix, xir)) in enumerate(zip(eachcol(xixcoarse), eachcol(xircoarse)))
    f = open("xiline$(i).dat", "w")
    for (x, r) in zip(xix, xir)
        write(f, "$(x) $(r)\n")
        if i == 2 && r >= etarcoarse[1,2] && r <= sum(etarcoarse[1, 3:4]) / 2.0
            write(g, "$(x) $(r)\n")
        end
    end
    close(f)
end
close(g)

# save grid horizontal lines
etaxcoarse = grid[1, :, 1:20:end]
etarcoarse = grid[2, :, 1:20:end]

g = open("xiarrow.dat", "w")
for (i, (etax, etar)) in enumerate(zip(eachcol(etaxcoarse), eachcol(etarcoarse)))
    f = open("etaline$(i).dat", "w")
    for (x, r) in zip(etax, etar)
        write(f, "$(x) $(r)\n")
        if i == 2 && x >= xixcoarse[1, 2] && x <= sum(xixcoarse[1, 3:4]) / 2.0
            write(g, "$(x) $(r)\n")
        end
    end
    close(f)
end
close(g)

# save eta arrow coordinates

# save xi arrow coordinates

plot(cb_xs, cb_rs; aspectratio=1, label="")
plot!(nacelle_xs, nacelle_rs; label="")
plot!(casing_xs, casing_rs; label="")
plot!(lower_wall_x, lower_wall_r; label="")
plot!(upper_wall_x, upper_wall_r; label="")
# plot!(grid[1, :, :], grid[2, :, :]; color=1, label="")
# plot!(grid[1, :, :]', grid[2, :, :]'; color=2, label="")

plot!(xixcoarse, xircoarse; linewidth=3, color=2, label="")

plot!(etaxcoarse, etarcoarse; linewidth=3, color=1, label="")
