#=

Attempt to make a Bspline parameterization for duct cross section

=#

include("../visualize/plots_default_new.jl")
include("../src/preliminary_design/1DModel_B.jl")
using Splines
using FLOWMath
const fm = FLOWMath

#=
Redo TODOs:

variables include Rhub, Rtip, and exit area, also maybe rotor axial position (in the form of inlet length as a percent of total enclosed length?)
constant parameters will be overall lengths, inlet area, etc.
=#

"""
cosine spacing
"""
function cosine_spacing(N)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
nosecone stop and tail cone start are ratio of internal length
tail cone tip radius is ratio of hub radius
internal length is from leading edge to duct trailing edge, assuming hub and duct LE are at same spot, and hub extends back further
hub chord is ratio to duct chord
nose_tip_cpx and tail_tip_cpx are also ratio of internal length
"""
function centerbody_geom(
    Rhub,
    internal_length;
    nosecone_stop=0.35,
    tailcone_start=0.65,
    nose_tip_cpx=0.0,
    tail_tip_cpx=0.85,
    tail_tip_radius=0.1,
    hub_chord=1.25,
    N=60,
)

    # - define flat portion first - #
    flatx = range(nosecone_stop * internal_length, tailcone_start * internal_length, 10)
    flatr = ones(10) * Rhub

    # - Nose Cone quadratic spline - #
    # knot vector
    knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # Nose control points
    #nose cpx2 is based on nose tip angle
    nose_cps = [
        [0.0, 0.0], [nose_tip_cpx * internal_length, flatr[1]], [flatx[1], flatr[1]]
    ]

    # spline
    nose_spline = Splines.NURBS(2, knots, ones(length(nose_cps)), nose_cps)

    # - tail cone quadratic spline - #
    # tail control points
    tail_cps = [
        [flatx[end], flatr[end]],
        [tail_tip_cpx * internal_length, flatr[end]],
        [hub_chord * internal_length, tail_tip_radius * Rhub],
    ]

    # spline
    tail_spline = Splines.NURBS(2, knots, ones(length(tail_cps)), tail_cps)

    # - Get nose and tail points - #
    u = range(0.0, 1.0; length=N + 1)
    n = length(nose_spline.ctrlpts) - 1
    Cw_nose = [[0.0; 0.0] for i in 1:length(u)]
    Cw_tail = [[0.0; 0.0] for i in 1:length(u)]

    #loop through parametric points to get curve points
    for i in 1:length(u)
        Cw_nose[i] = Splines.curvepoint(nose_spline, u[i])
        Cw_tail[i] = Splines.curvepoint(tail_spline, u[i])
    end

    # - Respline the whole geometry to make things smooth - #
    # assemble array
    xs = [getindex.(Cw_nose, 1); flatx[2:(end - 1)]; getindex.(Cw_tail, 1)]
    rs = [getindex.(Cw_nose, 2); flatr[2:(end - 1)]; getindex.(Cw_tail, 2)]

    # re-spline with akima spline and cosine spacing
    cbx = cosine_spacing(N) * internal_length * hub_chord
    cbsp = fm.Akima(xs, rs)
    cbr = cbsp(cbx)

    # return cbx, cbr, nose_cps, tail_cps
    return cbx, cbr, cbsp, nose_cps, tail_cps
end

"""
nacell_var is ratio relative to Rtip
"""
function duct_geom(
    inlet_area,
    exit_area,
    Rtip,
    chord,
    cbspline;
    nosecone_stop=0.35,
    tailcone_start=0.65,
    nacell_var=1.1,
    N=60,
)
    le_cpr = sqrt(inlet_area / pi)
    # - define flat portion first - #
    flatx = range(nosecone_stop * chord, tailcone_start * chord, 10)
    flatr = ones(10) * Rtip

    # - Get TE radial control point - #
    # find hub radius at the exit plane
    hubr = cbspline(chord)
    # solve for duct radial point based on annulus area
    # a = pi * (ductr^2 - hubr^2)
    ductter = sqrt(exit_area / pi + hubr^2)

    # - Define outlet spline - #
    # knot vector
    knots2 = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # control points
    outlet_cps = [
        [flatx[end], flatr[end]],
        [chord * (1.0 + tailcone_start) / 2.0, flatr[end]],
        [chord, ductter],
    ]

    # spline
    outlet_spline = Splines.NURBS(2, knots2, ones(length(outlet_cps)), outlet_cps)

    # - Define inlet spline - #
    # control points
    inlet_cps = [[0.0, le_cpr], [0.0, flatr[end]], [flatx[1], flatr[1]]]

    # spline
    inlet_spline = Splines.NURBS(2, knots2, ones(length(inlet_cps)), inlet_cps)

    # - Define nacelle spline - #
    # knot vector
    knots3 = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    # control points
    nacelle_cps = [
        [0.0, le_cpr],
        [0.0, le_cpr + 1.25 * nacell_var],
        [chord / 2.0, le_cpr + nacell_var],
        [chord, ductter],
    ]

    # spline
    nacelle_spline = Splines.NURBS(3, knots3, ones(length(nacelle_cps)), nacelle_cps)

    # - Get Casing Points - #
    u = range(0.0, 1.0; length=N + 1)
    n = length(inlet_spline.ctrlpts) - 1
    Cw_inlet = [[0.0; 0.0] for i in 1:length(u)]
    Cw_outlet = [[0.0; 0.0] for i in 1:length(u)]

    #loop through parametric points to get curve points
    for i in 1:length(u)
        Cw_inlet[i] = Splines.curvepoint(inlet_spline, u[i])
        Cw_outlet[i] = Splines.curvepoint(outlet_spline, u[i])
    end

    # - Respline the whole geometry to make things smooth - #
    # assemble array
    casing_xs = [getindex.(Cw_inlet, 1); flatx[2:(end - 1)]; getindex.(Cw_outlet, 1)]
    casing_rs = [getindex.(Cw_inlet, 2); flatr[2:(end - 1)]; getindex.(Cw_outlet, 2)]

    # re-spline with akima spline and cosine spacing
    cx = cosine_spacing(N) * chord
    csp = fm.Akima(casing_xs, casing_rs)
    cr = csp(cx)

    # - Get Nacelle Points - #
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

    # re-spline with akima spline and cosine spacing
    nx = cosine_spacing(N) * chord
    nsp = fm.Akima(nacelle_xs, nacelle_rs)
    nr = nsp(nx)

    return cx, cr, nx, nr, inlet_cps, outlet_cps, nacelle_cps
    # return dx, dr, dsp, cps
end

Rtip = 5.0 # inches
Rhub = 1.25 / 2.0 # inches
chord = 2 * Rtip
xrotor = Rtip
nosecone_stop = 0.45
tailcone_start = 0.65
nose_tip_cpx = 0.1
tail_tip_cpx = 1.0
tail_tip_radius = 1.0

cbx, cbr, cbspline, ncp, tcp = centerbody_geom(
    Rhub,
    chord;
    nosecone_stop=nosecone_stop,
    tailcone_start=tailcone_start,
    nose_tip_cpx=nose_tip_cpx,
    tail_tip_cpx=tail_tip_cpx,
    tail_tip_radius=tail_tip_radius,
    hub_chord=1.25,
    N=20,
)

# - APC 10x7 DATA at 4011 RPM - #
# J       CT       CP       eta
#
rotor_data = [
    0.144 0.1389 0.0726 0.276
    0.180 0.1339 0.0719 0.335
    0.214 0.1289 0.0710 0.389
    0.251 0.1229 0.0699 0.442
    0.287 0.1174 0.0686 0.491
    0.327 0.1102 0.0666 0.541
    0.361 0.1039 0.0649 0.578
    0.390 0.0984 0.0632 0.606
    0.437 0.0903 0.0610 0.648
    0.468 0.0849 0.0591 0.672
    0.501 0.0789 0.0571 0.692
    0.539 0.0724 0.0546 0.714
    0.568 0.0659 0.0520 0.720
    0.611 0.0576 0.0487 0.723
    0.647 0.0498 0.0453 0.711
    0.674 0.0438 0.0427 0.691
    0.718 0.0326 0.0374 0.627
]

Jfine = range(rotor_data[1, 1], rotor_data[end, 1], 50)
etafine = fm.akima(rotor_data[:, 1], rotor_data[:, end], Jfine)
ctfine = fm.akima(rotor_data[:, 1], rotor_data[:, 2], Jfine)
Jop = Jfine[findmax(etafine)[2]]
ctop = ctfine[findmax(etafine)[2]]
RPM = 4000 #change
Vinf = calculate_vstar(Rtip * 0.0254, RPM, 0.5)

exit_area, debug = size_exit(Rtip * 0.0254, RPM, Jop, ctop, Vinf; rho=1.225)
exit_area /= 0.0254^2 # convert from m^2 to in^2

inlet_area = pi * 5.875^2

cx, cr, nx, nr, cpi, cpo, cpn = duct_geom(
    inlet_area,
    exit_area,
    Rtip,
    chord,
    cbspline;
    nosecone_stop=nosecone_stop,
    tailcone_start=tailcone_start,
    nacell_var=1.25,
)

ductx = [reverse(nx); cx[2:end]] ./ chord
ductr = ([reverse(nr); cr[2:end]] .- nr[1]) ./ chord

f = open("d.dat", "w")
write(f, "duct\n")
for (x, r) in zip(eachrow(ductx), eachrow(ductr))
    write(f, "$(x[1]) $(r[1])\n")
end
close(f)

plot(; xlabel="x (in)", ylabel="r (in)", aspectratio=1, label="Center Body Geometry")
plot!(cbx, cbr; color=myblue[1], label="Center Body Geometry")
plot!(
    getindex.(ncp, 1),
    getindex.(ncp, 2);
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=myred[2],
    label="Nose Cone Control Points",
)
plot!(
    getindex.(tcp, 1),
    getindex.(tcp, 2);
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=myred[3],
    label="Tail Cone Control Points",
)
savefig("hubtest.pdf")

plot!(cx, cr; color=myblue[2], label="Duct Casing")
plot!(
    getindex.(cpi, 1),
    getindex.(cpi, 2);
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=mygreen[2],
    label="Inlet Control Points",
)
plot!(
    getindex.(cpo, 1),
    getindex.(cpo, 2);
    seriestype=:scatter,
    markersize=3,
    markershape=:square,
    color=mygreen[3],
    label="Outlet Control Points",
)

plot!(nx, nr; color=myblue[3], label="Duct Nacelle")
plot!(
    getindex.(cpn, 1),
    getindex.(cpn, 2);
    seriestype=:scatter,
    markersize=3,
    markershape=:diamond,
    color=mygray[1],
    label="Nacelle Control Points",
)

plot!(
    xrotor * ones(2), [maximum(nr), 0.0]; linestyle=:dash, color=:black, label="Rotor Plane"
)
savefig("ducttest.pdf")

#=
TODOs for exit area sizing

go get APC 10x7 data and get all the inputs for the area calcuation
get area.
=#

# script for duct geometry generation and plotting
# area will determine TE control point locations
# probably use quadratic from flat to TE
# probably use quadratic from flat to LE
# probably use cubic on top
# probably keep it a bit thin if possible
#
# after you get the geometry, run it through xfoil as a check on separation.
# probably want to run it at a slight negative angle of attack to simulate the affects of the rotor
#
#
#
#
#

#function generate_duct_section(
#    Rin,
#    Rd,
#    Rout,
#    Lin,
#    Ltot,
#    r_le,
#    wedge_angle,
#    te_camber_angle,
#    max_thick;
#    N=161,
#    disk_ratio=0.125,
#)
#    boattailangle = pi * wedge_angle / 360 #boattail angle is half of wedge angle

#    te_camber_angle *= pi / 180

#    #--Initialize knot vector
#    knots_inner = [
#        0.0,
#        0.0,
#        0.0,
#        Lin * (1.0 - disk_ratio) / Ltot,
#        Lin * (1.0 - disk_ratio) / Ltot,
#        Lin * (1.0 + disk_ratio) / Ltot,
#        Lin * (1.0 + disk_ratio) / Ltot,
#        1.0,
#        1.0,
#        1.0,
#    ]

#    knots_outer = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

#    #calculate unweighted controlpoints
#    Pouter = [
#        [0.0, Rin], [0.0, Rin + (Rin - Rd)], [Lin, Rd + 2.0 * max_thick], [Ltot, Rout]
#    ]

#    Pinner = [
#        [0.0, Rin], #sqrt(2.0 * leradius) / 3.0],
#        [0.0, Rd],
#        [Lin * (1.0 - disk_ratio), Rd],
#        [Lin * (1.0 - disk_ratio), Rd],
#        [Lin * (1.0 + disk_ratio), Rd],
#        # [Lin * 1.1, Rd],
#        [(Ltot + Lin) / 2.0, Rd],
#        [Ltot, Rout],
#    ]

#    inner_spline = Splines.NURBS(2, knots_inner, ones(length(Pinner)), Pinner)
#    outer_spline = Splines.NURBS(3, knots_outer, ones(length(Pouter)), Pouter)

#    u = range(0.0, 1.0; length=N + 1)
#    n = length(inner_spline.ctrlpts) - 1
#    Cw_inner = [[0.0; 0.0] for i in 1:length(u)]
#    Cw_outer = [[0.0; 0.0] for i in 1:length(u)]

#    #loop through parametric points to get curve points
#    for i in 1:length(u)
#        Cw_inner[i] = Splines.curvepoint(inner_spline, u[i])
#        Cw_outer[i] = Splines.curvepoint(outer_spline, u[i])
#    end

#    return Pinner, Pouter, Cw_inner, Cw_outer
#end

#Rin = 2.75
#Rd = 2.5
#Rout = 2.35
#Ltot = Rd * 2.0
#Lin = 0.4 * Ltot
#r_le = 0.25
#wedge_angle = 20.0
#te_camber_angle = 5.0
#max_thick = 0.5
#N = 161
#disk_ratio = 0.125

#pin, pout, xrin, xrout = generate_duct_section(
#    Rin,
#    Rd,
#    Rout,
#    Lin,
#    Ltot,
#    r_le,
#    wedge_angle,
#    te_camber_angle,
#    max_thick;
#    N=N,
#    disk_ratio=disk_ratio,
#)

#plot(; aspectratio=:equal)
#plot!(getindex.(pin, 1), getindex.(pin, 2); seriestype=:scatter)
#plot!(getindex.(pout, 1), getindex.(pout, 2); seriestype=:scatter)
#plot!(getindex.(xrout, 1), getindex.(xrout, 2))
#plot!(getindex.(xrin, 1), getindex.(xrin, 2))

#using FLOWMath
#using FLOWFoil
#t = Ltot .* FLOWFoil.cosine_spacing(161)
#osp = akima(getindex.(xrout, 1), getindex.(xrout, 2), t)
#isp = akima(getindex.(xrin, 1), getindex.(xrin, 2), t)

#xx = [reverse(t); t[2:end]]
#rr = [reverse(isp); osp[2:end]]
#if rr[1] != rr[end]
#    rr[1] = rr[end]
#end

#f = open("jenky.dat", "w")
#for i in 1:length(xx)
#    write(f, "$(xx[i]),  $(rr[i])\n")
#end
#close(f)
