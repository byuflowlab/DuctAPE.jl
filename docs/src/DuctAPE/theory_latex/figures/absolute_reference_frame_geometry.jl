#=

Attempt to make a Bspline parameterization for duct cross section

=#

include("../../../../../Code_Development/DuctTAPE/src/preliminary_design/1DModel_B.jl")
using Splines
using FLOWMath
const fm = FLOWMath

"""
cosine spacing
"""
function cosine_spacing(N)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
cosine spacing, but also scales and transforms
"""
function scaled_cosine_spacing(N, scale, translate; mypi=pi)
    return translate .+ scale * [0.5 * (1 - cos(mypi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
nosecone stop and tail cone start are ratio of internal length
tail cone tip radius is ratio of hub radius
internal length is from leading edge to duct trailing edge, assuming hub and duct LE are at same spot, and hub extends back further
hub chord is ratio to duct chord
nose_tip_cpx and tail_tip_cpx are also ratio of internal length
"""
function center_body_geom(
    Rhub,
    duct_chord;
    nosecone_start=0.125,
    nosecone_stop=0.35,
    tailcone_start=0.65,
    nose_tip_cpx=0.1,
    tail_tip_cpx=0.85,
    cb_te_radius=0.1,
    hub_chord=1.25,
    N=60,
)

    # - define flat portion first - #
    flatx = range(nosecone_stop * duct_chord, tailcone_start * duct_chord, 10)
    flatr = ones(10) * Rhub

    # - Nose Cone quadratic spline - #
    # knot vector
    knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # Nose control points
    #nose cpx2 is based on nose tip angle
    nosecone_cps = [
        [nosecone_start * duct_chord, 0.0],
        [(2.0 * nosecone_start + nosecone_stop) * duct_chord / 3.0, flatr[1]],
        [flatx[1], flatr[1]],
    ]

    # spline
    nose_spline = Splines.NURBS(2, knots, ones(length(nosecone_cps)), nosecone_cps)

    # - tail cone quadratic spline - #
    # tail control points
    tail_cps = [
        [flatx[end], flatr[end]],
        [tail_tip_cpx * duct_chord, flatr[end]],
        [hub_chord * duct_chord, cb_te_radius * Rhub],
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
    scale = tail_cps[end][1] - nosecone_cps[1][1]
    translate = nosecone_start * duct_chord
    cbx = scaled_cosine_spacing(N, scale, translate)
    cbsp = fm.Akima(xs, rs)
    cbr = cbsp(cbx)

    # return cbx, cbr, nosecone_cps, tail_cps
    return cbx, cbr, cbsp, nosecone_cps, tail_cps
end

"""
nacellevar is ratio relative to Rtip
"""
function duct_geom(
    duct_le_radius,
    duct_te_radius,
    Rtip,
    chord,
    cbspline;
    nosecone_stop=0.35,
    tailcone_start=0.65,
    nacellevar=1.0,
    N=60,
)

    # - define flat portion first - #
    flatx = range(nosecone_stop * chord, tailcone_start * chord, 10)
    flatr = ones(10) * Rtip

    # - Get TE radial control point - #
    # find hub radius at the exit plane
    hubr = cbspline(chord)
    # solve for duct radial point based on annulus area
    # a = pi * (ductr^2 - hubr^2)
    duct_te_radius = sqrt(exit_area / pi + hubr^2)

    # - Define outlet spline - #
    # knot vector
    knots2 = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # control points
    outlet_cps = [
        [flatx[end], flatr[end]],
        [chord * (1.0 + tailcone_start) / 2.0, flatr[end]],
        [chord, duct_te_radius],
    ]

    # spline
    outlet_spline = Splines.NURBS(2, knots2, ones(length(outlet_cps)), outlet_cps)

    # - Define inlet spline - #
    # control points
    inlet_cps = [[0.0, duct_le_radius], [0.0, flatr[end]], [flatx[1], flatr[1]]]

    # spline
    inlet_spline = Splines.NURBS(2, knots2, ones(length(inlet_cps)), inlet_cps)

    # - Define nacelle spline - #
    # knot vector
    knots3 = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    # control points
    nacelle_cps = [
        [0.0, duct_le_radius],
        [0.0, duct_le_radius + nacellevar],
        [chord / 2.0, duct_le_radius + nacellevar],
        [chord, duct_te_radius],
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

Rtip = 2.0 # inches
Rhub = 0.25 * Rtip # inches
# Rtip = 5.0 # inches
# Rhub = 0.25 * Rtip / 2.0 # inches
chord = 2 * Rtip
xrotor = Rtip
nosecone_start = 0.25
nosecone_stop = 0.45
tailcone_start = 0.65
nose_tip_cpx = 0.1
tail_tip_cpx = 1.0
cb_te_radius = 0.0
# cb_te_radius = 1.0

cbx, cbr, cbspline, nosecone_cps, tailcone_cps = center_body_geom(
    Rhub,
    chord;
    nosecone_start=nosecone_start,
    nosecone_stop=nosecone_stop,
    tailcone_start=tailcone_start,
    nose_tip_cpx=nose_tip_cpx,
    tail_tip_cpx=tail_tip_cpx,
    cb_te_radius=cb_te_radius,
    hub_chord=1.25,
    N=20,
)

f = open("absolute-frame-hub.dat", "w")
for (x, r) in zip(eachrow(cbx), eachrow(cbr))
    write(f, "$(x[1]) $(r[1])\n")
end
close(f)

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
facearea = pi * Rtip^2
capturearea = debug.Vface * facearea / Vinf
captureradius = sqrt(capturearea / pi)
exit_area /= 0.0254^2 # convert from m^2 to in^2
duct_te_radius = sqrt(exit_area / pi)

# inlet_area = pi * 5.875^2
inlet_area = pi * captureradius^2
duct_le_radius = captureradius

nacellevar = 0.75

cx, cr, nx, nr, cpi, cpo, cpn = duct_geom(
    duct_le_radius,
    duct_te_radius,
    Rtip,
    chord,
    cbspline;
    nosecone_stop=nosecone_stop,
    tailcone_start=tailcone_start,
    nacellevar=nacellevar,
)

ductx = [reverse(nx); cx[2:end]] ./ chord
ductr = ([reverse(nr); cr[2:end]] .- nr[1]) ./ chord

f = open("absolute-frame-duct.dat", "w")
for (x, r) in zip(eachrow([reverse(nx); cx[2:end]]), eachrow([reverse(nr); cr[2:end]]))
    write(f, "$(x[1]) $(r[1])\n")
end
close(f)
