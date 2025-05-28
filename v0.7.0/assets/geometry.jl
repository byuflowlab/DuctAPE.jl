using NURBS
using StaticArrays
using FLOWMath
const fm = FLOWMath

"""
cosine spacing, but also scales and transforms
"""
function scaled_cosine_spacing(N, scale, translate; mypi=pi)
    return translate .+ scale * [0.5 * (1 - cos(mypi * (i - 1) / (N - 1))) for i in 1:N]
end

"""

# Arguments
- `Rhub::Float` : maximum hub radius
- `duct_chord::Float` : dimensional duct (annular airfoil) chord

# Keyword Arguments:
- `cb_nc_le::Float=0.125` : centerbody nosecone leading edge, non-dimensionalized by duct_chord
- `cb_nc_stop::Float=0.35` : centerbody nosecone stop point, non-dimensionalized by duct_chord
- `cb_tc_start::Float=1.0` : centerbody tailcone start point, non-dimensionalized by duct_chord
- `cb_ncp_z::Float=0.25` : centerbody nosecone spline 2nd control point axial location, non-dimensionalized by duct_chord
- `cb_tcp_z::Float=0.5` : centerbody tailcone spline 2nd control point axial location, non-dimensionalized by duct_chord
- `cb_te_r::Float=1.0` : centerbody trailing edge radial position, non-dimensionalized by Rhub
- `cb_te_z::Float=1.1` : centerbody trailing edge axial position, non-dimensionalized by duct_chord
- `N::Float=60` : number of points to use in intermediate splines
"""
function centerbody_geom(
    Rhub,
    duct_chord;
    cb_nc_le=0.0625,
    cb_nc_stop=0.35,
    cb_tc_start=1.0,
    cb_ncp_z=0.35,
    cb_tcp_z=0.5,
    cb_te_r=0.9,
    cb_te_z=1.1,
    N=60,
    fmspline=(x, y) -> FLOWMath.Akima(x, y, 1e-100, 1e-200),
    smooth=false,
)
    TF = promote_type(
        typeof(Rhub),
        typeof(duct_chord),
        typeof(cb_nc_le),
        typeof(cb_nc_stop),
        typeof(cb_tc_start),
        typeof(cb_ncp_z),
        typeof(cb_tcp_z),
        typeof(cb_te_r),
        typeof(cb_te_z),
    )

    # - Dimensionalize Everything - #
    cb_nc_le *= duct_chord
    cb_nc_stop *= duct_chord
    cb_tc_start *= duct_chord
    cb_te_r *= Rhub
    cb_te_z *= duct_chord

    # - define flat portion first - #
    flatx = range(cb_nc_stop, cb_tc_start, 10)
    flatr = ones(TF, 10) * Rhub

    # - Nose Cone quadratic spline - #
    # knot vector
    knots = TF[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # Nose control points
    #nose cpx2 is based on nose tip angle
    nosecone_cps = [
        SVector(cb_nc_le, 0.0, 0.0),
        SVector(cb_nc_le + cb_ncp_z * (cb_nc_stop - cb_nc_le), flatr[1], 0.0),
        SVector(flatx[1], flatr[1], 0.0),
    ]

    # spline
    nose_spline = NURBS.NURBScurve(
        NURBS.NURB(2, knots, ones(TF, size(nosecone_cps, 1))), nosecone_cps
    )

    # - tail cone quadratic spline - #
    # tail control points
    tail_cps = [
        SVector(flatx[end], flatr[end], 0.0),
        SVector(flatx[end] + cb_tcp_z * (cb_te_z - cb_tc_start), flatr[end], 0.0),
        SVector(cb_te_z, cb_te_r, 0.0),
    ]

    # spline
    tail_spline = NURBS.NURBScurve(
        NURBS.NURB(2, knots, ones(TF, size(tail_cps, 1))), tail_cps
    )

    # - Get nose and tail points - #
    u = collect(TF, range(0.0, 1.0; length=N + 1))
    Cw_nose = nose_spline(u)
    Cw_tail = tail_spline(u)

    # - Respline the whole geometry to make things smooth - #
    # assemble array
    zs = [getindex.(Cw_nose, 1); flatx[2:(end - 1)]; getindex.(Cw_tail, 1)]
    rs = [getindex.(Cw_nose, 2); flatr[2:(end - 1)]; getindex.(Cw_tail, 2)]

    if smooth
        # re-spline with akima spline and cosine spacing
        scale = tail_cps[end][1] - nosecone_cps[1][1]
        translate = cb_nc_le
        cbz = scaled_cosine_spacing(N, scale, translate)
        cbsp = fmspline(zs, rs)
        cbr = cbsp(cbz)

        return cbz, cbr, cbsp, nosecone_cps, tail_cps
    else
        return zs, rs, nosecone_cps, tail_cps
    end
end

function duct_geom(
    Rtip, duct_chord, duct_le_radius, te_camber_angle, wedge_angle; duct_alpha=0, N=60
)
    TF = promote_type(
        eltype(Rtip),
        eltype(duct_chord),
        eltype(duct_le_radius),
        eltype(te_camber_angle),
        eltype(wedge_angle),
    )

    boat_tail_angle = wedge_angle * pi / 360.0
    te_camber_angle *= pi / 180.0

    knots3 = TF[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    nacelle_cps = [
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, sqrt(duct_le_radius) / 3.0, 0.0),
        SVector(1.0 / 3.0, 2.0 * tan(te_camber_angle + boat_tail_angle) / 3.0, 0.0),
        SVector(1.0, 0.0, 0.0),
    ]

    casing_cps = [
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, -sqrt(duct_le_radius) / 3.0, 0.0),
        SVector(1.0 / 3.0, 2.0 * tan(te_camber_angle - boat_tail_angle) / 3.0, 0.0),
        SVector(1.0, 0.0, 0.0),
    ]

    # spline
    nacelle_spline = NURBS.NURBScurve(
        NURB(3, knots3, ones(TF, size(nacelle_cps, 1))), nacelle_cps
    )

    casing_spline = NURBS.NURBScurve(
        NURB(3, knots3, ones(TF, size(casing_cps, 1))), casing_cps
    )

    u = collect(TF, range(0.0, 1.0; length=N + 1))

    Cw_nacelle = nacelle_spline(u) * duct_chord
    Cw_casing = casing_spline(u) * duct_chord

    # get rotation matrix
    R = [
        cosd(-duct_alpha) -sind(-duct_alpha) 0.0
        sind(-duct_alpha) cosd(-duct_alpha) 0.0
        0.0 0.0 1.0
    ]
    # rotate and translate
    rotated_nacelle = [R * c for c in Cw_nacelle]
    rotated_casing = [R * c for c in Cw_casing]

    nacelle_zs = getindex.(rotated_nacelle, 1)
    nacelle_rs = getindex.(rotated_nacelle, 2)

    casing_zs = getindex.(rotated_casing, 1)
    casing_rs = getindex.(rotated_casing, 2)

    return nacelle_zs, nacelle_rs, casing_zs, casing_rs, nacelle_cps, casing_cps
end

function duct_geom(
    Rtip,
    duct_chord;
    duct_le_r=1.1,
    duct_te_r=1.0,
    duct_le_z=0.0,
    duct_inlet_stop=1.0 / 3.0,
    duct_outlet_start=2.0 / 3.0,
    duct_outlet_cpz=0.5,
    nacelle_cpz=0.5,
    nacelle_thickness_var=2.0,
    N=60,
    fmspline=(x, y) -> FLOWMath.Akima(x, y, 1e-100, 1e-200),
    smooth=false,
)
    TF = promote_type(
        typeof(Rtip),
        typeof(duct_chord),
        typeof(duct_le_z),
        typeof(duct_le_r),
        typeof(duct_te_r),
        typeof(duct_inlet_stop),
        typeof(duct_outlet_start),
    )

    # - Dimensionalize Everything - #
    duct_le_r *= Rtip
    duct_te_r *= Rtip
    duct_le_z *= duct_chord
    duct_inlet_stop *= duct_chord
    duct_outlet_start *= duct_chord
    nacelle_cpz *= duct_chord
    nacelle_cpr = (duct_le_r - Rtip) * nacelle_thickness_var + duct_le_r
    duct_outlet_cpz *= (duct_chord - duct_outlet_start)
    duct_outlet_cpz += duct_outlet_start
    delta_r_te = duct_te_r - Rtip

    # - define flat portion first - #
    flatx = range(duct_inlet_stop, duct_outlet_start, 10)
    flatr = ones(TF, 10) * Rtip

    # - Define outlet spline - #
    # knot vector
    knots2 = TF[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # control points
    outlet_cps = [
        SVector(flatx[end], flatr[end], 0.0),
        SVector(duct_outlet_cpz, flatr[end], 0.0),
        SVector(duct_chord, duct_te_r, 0.0),
    ]

    # spline
    outlet_spline = NURBS.NURBScurve(
        NURB(2, knots2, ones(TF, size(outlet_cps, 1))), outlet_cps
    )

    # - Define inlet spline - #
    # control points
    inlet_cps = [
        SVector(0.0, duct_le_r, 0.0),
        SVector(0.0, flatr[end], 0.0),
        SVector(flatx[1], flatr[1], 0.0),
    ]

    # spline
    inlet_spline = NURBS.NURBScurve(
        NURB(2, knots2, ones(TF, size(inlet_cps, 1))), inlet_cps
    )

    ## --- Nominal, does not preserve trailing edge angle --- ###
    # knot vector
    knots3 = TF[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    # control points
    nacelle_cps = [
        SVector(0.0, duct_le_r, 0.0),
        SVector(0.0, nacelle_cpr, 0.0),
        SVector(nacelle_cpz, nacelle_cpr, 0.0),
        SVector(duct_chord, duct_te_r, 0.0),
    ]

    # spline
    nacelle_spline = NURBS.NURBScurve(
        NURB(3, knots3, ones(TF, size(nacelle_cps, 1))), nacelle_cps
    )

    # - Get Casing Points - #
    u = collect(TF, range(0.0, 1.0; length=N + 1))
    Cw_inlet = inlet_spline(u)
    Cw_outlet = outlet_spline(u)

    # assemble array
    casing_zs = [getindex.(Cw_inlet, 1); flatx[2:(end - 1)]; getindex.(Cw_outlet, 1)]
    casing_rs = [getindex.(Cw_inlet, 2); flatr[2:(end - 1)]; getindex.(Cw_outlet, 2)]

    # - Get Nacelle Points - #
    u = collect(TF, range(0.0, 1.0; length=N + 1))
    Cw_nacelle = nacelle_spline(u)

    # assemble array
    nacelle_zs = getindex.(Cw_nacelle, 1)
    nacelle_rs = getindex.(Cw_nacelle, 2)

    if smooth
        # - Respline the whole geometry to make things smooth - #
        # re-spline with akima spline and cosine spacing
        cz = scaled_cosine_spacing(N, duct_chord, duct_le_z)
        csp = fmspline(casing_zs, casing_rs)
        cr = csp(cz)

        # re-spline with akima spline and cosine spacing
        nz = scaled_cosine_spacing(N, duct_chord, duct_le_z)
        nsp = fmspline(nacelle_zs, nacelle_rs)
        nr = nsp(nz)

        return cz, cr, nz, nr, inlet_cps, outlet_cps, nacelle_cps
    else
        return casing_zs,
        casing_rs, nacelle_zs, nacelle_rs, inlet_cps, outlet_cps,
        nacelle_cps
    end
end
"""


# Arguments:

- `Rtip::Float` : rotor tip radius
- `duct_chord::Float` : duct (annular airfoil) chord length
- `duct_le_r::Float` : duct leading edge radial position, non-dimensionalized by Rtip
- `duct_te_r::Float` : duct trailing edge radial position, non-dimensionalized by Rtip

# Keyword Arguments:

- `duct_le_z::Float=0.0` : duct leading edge axial position
- `duct_inlet_stop::Float=0.35` : duct inlet stop point (analogous to centerbody nosecone stop point, also non-dimensionalized by duct_chord)
- `duct_outlet_start::Float=0.65` : duct outlet start point (analogous to centerbody tailcone start point, also non-dimensionalized by duct_chord)
- `duct_outlet_cpz::Float=0.5` : duct casing outlet spline 2nd control point axial location, non-dimensionalized by outlet length
- `nacelle_cpz::Float=0.5` : duct nacelle (outer side) spline 3rd control point axial location, non-dimensionalized by duct_chord
- `nacelle_thickness_var::Float=3.0` : duct nacelle (outer side) 2nd and 3rd control point radial position, scaled relative to duct_le_r and non-dimensionalized by Rtip. In other words, setting this to 1.0 will give the same distance from the rotor tip to the duct leading edge as from the duct leading edge to the 2nd and 3rd nacelle control points.  Although this does control the nacelle thickness, it does not direclty set the max thickness point of the nacelle.
"""
function duct_geom(
    Rtip,
    duct_chord,
    duct_le_r,
    duct_le_radius;
    duct_te_r=1.0,
    wedge_angle=14.0,
    duct_le_z=0.0,
    duct_inlet_stop=1.0 / 3.0,
    duct_outlet_start=2.0 / 3.0,
    duct_outlet_cpz=0.5,
    nacelle_cpz=0.5,
    nacelle_thickness_var=2.0,
    N=60,
    fmspline=(x, y) -> FLOWMath.Akima(x, y, 1e-100, 1e-200),
    smooth=false,
)
    TF = promote_type(
        typeof(Rtip),
        typeof(duct_chord),
        typeof(duct_le_z),
        typeof(duct_le_r),
        typeof(duct_te_r),
        typeof(duct_inlet_stop),
        typeof(duct_outlet_start),
    )

    # - Dimensionalize Everything - #
    duct_le_r *= Rtip
    duct_te_r *= Rtip
    duct_le_z *= duct_chord
    duct_inlet_stop *= duct_chord
    duct_outlet_start *= duct_chord
    nacelle_cpz *= duct_chord
    nacelle_cpr = (duct_le_r - Rtip) * nacelle_thickness_var + duct_le_r
    duct_outlet_cpz *= (duct_chord - duct_outlet_start)
    duct_outlet_cpz += duct_outlet_start
    wedge_angle *= pi / 180.0
    delta_r_te = duct_te_r - Rtip

    # - define flat portion first - #
    flatx = range(duct_inlet_stop, duct_outlet_start, 10)
    flatr = ones(TF, 10) * Rtip

    # - Define outlet spline - #
    # knot vector
    knots2 = TF[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    # control points
    # if duct_te_r <= Rtip
    outlet_cps = [
        SVector(flatx[end], flatr[end], 0.0),
        SVector(duct_outlet_cpz, flatr[end], 0.0),
        SVector(duct_chord, duct_te_r, 0.0),
    ]
    # else
    #     ncp3z0 = duct_chord - (nacelle_cpr - Rtip) / tan(wedge_angle)
    #     upper_theta = atan((nacelle_cpr - duct_te_r), duct_chord - ncp3z0)
    #     duct_outlet_cpz = duct_chord - (duct_te_r - Rtip) / tan(wedge_angle - upper_theta)
    #     println("duct_outlet_cpz = ", duct_outlet_cpz)
    #     outlet_cps = [
    #         SVector(flatx[end], flatr[end], 0.0),
    #         SVector(duct_outlet_cpz, flatr[end], 0.0),
    #         SVector(duct_chord, duct_te_r, 0.0),
    #     ]
    # end

    # spline
    outlet_spline = NURBS.NURBScurve(
        NURB(2, knots2, ones(TF, size(outlet_cps, 1))), outlet_cps
    )

    # - Define inlet spline - #
    # control points
    inlet_cps = [
        SVector(0.0, duct_le_r, 0.0),
        SVector(0.0, flatr[end], 0.0),
        SVector(flatx[1], flatr[1], 0.0),
    ]

    # spline
    inlet_spline = NURBS.NURBScurve(
        NURB(2, knots2, ones(TF, size(inlet_cps, 1))), inlet_cps
    )

    # # - Define nacelle spline - #
    # ### --- NURBS with only 3 pts (see Rajnarayan paper) --- ###
    # lower_theta = atan((Rtip-duct_te_r)/(duct_chord - duct_outlet_cpz))
    # println("beta = ", (lower_theta+wedge_angle)*180/pi)
    # delta_r_te = (duct_te_r - duct_le_r)
    # knotsN = TF[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    # nacelle_cps = [
    #     SVector(0.0, duct_le_r, 0.0),
    #     SVector(0.0, duct_le_r + delta_r_te + tan(lower_theta+wedge_angle), 0.0),
    #     SVector(duct_chord, duct_te_r, 0.0),
    # ]
    # nacelle_weights = [1.0; sqrt(duct_le_radius / 2.0) / (delta_r_te + tan(boat_tail_angle)); 1.0]
    # nacelle_spline = NURBS.NURBScurve(NURB(2, knotsN, nacelle_weights), nacelle_cps)

    ### --- vary z-pos of 3rd control point to maintain angle --- ###
    knots3 = TF[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    if duct_te_r <= Rtip
        lower_theta = -atan(delta_r_te, (duct_chord - duct_outlet_cpz))
        ncp3z = duct_chord - (nacelle_cpr - duct_te_r) / tan(wedge_angle + lower_theta)
        nacelle_cps = [
            SVector(0.0, duct_le_r, 0.0),
            SVector(0.0, duct_le_r + sqrt(duct_le_radius) * ncp3z, 0.0),
            SVector(ncp3z, nacelle_cpr, 0.0),
            SVector(duct_chord, duct_te_r, 0.0),
        ]

    else
        # ncp3z0 = duct_chord - (nacelle_cpr - Rtip) / tan(wedge_angle)
        # upper_theta = atan((nacelle_cpr - duct_te_r), duct_chord - ncp3z0)
        # println("upper_theta = ", upper_theta * 180 / pi)
        # nacelle_cpr += tan(upper_theta - wedge_angle) * (duct_chord - ncp3z0)
        # println("nacelle_cpr = ", nacelle_cpr)
        # nacelle_cps = [
        #     SVector(0.0, duct_le_r, 0.0),
        #     SVector(0.0, duct_le_r + sqrt(duct_le_radius) * ncp3z0, 0.0),
        #     SVector(ncp3z0, nacelle_cpr, 0.0),
        #     SVector(duct_chord, duct_te_r, 0.0),
        # ]
        # nacelle_spline = NURBS.NURBScurve(
        #     NURB(3, knots3, ones(TF, size(nacelle_cps, 1))), nacelle_cps
        # )

        ## --- Nominal, does not preserve trailing edge angle --- ###
        # knot vector
        # knots3 = TF[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
        # control points
        nacelle_cps = [
            SVector(0.0, duct_le_r, 0.0),
            SVector(0.0, nacelle_cpr, 0.0),
            SVector(nacelle_cpz, nacelle_cpr, 0.0),
            SVector(duct_chord, duct_te_r, 0.0),
        ]
    end
    # spline
    nacelle_spline = NURBS.NURBScurve(
        NURB(3, knots3, ones(TF, size(nacelle_cps, 1))), nacelle_cps
    )

    # - Get Casing Points - #
    u = collect(TF, range(0.0, 1.0; length=N + 1))
    Cw_inlet = inlet_spline(u)
    Cw_outlet = outlet_spline(u)

    # assemble array
    casing_zs = [getindex.(Cw_inlet, 1); flatx[2:(end - 1)]; getindex.(Cw_outlet, 1)]
    casing_rs = [getindex.(Cw_inlet, 2); flatr[2:(end - 1)]; getindex.(Cw_outlet, 2)]

    # - Get Nacelle Points - #
    u = collect(TF, range(0.0, 1.0; length=N + 1))
    Cw_nacelle = nacelle_spline(u)

    # assemble array
    nacelle_zs = getindex.(Cw_nacelle, 1)
    nacelle_rs = getindex.(Cw_nacelle, 2)

    if smooth
        # - Respline the whole geometry to make things smooth - #
        # re-spline with akima spline and cosine spacing
        cz = scaled_cosine_spacing(N, duct_chord, duct_le_z)
        csp = fmspline(casing_zs, casing_rs)
        cr = csp(cz)

        # re-spline with akima spline and cosine spacing
        nz = scaled_cosine_spacing(N, duct_chord, duct_le_z)
        nsp = fmspline(nacelle_zs, nacelle_rs)
        nr = nsp(nz)

        return cz, cr, nz, nr, inlet_cps, outlet_cps, nacelle_cps
    else
        return casing_zs,
        casing_rs, nacelle_zs, nacelle_rs, inlet_cps, outlet_cps,
        nacelle_cps
    end
end

"""

see centerbody_geom and duct_geom for input descriptions
"""
function bodies_from_params(
    Rhub,
    Rtip,
    duct_chord,
    duct_te_r,
    duct_outlet_start,
    duct_le_z,
    duct_le_r,
    # duct_le_radius,
    # duct_te_camber_angle,
    # duct_wedge_angle,
    # duct_alpha,
    duct_inlet_stop,
    duct_outlet_cpz,
    nacelle_cpz,
    nacelle_thickness_var,
    cb_nc_le,
    cb_ncp_z,
    cb_nc_stop,
    cb_tc_start,
    cb_tcp_z,
    cb_te_z,
    cb_te_r;
    N=60,
    smooth=false,
)

    # - Get Centerbody Coordinates - #
    cbz, cbr, _, _ = centerbody_geom(
        Rhub,
        duct_chord;
        cb_nc_le=cb_nc_le,
        cb_nc_stop=cb_nc_stop,
        cb_tc_start=cb_tc_start,
        cb_ncp_z=cb_ncp_z,
        cb_tcp_z=cb_tcp_z,
        cb_te_r=cb_te_r,
        cb_te_z=cb_te_z,
        N=N,
        smooth=smooth,
    )

    # # - Calculate duct_le_r and duct_te_r
    # if cb_nc_le <= 0.0
    #     r_hub_inlet = cbsp(cb_nc_le * duct_chord)
    # else
    #     r_hub_inlet = 0.0
    # end
    # r_hub_exit = cbsp(duct_chord)
    # duct_le_r = sqrt(duct_inlet_area / pi + r_hub_inlet^2)
    # duct_te_r = sqrt(duct_outlet_area / pi + r_hub_exit^2)

    # # - Get duct geometry coordinates - #
    cz, cr, nz, nr, _, _, _ = duct_geom(
        Rtip,
        duct_chord;
        duct_le_r=duct_le_r,
        duct_te_r=duct_te_r,
        duct_le_z=duct_le_z,
        duct_inlet_stop=duct_inlet_stop,
        duct_outlet_start=duct_outlet_start,
        duct_outlet_cpz=duct_outlet_cpz,
        nacelle_cpz=nacelle_cpz,
        nacelle_thickness_var=nacelle_thickness_var,
        N=N,
        smooth=smooth,
    )

    # cz, cr, nz, nr, _, _, _ = duct_geom(
    #     Rtip,
    #     duct_chord,
    #     duct_le_r,
    #     duct_le_radius;
    #     duct_te_r=duct_te_r,
    #     wedge_angle=wedge_angle,
    #     duct_le_z=duct_le_z,
    #     duct_inlet_stop=duct_inlet_stop,
    #     duct_outlet_start=duct_outlet_start,
    #     duct_outlet_cpz=duct_outlet_cpz,
    #     nacelle_cpz=nacelle_cpz,
    #     nacelle_thickness_var=nacelle_thickness_var,
    #     N=N,
    #     smooth=smooth,
    # )

    # nz, nr, cz, cr, _, _ = duct_geom(
    #     Rtip,
    #     duct_chord,
    #     duct_le_radius,
    #     duct_te_camber_angle,
    #     duct_wedge_angle;
    #     duct_alpha=duct_alpha,
    # )

    # assemble full duct geometry
    dz = [reverse(cz); nz[2:end]]
    dr = [reverse(cr); nr[2:end]]

    return [dz dr], [cbz cbr]
end
