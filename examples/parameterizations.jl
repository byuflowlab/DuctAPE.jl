#=

NURBS Parameterizations for Duct and Hub Geometries

=#

using Splines
using FLOWMath
const fm = FLOWMath
using FLOWFoil
const ff = FLOWFoil
using DuctTAPE
const dt = DuctTAPE

"""
"""
function generate_duct_coordinates(;
    r_le=0.01, #0.02
    le_x=0.0,
    te_x=0.75,
    inlet_radius=0.135,
    outlet_radius=0.08,
    wedge_angle=2.0,
    te_camber_angle=15.0,
    ctrlpt_xpos=1.0 / 6.0,
    weights=nothing,
    degree=3,
    N=80,
)
    boat_tail_angle = pi * wedge_angle / 360 #boattail angle is half of wedge angle
    te_camber_angle *= pi / 180

    Pl = [
        [le_x, inlet_radius],
        [le_x, inlet_radius - sqrt(2.0 * r_le) / 3.0],
        [
            ctrlpt_xpos,
            (outlet_radius + (2.0 * tan(te_camber_angle - boat_tail_angle) / 3.0)),
        ],
        [te_x, outlet_radius],
    ]

    Pu = [
        [le_x, inlet_radius],
        [le_x, inlet_radius + sqrt(2.0 * r_le) / 3.0],
        [
            ctrlpt_xpos,
            (outlet_radius + (2.0 * tan(te_camber_angle + boat_tail_angle) / 3.0)),
        ],
        [te_x, outlet_radius],
    ]

    #--Initialize knot vector
    knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    if weights == nothing
        wu = [1.0 for i in 1:4]
        wl = [1.0 for i in 1:4]
    else
        wu = getindex(weights, 1)
        wl = getindex(weights, 2)
    end

    upper_spline = Splines.NURBS(degree, knots, wu, Pu)
    lower_spline = Splines.NURBS(degree, knots, wl, Pl)

    xu, ru = get_coordinates(upper_spline; N=N)
    xl, rl = get_coordinates(lower_spline; N=N)

    return xu, ru, xl, rl
end

"""
    get_coordinates()
Get plotable coordinates from the weighted control points.

**Arguments:**
- `knots::Array{Float}` : knot vector
- `controlpoints::Array{Float,2} `: Control Point matrix
- `degree::Int` : degree of the NURBS curve

**Keyword Arguments:**
- `N::Int` : number of coordinates, or panels that will be generated

**Returns:**
- `x::Array{Float}` : x Airfoil coordinates
- `z::Array{Float}` : z Airfoil coordinates
"""
function get_coordinates(nurbs; N=160)

    #create parametric point array
    u = range(0.0, 1.0; length=N + 1)

    #n = number of control points - 1
    n = length(nurbs.ctrlpts) - 1

    #initialize curve point array
    TF = eltype(nurbs.ctrlpts[1])
    Cw = [zeros(TF, 2) for i in 1:length(u)]

    #loop through parametric points to get curve points
    for i in 1:length(u)
        Cw[i] = Splines.curvepoint(nurbs, u[i])
    end

    return getindex.(Cw, 1), getindex.(Cw, 2)
end

"""
    generate_hub_coordinates(
    hub_le=0.0,
    hub_te=1.0,
    rmax=0.2,
    nose_angle=pi / 4.0,
    tail_angle=pi / 6.0,
    cyl_x=nothing,
    degree=3,
    weights=nothing,
    N=60,
    debug=false,
)

"""
function generate_hub_coordinates(;
    hub_le=0.05,
    hub_te=2.0 / 3.0,
    rmax=0.055,
    nose_angle=pi / 4.0,
    tail_angle=pi / 9.0,
    cyl_x=[0.2125, 0.26],
    N=80,
    debug=false,
)

    # calculate postitions of unknown x control point coordinates
    xnose = hub_le + rmax / tan(nose_angle)
    xtail = hub_te - rmax / tan(tail_angle)

    TF = typeof(xnose)

    # Assemble Control Points
    if cyl_x == nothing
        degree = 3
        ctrlpts = [[hub_le, 0.0], [xnose, rmax], [xtail, rmax], [hub_te, 0.0]]

        knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    else
        degree = 2
        nose_pts = [[hub_le, 0.0], [xnose, rmax], [cyl_x[1], rmax]]

        tail_pts = [[cyl_x[2], rmax], [xtail, rmax], [hub_te, 0.0]]

        # Knot vector set such that flat portion at rmax is at rotor root leading and trailing edge.
        knots = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    end

    # Keep weights as ones unless otherwise specified
    nose_weights = ones(length(nose_pts))
    tail_weights = ones(length(tail_pts))

    # Set up NURBS object
    nose_sp = NURBS(degree, knots, nose_weights, nose_pts)
    tail_sp = NURBS(degree, knots, tail_weights, tail_pts)

    # Extract coordinates
    n_x, n_r = get_coordinates(nose_sp; N=N)
    t_x, t_r = get_coordinates(tail_sp; N=N)

    # - Add in cylinder section - #
    c_x = range(cyl_x[1], cyl_x[2]; step=n_x[end] - n_x[end - 1])
    c_r = rmax .* ones(TF, length(c_x))

    # - Put everything together - #
    h_x = [n_x; c_x[2:(end - 1)]; t_x]
    h_r = [n_r; c_r[2:(end - 1)]; t_r]

    # - Smooth things out - #

    x = dt.lintran([hub_le; hub_te], [0.0; 1.0], ff.cosine_spacing(N))
    r = fm.akima(h_x, h_r, x)

    return x, r
end
