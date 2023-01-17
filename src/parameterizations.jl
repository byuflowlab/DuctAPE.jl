#=

NURBS Parameterizations for Duct and Hub Geometries

=#

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
    xpos=1.0 / 6.0,
    weights=nothing,
    degree=3,
    N=60,
)
    boat_tail_angle = pi * wedge_angle / 360 #boattail angle is half of wedge angle
    te_camber_angle *= pi / 180

    Pl = [
        [le_x, inlet_radius],
        [le_x, inlet_radius - sqrt(2.0 * r_le) / 3.0],
        [xpos, (outlet_radius + (2.0 * tan(te_camber_angle - boat_tail_angle) / 3.0))],
        [te_x, outlet_radius],
    ]

    Pu = [
        [le_x, inlet_radius],
        [le_x, inlet_radius + sqrt(2.0 * r_le) / 3.0],
        [xpos, (outlet_radius + (2.0 * tan(te_camber_angle + boat_tail_angle) / 3.0))],
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
    Cw = [[0.0; 0.0] for i in 1:length(u)]

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
    degree=3,
    weights=nothing,
    N=80,
    debug=false,
)

    # calculate postitions of unknown x control point coordinates
    x1 = hub_le + rmax / tan(nose_angle)
    x2 = hub_te - rmax / tan(tail_angle)

    # Assemble Control Points
    if cyl_x == nothing
        ctrlpts = [[hub_le, 0.0], [x1, rmax], [x2, rmax], [hub_te, 0.0]]

        knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]
    else
        ctrlpts = [
            [hub_le, 0.0],
            [x1, rmax],
            [cyl_x[1], rmax],
            [cyl_x[1], rmax],
            [cyl_x[2], rmax],
            [cyl_x[2], rmax],
            [x2, rmax],
            [hub_te, 0.0],
        ]

        #Calculate knot ratios:
        knot_ratio1 = (cyl_x[1] - hub_le) / (hub_te - cyl_x[1])
        knot_ratio2 = (cyl_x[2] - hub_le) / (hub_te - cyl_x[2])

        # Knot vector set such that flat portion at rmax is at rotor root leading and trailing edge.
        knots = [
            0.0,
            0.0,
            0.0,
            0.0,
            knot_ratio1,
            knot_ratio1,
            knot_ratio2,
            knot_ratio2,
            1.0,
            1.0,
            1.0,
            1.0,
        ]
    end

    # Keep weights as ones unless otherwise specified
    if isnothing(weights)
        weights = ones(length(ctrlpts))
    end

    # Set up NURBS object
    hubsp = NURBS(degree, knots, weights, ctrlpts)

    h_x, h_r = get_coordinates(hubsp; N=N)

    # return control points and knots if debugging
    if debug

        #initialize vector
        knot_points = [[0.0; 0.0] for i in 1:length(knots)]

        # keep track of repeated knots
        count = 1

        #loop through knots
        for i in 1:length(knots)

            #get point
            knot_points[i] = Splines.curvepoint(hubsp, knots[i])

            #if knot is repeated shift it up for display purposes
            if i > 1
                if knots[i] == knots[i - 1]
                    knot_points[i][2] += count * 0.01
                    count += 1
                else
                    count = 1
                end
            end
        end

        # return x and r coordiates
        return h_x, h_r, ctrlpts, knot_points
    else
        return h_x, h_r
    end
end
