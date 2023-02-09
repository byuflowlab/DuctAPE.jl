######################################################################
#                                                                    #
#              GEOMETRY FROM ORIGINAL MARIN PUBLICATIONS             #
#                                                                    #
######################################################################

#=
The MARIN 19A duct has a length to diameter ratio of 0.5.
The r-coordinates are given with the inner most coordinates at zero, so we add the inverse of that ratio divided by two (for the radius) to all r-coordinates as a radial offset.
Scaling the length of the duct should then scale everything else as required to maintain the proper length to diameter ratio.
=#

L_D = 0.5
R_L = 1.0 / (2.0 * L_D)

# - Non-dimensional x-coordinates - #
x = [
    # 0.0
    0.0125
    0.025
    0.05
    0.075
    0.1
    0.15
    0.2
    0.25
    0.3
    0.4
    0.5
    0.6
    0.7
    0.8
    0.9
    0.95
    1.0
]

# - Non-dimensional inner surface r-coordinates - #
# Note that in the original, there is a circular cylidric portion of the inner r-coordinates which does not have any numeric values set, we set these values to zero (before adding the radial offset).
r_inner =
    [
        # 0.1825
        0.1466
        0.1280
        0.10#87
        0.08
        0.0634
        0.0387
        0.0217
        0.011
        0.0048
        0.0
        0.0
        0.0
        0.0029
        0.0082
        0.0145
        0.0186
        0.0236
    ] .+ R_L

#=
Note: the original geometry just says "straight line" for most of the coordinates,
so we need to provide a function to calculate the missing information.
Since the original x-coordinates have irregular spacing in this region, we'll add some interpolated points in this linear region to make things work nicely.
=#

# - Get circle coordinates at nose - #
le_radius = 5.57 / 200.0
center = [le_radius; R_L + 0.2107 - le_radius]
theta = range(255.0, 360.0; length=91)
x_circle = le_radius * sind.(theta) .+ center[1]
r_circle = le_radius * cosd.(theta) .+ center[2]

# - Non-dimensional outer surface x-coordinates - #
straight_line_x = range(0.05, 1.0; step=0.025)
x_outer = [0.0; 0.0125; 0.025; straight_line_x]

# - Non-dimensional outer surface r-coordinates - #
straight_line_r = range(0.208, 0.0636, length(straight_line_x))
r_outer = [0.1825; 0.2072; 0.2107; straight_line_r] .+ R_L

full_x = [reverse(x); x_circle; straight_line_x]
full_r = [reverse(r_inner); r_circle; straight_line_r .+ R_L]

