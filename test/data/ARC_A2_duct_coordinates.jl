x = [
    0.0
    0.06
    0.3
    0.6
    1.2
    1.8
    2.4
    3.0
    3.6
    4.2
    4.8
    5.4
    6.0
    6.6
    7.2
    7.8
    8.4
    9.0
    9.6
    10.2
    10.8
    11.4
    12.0
]

ri = [
    5.69
    5.66
    5.62
    5.62
    5.64
    5.66
    5.69
    5.71
    5.72
    5.75
    5.77
    5.79
    5.8
    5.82
    5.83
    5.83
    5.83
    5.83
    5.82
    5.8
    5.77
    5.69
    5.6
]

ro = [
    5.69
    5.79
    5.93
    6.04
    6.19
    6.3
    6.37
    6.42
    6.45
    6.47
    6.47
    6.46
    6.44
    6.4
    6.37
    6.33
    6.27
    6.21
    6.13
    6.05
    5.95
    5.85
    5.66
]

x = x ./ x[end]
ri = ri ./ x[end]
ro = ro ./ x[end]

# 4 inch inner body 0 for inner body x is 16 inches in front of duct leading edge.
#nose is ellipse 2:1 length to radius.
xnrange = [0; 4.0] ./ x[end]
#most of the body is a cylinder
xcylrange = [4; 60.0] ./ x[end]
#the tail is just a wedge to the midline
xtrange = [60; 72.0] ./ x[end]

