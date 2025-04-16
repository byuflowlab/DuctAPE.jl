######################################################################
#                                                                    #
#               NACA 65-series Compressor Airfoil                    #
#                                                                    #
######################################################################
"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 1
"""
function thickness65c(x; method="scaled")
    if method == "derived"
        tr =
            [
                0.0
                0.772
                0.932
                1.169
                1.574
                2.177
                2.647
                3.040
                3.666
                4.143
                4.503
                4.760
                4.924
                4.996
                4.963
                4.812
                4.530
                4.146
                3.682
                3.156
                2.584
                1.987
                1.385
                0.810
                0.306
                0.0
            ] * 1e-2
        ler = 0.687
    elseif method == "scaled"
        tr =
            [
                0.0
                0.752
                0.890
                1.124
                1.571
                2.222
                2.709
                3.111
                3.746
                4.218
                4.570
                4.824
                4.982
                5.057
                5.029
                4.870
                4.570
                4.151
                3.627
                3.038
                2.451
                1.847
                1.251
                0.749
                0.354
                0.150
            ] * 1e-2
        ler = 0.666
    else
        @error "no method $method, please choose scaled or derived."
    end
    tx =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2

    return FLOWMath.akima(tx, tr, x)
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function camber65c(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    a1y =
        [
            0.0
            0.25
            0.35
            0.535
            0.93
            1.580
            2.120
            2.585
            3.365
            3.98
            4.475
            4.86
            5.15
            5.355
            5.475
            5.515
            5.475
            5.355
            5.15
            4.86
            4.475
            3.98
            3.365
            2.585
            1.58
            0.0
        ] * 1e-2

    a1fine = FLOWMath.akima(a1x, a1y, x)

    return a1fine * clo
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function slopes65c(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    dydx = [
        0.4212
        0.38875
        0.3477
        0.29155
        0.2343
        0.19995
        0.17485
        0.13805
        0.11030
        0.08745
        0.06745
        0.04925
        0.03225
        0.01595
        0.0
        -0.01595
        -0.03225
        -0.04925
        -0.06745
        -0.08745
        -0.11030
        -0.13805
        -0.17485
        -0.23430
    ]

    dydxfine = FLOWMath.akima(a1x[2:(end - 1)], dydx, x)

    return dydxfine * clo
end

"""
Assumes x in non-dimensional range [0.0,1.0]


Description from NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds:"

The 65-series compressor blade family is formed by combining a basic thickness form with cambered mean lines.
The basic thickness form used is the NACA 65(216)-010 thickness form with the ordinates increased by 0.0015 times the chordwise stations to provide slightly, increased thickness toward the trailing edge.
In the scaled case, it was not derived for 10-percent thickness but was scaled down from the NACA 65,2-016 airfoil.
The scaling procedure gives the best results whep it is restricted to maximum thickness changes of a few percent.
The NACA 65-010 basic thickness has also been derived.
These thickness forms differ slightly but are considered to be interchangeable.

The basic mean line used is the a=1.0 mean line.
The amount of camber is for the design lift coefficient for the isolated airfoil with cl_o of 1.0.
Both ordinates and slopes are scaled directly to obtain other cambers.
Cambered blade sections are obtained by applying the thickness perpendicular to the mean line at stations laid out along the chord line.
In the designation the camber is given by the first number after the dash in tenths of cl_o.
For example, the NACA 65-810 and NACA 65-(12)10 blade sections are cambered for cl_o = 0.8 and cl_o = 1.2, respectively.
"""
function naca65c(clo; smoothed_leading_edge=true, method="scaled", N=161, x=nothing, split=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = cosine_spacing(N)
    end

    leading_edge_radius = 0.666*1e-2
    theta_initial = slopes65c(clo, 0.005)
    cos_theta_initial = cos(theta_initial)
    leading_edge_circle_center_x = leading_edge_radius * cos_theta_initial
    x = linear_transform(
        (0, 1),
        (-leading_edge_radius * cos_theta_initial + leading_edge_circle_center_x, 1),
        x,
    )

    t = thickness65c(x; method=method)
    c = camber65c(clo, x)
    s = slopes65c(clo, x)

    if !smoothed_leading_edge

        #y-positions at chordwise stations
        yl = c .- t
        yu = c .+ t

        if split
            return reverse(x), x, reverse(yl), yu
        else
            return [reverse(x); x[2:end]], [reverse(yl); yu[2:end]]
        end

    else

        #define output vectors
        y_upper = similar(x) .= 0.0
        y_lower = similar(x) .= 0.0
        x_upper = similar(x) .= 0.0
        x_lower = similar(x) .= 0.0

        cos_theta = cos.(s)
        sin_theta = sin.(s)
        tan_theta = tan.(s)

        for i in 1:N
            if x[i] < 0.005
                c[i] = camber65c(clo, 0.005)
                tan_theta[i] = tan(theta_initial)
                sin_theta[i] = sin(theta_initial)
                cos_theta[i] = cos(theta_initial)
                rad_a =
                    leading_edge_radius^2 -
                    ((leading_edge_circle_center_x - x[i]) / cos_theta[2])^2
                if rad_a < 0
                    a = 0.0
                else
                    a = sqrt(rad_a)
                end
                x_upper[i] = x[i] - a * sin_theta[i]
                x_lower[i] = x[i] + a * sin_theta[i]
                y_upper_1 = x[i] * tan_theta[i] + a * cos_theta[i]
                y_upper_2 = c[i] + t[i] * cos_theta[i]
                y_lower_1 = x[i] * tan_theta[i] - a * cos_theta[i]
                y_lower_2 = c[i] - t[i] * cos_theta[i]
                if y_upper_1 >= y_upper_2
                    y_upper[i] = y_upper_1
                else
                    y_upper[i] = y_upper_2
                end
                if y_lower_1 <= y_lower_2
                    y_lower[i] = y_lower_1
                else
                    y_lower[i] = y_lower_2
                end
                transition_index = i
            else
                x_upper[i] = x[i] - t[i] * sin_theta[i]
                y_upper[i] = c[i] + t[i] * cos_theta[i]
                x_lower[i] = x[i] + t[i] * sin_theta[i]
                y_lower[i] = c[i] - t[i] * cos_theta[i]
            end
        end

        if split
            return reverse(x_lower), x_upper, reverse(y_lower), y_upper
        else
            return [reverse(x_lower); x_upper[2:end]], [reverse(y_lower); y_upper[2:end]]
        end
    end
end
