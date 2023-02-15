#from web plot digitizer on figure 6.16 of marine propellers book (also figure 2.20 from original work: wake adapted ducted propellers that the other book took it from)

# rotor_raw_les = [0.206987553534413  0.1943127962085307
# 0.15090540503520256  0.29541864139020535
# 0.10062186669588846  0.3949447077409163
# 0.053210238427892076  0.4960505529225908
# 0.014450867052023364  0.5987361769352291
# -0.018455104146690937  0.6951026856240126
# -0.04276360847053651  0.7977883096366508
# -0.05256645572509977  0.895734597156398
# -0.050808609338045096  0.9936808846761453
# ]

# rotor_raw_tes = [0.8717274379274762  0.1943127962085307
# 0.9081308385612139  0.29541864139020535
# 0.9416623291236335  0.3949447077409163
# 0.9722853829365621  0.4960505529225908
# 0.9999999999999998  0.5987361769352291
# 1.0220073235989735  0.6951026856240126
# 1.0410514204311976  0.7977883096366508
# 1.0514615237103795  0.8973143759873617
# 1.0561095435078396  0.995260663507109
# ]
using Xfoil
using Splines
include("../../plots_default.jl")

rotor_diameter = 240.00 #mm
num_blades = 4 # assume that 4-55 is 4 blades at 55 degrees twist?
radial_positions = [0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9]

# r/R = 0.6 chord length based on known 240mm diameter
c6 = 65.0 #mm

chords = [67.15; 76.59; 85.19; 93.01; 100.0; 105.86; 110.08; 112.66; 112.88] ./ 100 .* c6
twists = 55.0 .* ones(length(chords)) #twist is constant, but what is it?

# Airfoil Geoemetries
# max thickness values are percentages of rotor diameter
max_thickness =
    [4.0; 3.52; 3.0; 2.45; 1.9; 1.38; 0.92; 0.61; 0.5] ./ 100.0 .* rotor_diameter

# location of max thickness from leading edge is in percentage of the chord length
max_thickness_x =
    [34.98; 39.76; 46.02; 49.13; 49.98; 50.0; 50.0; 50.0; 50.0] ./ 100.0 .* chords

# x coordinates
xu_table = zeros(length(chords), 13)
xl_table = zeros(length(chords), 13)

for i in 1:length(chords)

    # in percentages from the leading/trailing edges to the location of maximum thickness
    LE_lengths = [0.0; 0.05; 0.1; 0.2; 0.4; 0.6; 0.8; 1.0] * max_thickness_x[i]
    TE_lengths =
        max_thickness_x[i] .+ [0.2; 0.4; 0.6; 0.8; 1.0] * (chords[i] .- max_thickness_x[i])
    xu_table[i, :] = [LE_lengths; TE_lengths] / chords[i] #divide by chord to normalize
    xl_table[i, :] = [LE_lengths; TE_lengths] / chords[i]
end

# r coordinates
# in percentages of maximum thickness.
# also divide by chord to normalize

rl_frac =
    [
        [33.33 20.62 16.04 10.52 04.37 01.46 00.21 00.00 00.00 00.10 01.77 07.29 20.21]
        [21.18 10.30 08.28 06.15 02.72 00.83 00.12 00.00 00.00 00.00 01.07 04.62 13.85]
        [13.47 04.44 03.89 02.92 01.39 00.42 00.00 00.00 00.00 00.00 00.56 02.36 09.17]
        [07.81 01.53 01.36 01.02 00.51 00.17 00.00 00.00 00.00 00.00 00.17 00.68 06.62]
        [00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00]
        [00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00]
        [00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00]
        [00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00]
        [00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00 00.00]
    ] ./ 100.0

ru_frac =
    [
        [00.00 27.40 38.75 55.00 77.19 90.83 97.92 100.0 95.00 82.40 63.65 38.23 00.00]
        [00.00 27.57 37.87 53.02 75.62 90.06 97.63 100.0 95.86 84.14 66.63 39.05 00.00]
        [00.00 25.83 34.72 50.00 73.61 88.89 97.22 100.0 96.25 85.69 66.94 40.56 00.00]
        [00.00 22.24 30.22 45.84 70.46 87.10 96.77 100.0 96.60 86.42 68.59 41.77 00.00]
        [00.00 20.44 28.59 43.58 68.26 85.89 96.47 100.0 96.47 85.89 68.26 43.58 00.00]
        [00.00 22.88 30.79 45.31 69.24 86.33 96.58 100.0 96.58 86.33 69.24 45.31 00.00]
        [00.00 26.90 34.39 48.16 70.84 87.04 96.76 100.0 96.76 87.04 70.84 48.16 00.00]
        [00.00 31.87 38.87 51.75 72.94 88.09 97.17 100.0 97.17 88.09 72.94 51.75 00.00]
        [00.00 32.31 39.25 52.00 73.00 88.00 97.00 100.0 97.00 88.00 73.00 52.00 00.00]
    ] ./ 100.0

rl_table = similar(rl_frac)
ru_table = similar(ru_frac)
for i in 1:length(chords)
    rl_table[i, :] = rl_frac[i, :] ./ chords[i] .* max_thickness[i]
    ru_table[i, :] = (ru_frac[i, :] .+ rl_frac[i, :]) ./ chords[i] .* max_thickness[i]
end

plot(; aspectratio=:equal)
plot!(
    xu_table[1, :],
    ru_table[1, :];
    seriestype=:scatter,
    markersize=2,
    label="Upper Ordiantes",
)
plot!(
    xl_table[1, :],
    rl_table[1, :];
    seriestype=:scatter,
    markersize=2,
    label="Lower Ordinates",
)
savefig("Ka-455_airfoil_test.pdf")

#---------------------------------#
#           SMOOTH DATA           #
#---------------------------------#

pts = [[0.0 0.0] for i in 1:25]
u = range(0.0, 1.0; length=160)
coordinates = zeros(160, 2, length(chords))
for i in 1:length(chords)
    xx = [reverse(xl_table[i, :]); xu_table[i, 2:end]]
    rr = [reverse(rl_table[i, :]); ru_table[i, 2:end]]

    for j in 1:25
        pts[j] = [xx[j] rr[j]]
    end

    nctrl = 19
    deg = 3

    bspline = Splines.leastsquarescurve(pts, nctrl, deg)
    for j in 1:length(u)
        coordinates[j, :, i] = Splines.curvepoint(bspline, u[j])
    end

    plot(; aspectratio=:equal)
    plot!(
        xu_table[i, :],
        ru_table[i, :];
        seriestype=:scatter,
        markersize=2,
        label="Upper Ordiantes",
    )
    plot!(
        xl_table[i, :],
        rl_table[i, :];
        seriestype=:scatter,
        markersize=2,
        label="Lower Ordinates",
    )
    plot!(
        coordinates[:, 1, i],
        coordinates[:, 2, i];
        seriestype=:scatter,
        markersize=1,
        label="BSpline Fit",
    )
    savefig("Ka-455_airfoil_test_bspline_smooth$(i*10+10).pdf")
    println("$(i*10+10) Airfoil:")
    println("Max Thickness = ", maximum(coordinates[:, 2, i]))
    println("Max Camber = ", maximum(coordinates[:, 2, i]) / 2.0)

    # - xfoil test - #
    # Xfoil.set_coordinates(reverse(coordinates[:, 1, i]), reverse(coordinates[:, 2, i]))
    # g = Xfoil.get_globals()
    # plot(; aspectratio=:equal)
    # plot!(
    #     xu_table[1, :],
    #     ru_table[1, :];
    #     seriestype=:scatter,
    #     markersize=2,
    #     label="Upper Ordiantes",
    # )
    # plot!(
    #     xl_table[1, :],
    #     rl_table[1, :];
    #     seriestype=:scatter,
    #     markersize=2,
    #     label="Lower Ordinates",
    # )
    # plot!(g.x[1:140], g.y[1:140]; seriestype=:scatter, markersize=1, label="PANE")
    # savefig("Ka-455_airfoil_test_bspline_smooth_xfoilpane.pdf")
    # cl, cd, cdp, cm, conv = Xfoil.solve_alpha(5.0, 1e6)
    # println("$(i*10+10) r/R converged: ", conv)
    # println("visc cl: ", cl)
    # cl, cm = Xfoil.solve_alpha(5.0)
    # println("inv cl: ", cl)
end

# using FLOWMath
# using FLOWFoil

# N = 80
# xu = FLOWFoil.cosine_spacing(N)
# xl = reverse(xu)

# ru = zeros(length(chords), N)
# rl = zeros(length(chords), N)

# for i in 1:length(chords)
#     ru[i, :] = FLOWMath.akima(xu_table[i, :], ru_table[i, :], xu)
#     rl[i, :] = reverse(FLOWMath.akima(xl_table[i, :], rl_table[i, :], xu))
# end

# plot!(
#     xu,
#     ru[1, :];
#     seriestype=:scatter,
#     color=:black,
#     markersize=1,
#     label="Smoothed Coordinates",
# )
# plot!(xl, rl[1, :]; seriestype=:scatter, color=:black, markersize=1, label="")
# savefig("Ka-455_airfoil_test_smooth.pdf")
# for i in 1:length(chords)
#     plot(; aspectratio=:equal)
#     plot!(xu, ru[i, :]; seriestype=:scatter, markersize=1, label="Smoothed Coordinates")
#     plot!(xl, rl[i, :]; seriestype=:scatter, markersize=1, label="")
#     savefig("Ka-455_airfoil_test_smooth_rR$(i*10+10).pdf")
# end

