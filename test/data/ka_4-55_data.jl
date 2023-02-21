#from web plot digitizer on figure 6.16 of marine propellers book (also figure 2.20 from original work: wake adapted ducted propellers that the other book took it from)

using Xfoil
using Splines
using CCBlade
include("../../plots_default.jl")

rotor_diameter = 240.00 #mm
num_blades = 4 # assume that 4-55 is 4 blades at 55 degrees twist?
radial_positions = [0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]

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

    # shift things down so the trailing edge it one the x-axis
    coordinates[:, 2, i] .-= coordinates[end, 2, i]

    f = open("test/data/ka_4-55_geometry_section_0.$(i+1)R.dat", "w")
    for j in 1:length(u)
        write(f, "$(coordinates[j,1,i]) $(coordinates[j,2,i])\n")
    end
    close(f)

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
end

### --- Figure out Reynolds Number to run thick sections at.

rho = 1.225 #air
# rho = 1000 #water
mu = 1.81e-5 #air
# mu = 1e-3 #water
c = chords
r = radial_positions * rotor_diameter / 2.0

J = 0.36
Vinf = 10.0
Omega = Vinf / (J * r[end] / pi)
for i in 1:length(chords)
    re = rho * c[i] * sqrt(Vinf^2 + (Omega * r[i])^2) / mu
    println("Re @ $(i*10+10)% R = ", re)
    inflow_angle = atan(Vinf / (Omega * r[i])) * 180 / pi
    aoa = 55.0 - inflow_angle
    println("AoA @ $(i*10+10)% R = ", aoa)
end

#
#
#
#
#
#
#

# - xfoil test - #
alpha = range(-15.0, 20.0; step=0.25)
for i in 1:4#length(chords)
    Xfoil.set_coordinates(reverse(coordinates[:, 1, i]), reverse(coordinates[:, 2, i]))

    f = open("test/data/ka4-55_0.$(i*10+10)R_polar.dat", "w")
    write(f, "")
    write(f, "Ka 4-55 0.$(i*10+10) percent radius section\n")
    write(f, "100000000\n")
    write(f, "0\n")
    println("foil $(i*10+10)")
    clvec = []
    cdvec = []
    idxs = []
    for j in 1:length(alpha)
        cl, cd, cdp, cm, conv = Xfoil.solve_alpha(alpha[j], 1e8; reinit=true)
        if conv
            println("Alpha: $(alpha[j]) converged")
            write(f, "$(alpha[j]*pi/180.0) $cl $cd\n")
            push!(clvec, cl)
            push!(cdvec, cd)
            push!(idxs, j)
        end
    end
    close(f)

    # apply corrections

    cr75 = 0.57
    rR = 0.75  # r/R = 75%
    tsr = 9.0  # representative tip-speed ratio
    alpha_ext, cl_ext, cd_ext = viterna(alpha[idxs] * pi / 180.0, clvec, cdvec, cr75)
    plot(alpha_ext, cl_ext)
    savefig("clextrap_test$i.pdf")

    f = open("test/data/ka4-55_0.$(i*10+10)R_polar_extrot.dat", "w")
    write(f, "")
    write(
        f,
        "Ka 4-55 0.$(i*10+10) percent radius section extrapolated and rotation coreected\n",
    )
    write(f, "100000000\n")
    write(f, "0\n")
    cl_rot = similar(cl_ext)
    cd_rot = similar(cd_ext)
    for j in 1:length(alpha_ext)
        cl_rot[j], cd_rot[j] = rotation_correction(
            DuSeligEggers(), cl_ext[j], cd_ext[j], cr75, rR, tsr, alpha_ext[j]
        )
        write(f, "$(alpha_ext[j]) $(cl_rot[j]) $(cd_rot[j])\n")
    end
    close(f)

    plot(alpha_ext, cl_rot)
    savefig("clrot_test$i.pdf")
end
