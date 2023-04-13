#=
APC 5 x 5.75 Propeller data

From APC geometry database

units: inches, degrees
=#
using FLOWMath
include("../../plots_default.jl")

i2m = 0.0254 #inch to meter conversion

Rtip = 2.5 * i2m
Rhub = 0.66 * i2m
B = 2

radial_positions =
    [
        0.6630
        0.6932
        0.7234
        0.7535
        0.7837
        0.8138
        0.8490
        0.9029
        0.9613
        1.0200
        1.0786
        1.1373
        1.1959
        1.2545
        1.3132
        1.3719
        1.4305
        1.4892
        1.5478
        1.6064
        1.6651
        1.7237
        1.7824
        1.8411
        1.8997
        1.9583
        2.0170
        2.0756
        2.1343
        2.1930
        2.2516
        2.3103
        2.3688
        2.4209
        2.4672
        # 2.5000
    ] * i2m

chords =
    [
        0.5499
        0.5485
        0.5470
        0.5456
        0.5442
        0.5428
        0.5411
        0.5384
        0.5354
        0.5321
        0.5285
        0.5247
        0.5204
        0.5157
        0.5105
        0.5047
        0.4983
        0.4913
        0.4835
        0.4750
        0.4657
        0.4554
        0.4443
        0.4321
        0.4189
        0.4046
        0.3891
        0.3724
        0.3544
        0.3351
        0.3145
        0.2923
        0.2672
        0.2242
        0.1466
        # 0.0001
    ] * i2m

#sweep in inches
sweep =
    [
        0.2049
        0.2061
        0.2072
        0.2083
        0.2092
        0.2101
        0.2109
        0.2121
        0.2129
        0.2135
        0.2136
        0.2135
        0.2130
        0.2121
        0.2108
        0.2092
        0.2073
        0.2049
        0.2022
        0.1991
        0.1957
        0.1918
        0.1876
        0.1830
        0.1780
        0.1726
        0.1668
        0.1606
        0.1540
        0.1470
        0.1396
        0.1319
        0.1222
        0.0960
        0.0378
        # -0.0859
    ] * i2m

thickness_ratios = [
    0.1249
    0.1232
    0.1216
    0.1201
    0.1186
    0.1172
    0.1156
    0.1134
    0.1112
    0.1092
    0.1075
    0.1060
    0.1047
    0.1037
    0.1029
    0.1024
    0.1021
    0.1020
    0.1019
    0.1017
    0.1016
    0.1015
    0.1014
    0.1013
    0.1012
    0.1011
    0.1010
    0.1008
    0.1007
    0.1006
    0.1005
    0.1004
    0.1003
    0.1002
    0.1001
    # 0.1000
]

twist = [
    56.6834
    56.9056
    56.8808
    56.6379
    56.1942
    55.5586
    54.5788
    52.8956
    51.1540
    49.4871
    47.8988
    46.3864
    44.9463
    43.5753
    42.2700
    41.0269
    39.8429
    38.7147
    37.6394
    36.6139
    35.6356
    34.7018
    33.8100
    32.9577
    32.1429
    31.3633
    30.6171
    29.9023
    29.2172
    28.5602
    27.9298
    27.3245
    26.7439
    26.2465
    25.8185
    # 25.6703
]

max_thickness =
    [
        0.0687
        0.0676
        0.0665
        0.0655
        0.0645
        0.0636
        0.0626
        0.0611
        0.0595
        0.0581
        0.0568
        0.0556
        0.0545
        0.0535
        0.0525
        0.0517
        0.0509
        0.0501
        0.0493
        0.0483
        0.0473
        0.0462
        0.0451
        0.0438
        0.0424
        0.0409
        0.0393
        0.0376
        0.0357
        0.0337
        0.0316
        0.0294
        0.0268
        0.0225
        0.0147
        # 0.0000
    ] * i2m

airfoil_ranges = [radial_positions[1]; 0.68 * i2m; 1.5 * i2m; radial_positions[end]]

airfoil_files = ["e63_coordinates.dat", "naca4412_coordinates.dat"]
airfoil_maxtc = [4.5; 4.5; 12.0; 12.0] / 100.0

save_path = "APC-5x5.75_airfoils/"

function write_apc_airfoil_geometries(
    airfoil_files,
    airfoil_ranges,
    radial_positions,
    max_thickness,
    chords,
    rotor_name="APC_";
    save_path="",
)

    # - Read in Raw Airfoil Coordinates - #
    f = open(airfoil_files[1])
    nx = countlines(f)
    close(f)
    x = zeros(nx)
    z = zeros(nx, length(airfoil_ranges))
    for iaf in 1:length(airfoil_files)
        open(airfoil_files[iaf]) do f
            ix = 1
            for line in eachline(f)
                parts = split(line)
                x[ix] = parse(Float64, parts[1])
                z[ix, iaf + 1] = parse(Float64, parts[2])
                ix += 1
            end
        end
    end

    # repeat first and last airfoils
    z[:, 1] .= z[:, 2]
    z[:, end] .= z[:, end - 1]

    # - Spline Airfoil z-coordinates - #
    # Assumes you use the same unit x-coordinates for all airofils
    coord_sp = [FLOWMath.Akima(airfoil_ranges, z[ix, :]) for ix in 1:nx]

    # - Spline Raw Airfoil Max Thicknesses - #
    maxt_sp = FLOWMath.Akima(radial_positions, max_thickness)
    # - Spline Airfoil Chords - #
    c_sp = FLOWMath.Akima(radial_positions, chords)

    saved_files = String[]
    for ir in 1:length(radial_positions)
        # - Get interpolated z-coordinates - #
        zi = [coord_sp[ix](radial_positions[ir]) for ix in 1:nx]

        # - Get scaling factor - #
        # split airfoil into top and bottom
        _, leidx = findmin(x)
        zbot = reverse(zi[1:leidx])
        ztop = zi[leidx:end]
        ts = (ztop .- zbot) * c_sp(radial_positions[ir])
        maxt, _ = findmax(ts)
        # find index of maximum thickness
        # scale to match thickness ratio
        thickness_scaling = maxt_sp(radial_positions[ir]) / maxt

        # - Apply Scaling Factor - #
        zi .*= thickness_scaling

        # - Write Geometry - #
        save_name = rotor_name * "r-$(radial_positions[ir]).coord"
        push!(saved_files, save_name)
        f = open(save_path * save_name, "w")
        for ix in 1:nx
            write(f, "$(x[ix]) $(zi[ix])\n")
        end
    end

    return saved_files
end

saved_files = write_apc_airfoil_geometries(
    airfoil_files,
    airfoil_ranges,
    radial_positions,
    max_thickness,
    chords,
    "APC-5x5.75_airfoil_";
    save_path=save_path,
)

function sanity_plot(saved_files)

    # generate colors
    gr = range(0.0, 255.0 - 93.0, length(saved_files)) ./ 255.0
    gg = range(46.0, 255.0 - 46.0, length(saved_files)) ./ 255.0
    gb = range(93.0, 255.0, length(saved_files)) ./ 255.0

    p = plot(; aspectratio=1, xlabel="x/c", ylabel="z/c")

    # - Read in Raw Airfoil Coordinates - #
    f = open(data_path * saved_files[1])
    nx = countlines(f)
    close(f)
    for iaf in 1:length(saved_files)
        x = zeros(nx)
        z = zeros(nx)
        open(data_path * saved_files[iaf]) do f
            ix = 1
            for line in eachline(f)
                parts = split(line)
                x[ix] = parse(Float64, parts[1])
                z[ix] = parse(Float64, parts[2])
                ix += 1
            end
        end
        # Plot
        plot!(p, x, z; label="station $iaf", color=RGB(gr[iaf], gg[iaf], gb[iaf]))
    end

    savefig(p, "APC-5x5.75_test_plot.pdf")

    return nothing
end

sanity_plot(saved_files)
