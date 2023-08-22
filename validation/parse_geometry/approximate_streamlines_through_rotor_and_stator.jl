include("geometry_parsing_functions.jl")

# conversion factor inch to meter
i2m = 0.0254

# - User Inputs - #
# number of interpolation sections
ninterp = 25

#number of streamline x coordinates per section (between blades)
N = 100

datapath = "NASA_geometry_files/"
savepath = "CSVs/"

## -- LOAD GEOMETRIES -- ##
# - Load body geometries - #
file = datapath * "casing.xr" #in inches
x, r = parse2d(file)
duct_inner = [x r] .* i2m

file = datapath * "centerbody.xr" #in inches
x, r = parse2d(file)
hub = [x r] .* i2m

# - Load rotor geometry - #
file = datapath * "rotor.xyz"
rotor_sections = parse3d(file) .* i2m

# - Load stator geometry - #
file = datapath * "stator.xyz"
stator_sections = parse3d(file) .* i2m
rotor_x, rotor_r, rotor_id = rotor_edges(stator_sections)

## -- Interpolate Rotor Sections -- ##
interp_rotor_sections = interpolate_sections(
    rotor_sections,
    rotor_sections[:, 1, 3];
    nsections=ninterp,
    dimrange=nothing,
    cutoff=0.03,
)

# - Get the LE coordinates for the rotor sections - #
rotor_lex, rotor_ler, rotor_leid = rotor_edges(interp_rotor_sections)

# - Get the x stations for the inlet - #
xinlet = zeros(length(rotor_lex), N)
for is in 1:length(rotor_lex)
    xinlet[is, :] = range(hub[1, 1], rotor_lex[is], N)
end

# - Get the TE coordinates for the rotor sections - #
rotor_tex, rotor_ter, rotor_teid = rotor_edges(interp_rotor_sections; edge="TE")

# - Get the x stations for the interior - #
xinterior = zeros(length(rotor_tex), N)
for is in 1:length(rotor_tex)
    xinterior[is, :] = range(rotor_tex[is], minimum(stator_sections[:, :, 1]), N)
end

## -- Calculate streamlines from LE to rotor -- ##
xinlet, rinlet = approximate_streamlines(
    duct_inner, hub, xinlet, rotor_lex, rotor_ler, rotor_leid
)

## -- Calculate streamlines from rotor to stator -- ##
xinterior, rinterior = approximate_streamlines(
    duct_inner, hub, xinterior, rotor_tex, rotor_ter, rotor_teid
)

## -- Interpolate stator Sections -- ##
_, statorleid = findmin(stator_sections[1, :, 1])
interp_stator_sections = interpolate_sections(
    stator_sections,
    sqrt.(stator_sections[:, statorleid, 3] .^ 2 .+ stator_sections[:, statorleid, 2] .^ 2);
    nsections=ninterp,
    dimrange=rinterior[:, end],
)

# - Get the TE coordinates for the stator sections - #
stator_lex, stator_ler, stator_leid = rotor_edges(interp_stator_sections; edge="LE")
stator_tex, stator_ter, stator_teid = rotor_edges(interp_stator_sections; edge="TE")

# - Get the x stations for the inlet - #
xoutlet = zeros(length(stator_tex), N)
for is in 1:length(rotor_lex)
    xoutlet[is, :] = range(stator_tex[is], duct_inner[end, 1], N)
end

## -- Calculate streamlines from rotor to stator -- ##
xoutlet, routlet = approximate_streamlines(
    duct_inner, hub, xoutlet, stator_tex, stator_ter, stator_teid
)

##### ----- ASSEMBLE STREAMLINES AND WRITE TO FILES ----- #####

for i in 1:ninterp #this is the number of streamlines we have
    streamline_x = [
        xinlet[i, :]
        interp_rotor_sections[i, rotor_leid[i]:end, 1]
        xinterior[i, :]
        interp_stator_sections[i, stator_leid[i]:end, 1]
        xoutlet[i, :]
    ]

    # set up for calclulating r
    ry = interp_rotor_sections[i, rotor_leid[i]:end, 2]
    rz = interp_rotor_sections[i, rotor_leid[i]:end, 3]
    sy = interp_stator_sections[i, stator_leid[i]:end, 2]
    sz = interp_stator_sections[i, stator_leid[i]:end, 3]

    streamline_r = reduce(
        vcat,
        [
            rinlet[i, :]
            [[sqrt(ry[c]^2 + rz[c]^2)] for c in 1:length(ry)]
            rinterior[i, :]
            [[sqrt(sy[c]^2 + sz[c]^2)] for c in 1:length(sy)]
            routlet[i, :]
        ],
    )

    # plot!(streamline_x, streamline_r; linewidth=2, color=:orange, label="")

    f = open(savepath * "streamline_section_$(i).csv", "w")
    write(f, "x,r\n")
    for j in 1:length(streamline_x)
        write(f, "$(streamline_x[j]),$(streamline_r[j])\n")
    end
    close(f)
end

#---------------------------------#
#             PLOTTING            #
#---------------------------------#
# plot(; aspectratio=1, legend=:outertopright)
plot(; xlabel="x (m)", ylabel="r (m)", aspectratio=1, legend=:outerright)

# - Plot Duct and Hub - #
plot!(duct_inner[:, 1], duct_inner[:, 2]; color=myblue[1], label="Casing")
plot!(hub[:, 1], hub[:, 2]; color=myred[1], label="Center Body")

# - Plot interpolated blade element streamlines (based on nasa geometry) - #
for i in 1:length(interp_rotor_sections[:, 1, 1])
    ry = interp_rotor_sections[i, rotor_leid[i]:end, 2]
    rz = interp_rotor_sections[i, rotor_leid[i]:end, 3]
    plot!(
        interp_rotor_sections[i, rotor_leid[i]:end, 1],
        reduce(vcat, [[sqrt(ry[c]^2 + rz[c]^2)] for c in 1:length(ry)]);
        color=myblue[3],
        linewidth=0.5,
        label="",
    )
end
for i in 1:length(interp_stator_sections[:, 1, 1])
    sy = interp_stator_sections[i, stator_leid[i]:end, 2]
    sz = interp_stator_sections[i, stator_leid[i]:end, 3]
    lab = i == 1 ? "Blade Element Streamlines" : ""
    plot!(
        interp_stator_sections[i, stator_leid[i]:end, 1],
        reduce(vcat, [[sqrt(sy[c]^2 + sz[c]^2)] for c in 1:length(sy)]);
        color=myblue[3],
        linewidth=0.5,
        label=lab,
    )
end

# - Plot approximated streamlines inside the rest of the duct, based on conservation of mass - #
for ir in 1:length(rinlet[:, 1])
    lab = ir == 1 ? "Approximated Streamlines" : ""
    plot!(xinlet[ir, :], rinlet[ir, :]; color=mygray[1], linewidth=0.25, label=lab)
end

for ir in 1:length(rinterior[:, 1])
    plot!(xinterior[ir, :], rinterior[ir, :]; color=mygray[1], linewidth=0.25, label="")
end

for ir in 1:length(routlet[:, 1])
    plot!(xoutlet[ir, :], routlet[ir, :]; color=mygray[1], linewidth=0.25, label="")
end

savefig(savepath * "../figures/approximated_streamlines.pdf")
# savefig(savepath * "../figures/approximated_streamlines.png")
