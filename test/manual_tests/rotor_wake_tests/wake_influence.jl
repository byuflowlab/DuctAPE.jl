#=

Look at the wake influence at arbitrary locations within the wake, outside the wake, and on the rotor plane.

=#
include("../../../plots_default.jl")
using DuctTAPE
const dt = DuctTAPE
using FLOWFoil
const ff = FLOWFoil

#---------------------------------#
# Constant Rotor Circulation Case #
#---------------------------------#
# include constant_circulation_rotor file to have the values for the constant vorticity in the wake.
include("./constant_circulation_rotor.jl")

## -- Create a wake geometry -- ##
# use data from constant_circulation_rotor
R = maximum(radial_stations)
D = 2.0 * R
# Just make it a rectangle that is 2D long.
xrange = range(0.0, 2.0 * D; step=0.1)
grid_points = [[x r] for x in xrange, r in radial_stations]

# - Plot grid to max sure you made it correctly - #
plot(
    getindex.(grid_points, 1),
    getindex.(grid_points, 2);
    color=mycolors[1],
    aspectratio=:equal,
    seriestype=:scatter,
    markersize=1,
    label="",
)
savefig("./test/manual_tests/rotor_wake_tests/basic_wake_grid.pdf")

## -- Create Wake Panel Objects -- ##
#use the fuctions you already have in ducttape
#requires x points and r points separately, since that's how they are generated normally
#also, x is the first dimension and r the second
x_grid_points = repeat(xrange; outer=(1, length(radial_stations)))
r_grid_points = repeat(radial_stations'; outer=(length(xrange), 1))
nr = length(radial_stations)
#returns the wake panels as vectors for each "row"
wake_panels = dt.generate_wake_panels(
    x_grid_points,
    r_grid_points,
    nr;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [true]),
)

# - Plot the panels to make sure you made them correctly - #
# plot(; aspectratio=:equal)
for i in 1:nr
    plot!(
        wake_panels[i].panel_center[:, 1],
        wake_panels[i].panel_center[:, 2];
        color=mycolors[2],
        seriestype=scatter,
        markersize=1,
        label="",
    )
end
savefig("./test/manual_tests/rotor_wake_tests/basic_wake_panels.pdf")

## -- Determine points at which to sample velocity field -- ##
#how about the radial averages of the grid points?
x_sample_points = repeat(xrange; outer=(1, length(radial_stations) - 1)) .- 0.05
avg_radial_stations = (radial_stations[2:end] .+ radial_stations[1:(end - 1)]) ./ 2.0
r_sample_points = repeat(avg_radial_stations'; outer=(length(xrange), 1))

nra = length(avg_radial_stations)
#returns the wake panels as vectors for each "row"
sample_panels = dt.generate_wake_panels(
    x_sample_points,
    r_sample_points,
    nra;
    method=ff.AxisymmetricProblem(Vortex(Constant()), Neumann(), [true]),
)

# - Plot to make sure they're in the right place - #
for i in 1:(nr - 1)
    plot!(
        sample_panels[i].panel_center[:, 1],
        sample_panels[i].panel_center[:, 2];
        color=mycolors[3],
        seriestype=scatter,
        markersize=1,
        label="",
    )
end
savefig("./test/manual_tests/rotor_wake_tests/basic_wake_sample_points.pdf")

## -- Generate the geometry meshes for finding the velocity at all the points based on the wake panels -- ##
wake_meshes = [
    dt.generate_one_way_mesh(wake_panels[i], sample_panels[j]) for
    i in 1:length(wake_panels), j in 1:length(sample_panels)
]

## -- Generate the induced unit induced velocities -- ##
A_wake_to_wake = [
    dt.assemble_induced_velocity_matrices(
        wake_meshes[i, j], wake_panels[i], sample_panels[j]
    ) for i in 1:length(wake_panels), j in 1:length(sample_panels)
]
vxd_wake_to_wake = [
    A_wake_to_wake[i, j][1] for i in 1:length(wake_panels), j in 1:length(sample_panels)
]
vrd_wake_to_wake = [
    A_wake_to_wake[i, j][2] for i in 1:length(wake_panels), j in 1:length(sample_panels)
]

## -- Calculate the induced velocities inside the wake and on the rotor plane -- ##

# - need to populate the panel strengths based on the constant_circulation_rotor outputs. - #
panel_strengths = repeat(gamma_thetas; inner=(1, length(xrange) - 1))

# - loop through all the induced velocities summing up the totals - #
#initialize the outputs
vxs = zeros(nr - 1, length(xrange) - 1)
vrs = zeros(nr - 1, length(xrange) - 1)

# loop through the rows of sample panels
for i in 1:(nr - 1)
    # loop through the rows of vortex panels
    for j in 1:nr
        # populate the ith sample row with the contributions from the jth vortex panel row
        vxs[i, :] .+= vxd_wake_to_wake[j, i] * panel_strengths[j, :]
        vrs[i, :] .+= vrd_wake_to_wake[j, i] * panel_strengths[j, :]
    end
end

# - Put together some contour plots - #

pyplot()
contour(
    sample_panels[1].panel_center[:, 1],
    avg_radial_stations,
    vxs;
    aspectratio=:equal,
    fill=true,
)
savefig("./test/manual_tests/rotor_wake_tests/vx_contour.pdf")

contour(
    sample_panels[1].panel_center[:, 1],
    avg_radial_stations,
    vrs;
    aspectratio=:equal,
    fill=true,
)
savefig("./test/manual_tests/rotor_wake_tests/vr_contour.pdf")
