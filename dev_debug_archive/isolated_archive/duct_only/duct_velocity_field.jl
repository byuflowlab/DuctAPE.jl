
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

using DuctAPE
const dt = DuctAPE

using Statistics

#---------------------------------#
#          Setup Geometry         #
#---------------------------------#

npan = 40
# Duct Geometry
include(project_dir * "/test/data/naca_662-015.jl")
scale = 0.5

# x_duct = x_duct .* scale
# r_duct = r_duct .* scale
# r_duct .-= 0.3

x_duct_raw = x_duct .* scale
r_duct_raw = r_duct .* scale
r_duct_raw .-= 0.3
#split into upper/lower
ducttrans, leidx = findmin(x_duct_raw)
innerx = x_duct_raw[1:leidx]
outerx = x_duct_raw[(leidx + 1):end]
innerr = r_duct_raw[1:leidx]
outerr = r_duct_raw[(leidx + 1):end]
x_duct_fine = dt.scaled_cosine_spacing(npan, scale, ducttrans)
r_duct_in = FLOWMath.akima(reverse(innerx), reverse(innerr), x_duct_fine)
r_duct_out = FLOWMath.akima(outerx, outerr, x_duct_fine)
x_duct = [reverse(x_duct_fine); x_duct_fine[2:end]]
r_duct = [reverse(r_duct_in); r_duct_out[2:end]]

duct_coordinates = [x_duct r_duct]

# - Hub Coordinates - #
Rhub = minimum(r_duct) * 0.1
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
x_hub_raw = x_hub[1:(end - 2)] .* scale
r_hub_raw = r_hub[1:(end - 2)] * Rhub / maximum(r_hub)

x_hub = dt.scaled_cosine_spacing(npan, x_hub_raw[end] - x_hub_raw[1], x_hub_raw[1])
r_hub = FLOWMath.akima(x_hub_raw, r_hub_raw, x_hub)

hub_coordinates = [x_hub r_hub]

#---------------------------------#
#            Solve Duct           #
#---------------------------------#

##### ----- Panels ----- #####
body_panels = dt.generate_body_panels(duct_coordinates, hub_coordinates)

##### ----- Mesh ----- #####
mesh_bb = dt.generate_one_way_mesh(body_panels, body_panels)

kutta_idxs = dt.get_kutta_indices([false, true], mesh_bb)

##### ----- Linear System ----- #####
A_bb = dt.assemble_induced_velocity_on_body_matrix(
    mesh_bb, body_panels, body_panels; singularity="vortex"
)

# apply back-diagonal correction to duct portions of coefficient matrix
dt.apply_back_diagonal_correction!(
    A_bb, body_panels[1], mesh_bb.affect_panel_indices[1], mesh_bb.mesh2panel_affect
)

# - freestream to body - #
Vinf = 5.0
b_bf = Vinf .* dt.assemble_body_freestream_boundary_conditions(body_panels, mesh_bb)

body_strengths = dt.solve_body_system(A_bb, b_bf, kutta_idxs) # get circulation strengths from solving body to body problem

#---------------------------------#
#     check at "rotor" plane      #
#---------------------------------#
# simply generte a "rotor" set of panels at the point you want to do a mass conservation ballpark check...
rotorzloc = range(0.1, 0.9; length=5) .* scale
ridxd = [findfirst(x -> x < rotorzloc[i], x_duct) for i in 1:length(rotorzloc)]
ridxh = [findfirst(x -> x > rotorzloc[i], x_hub) for i in 1:length(rotorzloc)]
rs = [range(r_hub[ridxh[i]], r_duct[ridxd[i]]; length=100) for i in 1:length(rotorzloc)]

# rotor source panel objects
rotor_source_panels = [dt.generate_rotor_panels(rotorzloc[i], rs[i]) for i in 1:length(rotorzloc)]

#rotor panel centers
rpc = [rotor_source_panels[i].panel_center[:, 2] for i in 1:length(rotorzloc)]

# - body to rotor unit induced velocities- #
mesh_rb = [
    dt.generate_one_way_mesh(body_panels, rotor_source_panels[i]) for
    i in 1:length(rotor_source_panels), j in 1:1
]

A_rb = [
    dt.assemble_induced_velocity_matrices(
        mesh_rb[i, j], body_panels, rotor_source_panels[i]
    ) for i in 1:length(rotor_source_panels), j in 1:1
]

# axial components
vz_rb = [A_rb[i, j][1] for i in 1:length(rotor_source_panels), j in 1:1]

# radial components
vr_rb = [A_rb[i, j][2] for i in 1:length(rotor_source_panels), j in 1:1]

#---------------------------------#
#  ballpark conservation of mass  #
#---------------------------------#
#inlet area
_, leidx = findmin(x_duct)
Ai = pi * r_duct[leidx]^2
#area of interest
As = pi .* [(r_duct[ridxd[i]] - r_hub[ridxh[i]])^2 for i in 1:length(rotorzloc)]

Vs = Vinf .* Ai ./ As

#---------------------------------#
#             PLOTS               #
#---------------------------------#
pgeom = plot(; aspectratio=1, xlabel="x", ylabel="r")
plot!(pgeom, x_duct, r_duct; color=:black, linewidth=0.5, label="")
plot!(pgeom, x_hub, r_hub; color=:black, linewidth=0.5, label="")
plot!(
    pgeom,
    body_panels[1].panel_center[:, 1],
    body_panels[1].panel_center[:, 2];
    seriestype=:scatter,
    color=:black,
    markersize=1,
    label="",
)
plot!(
    pgeom,
    body_panels[2].panel_center[:, 1],
    body_panels[2].panel_center[:, 2];
    seriestype=:scatter,
    color=:black,
    markersize=1,
    label="",
)
plot!(
    pgeom,
    zeros(2),
    [0.0; r_duct[leidx]];
    linestyle=:dash,
    color=:black,
    label="reference area",
)

pv = plot(; xlabel="induced x velocity + Vinf", ylabel="r")

for i in 1:length(rotorzloc)
    plot!(
        pgeom,
        rotorzloc[i] * ones(length(rpc[i])),
        rpc[i];
        # linestyle=:dash,
        color=mycolors[i],
        label="rotor location",
    )

    body_induced_x = vz_rb[i] * body_strengths .+ Vinf

    plot!(pv, body_induced_x, rpc[i]; color=mycolors[i], label="x location = $(rotorzloc[i])")

    plot!(
        pv,
        Statistics.mean(body_induced_x) .* ones(length(rpc[i])),
        rpc[i];
        linestyle=:dash,
        color=mycolors[i],
        label="Mean $(rotorzloc[i])",
    )

    plot!(
        pv,
        Vs[i] .* ones(length(rpc[i])),
        rpc[i];
        linestyle=:dot,
        color=mycolors[i],
        label="Cons Mass $(rotorzloc[i])",
    )
end

pg = plot(; xlabel="x", ylabel="Surface Velocity")
npd = length(body_panels[1].panel_center[:, 1])
_, splitidx = findmin(body_panels[1].panel_center[:, 1])
pdi = body_panels[1].panel_center[1:splitidx, 1]
pdo = body_panels[1].panel_center[splitidx:end, 1]
ph = body_panels[2].panel_center[:, 1]
vsi = body_strengths[1:splitidx]
vso = body_strengths[splitidx:npd]
vsh = body_strengths[(npd + 1):end]
plot!(pg, pdi, abs.(vsi); label="Inner Surface")
plot!(pg, pdo, abs.(vso); label="Outer Surface")
plot!(pg, ph, abs.(vsh); label="Hub Surface")
plot!(
    pg,
    pdi,
    Vinf .* ones(length(pdi));
    linestyle=:dash,
    color=:black,
    label="Vinf for Reference",
)

savefig(pgeom, project_dir * "/dev_debug_archive/duct_only/bodygeomforvx.pdf")
savefig(pv, project_dir * "/dev_debug_archive/duct_only/bodyvx.pdf")
savefig(pg, project_dir * "/dev_debug_archive/duct_only/bodyvs.pdf")
