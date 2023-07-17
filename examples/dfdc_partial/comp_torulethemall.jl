#---------------------------------#
#           Initial Bits          #
#---------------------------------#
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctTAPE
const dt = DuctTAPE

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath * "outputs/"

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/plots_default_new.jl")

#---------------------------------#
#        Set Up Comparion         #
#---------------------------------#

## -- OPTIONS -- ##
geomtype = "l"
comptype = "d"
plotexp = true

gc = geomtype * comptype

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 30.0))
pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-1.25, 1.0), yflip=true)

include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
# - get DFDC data - #
if gc == "ld"
    dfdcfile = "lewis_duct/"
    # actually generate the panels

    plot!(
        pc,
        pressurexupper,
        pressureupper;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:utriangle,
        markersize=3,
        label="exp duct outer",
    )
    plot!(
        pc,
        pressurexlower,
        pressurelower;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:dtriangle,
        markersize=3,
        label="exp duct inner",
    )

elseif gc == "lh"
    dfdcfile = "lewis_hub/"
    plot!(
        pv,
        Vs_over_Vinf_x,
        Vs_over_Vinf_vs * Vinf;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:utriangle,
        markersize=3,
        label="hub experimental",
    )
elseif gc == "lb"
    dfdcfile = "lewis_body/"
    plot!(
        pc,
        pressurexupper,
        pressureupper;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:utriangle,
        markersize=3,
        label="exp duct outer",
    )
    plot!(
        pc,
        pressurexlower,
        pressurelower;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:dtriangle,
        markersize=3,
        label="exp duct inner",
    )
    plot!(
        pv,
        Vs_over_Vinf_x,
        Vs_over_Vinf_vs * Vinf;
        seriestype=:scatter,
        color=myblue[1],
        markershape=:utriangle,
        markersize=3,
        label="hub experimental",
    )
elseif gc == "lc"
    dfdcfile = "lewis_coupled/"
else
    @error "no DFDC outputs exist for option $gc, need to run it first"
end

include(datapath * "dfdc_outs/" * dfdcfile * "DFDC_CPS.jl")

hubx = dfdc_hub_cp[:, 1]
hubr = dfdc_hub_cp[:, 2]
hubvs = dfdc_hub_cp[:, end - 1]
hubcp = dfdc_hub_cp[:, 4]

ductx = dfdc_duct_cp[:, 1]
ductr = dfdc_duct_cp[:, 2]
ductvs = dfdc_duct_cp[:, end - 1]
ductcp = dfdc_duct_cp[:, 4]

duct_coordinates = reverse([ductx ductr]; dims=1)
duct_coordinates[1,:] .= duct_coordinates[end,:]
hub_coordinates = reverse([hubx hubr]; dims=1)

# - Genrate Panels - #
if comptype == "d"
    _, leid = findmin(duct_coordinates[:, 1])
    coordinates = [duct_coordinates]
    prescribedpanels = [(leid, 0.0)]

    plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
    plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
elseif comptype == "h"
    coordinates = [hub_coordinates]
    prescribedpanels = [(1, 0.0)]

    plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
    plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
else
    _, leid = findmin(duct_coordinates[:, 1])
    coordinates = [duct_coordinates, hub_coordinates]
    prescribedpanels = [(leid, 0.0); (1, 0.0)]

    plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
    plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
    plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
    plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
end

panels = dt.generate_panels(coordinates)

# - Visualize paneling - #
visualize_paneling(;
    body_panels=panels,
    coordinates=[duct_coordinates, hub_coordinates],
    controlpoints=true,
    nodes=true,
    normals=true,
    normal_scaling=0.1,
    savepath=savepath,
    filename="$gc-body-geometry.pdf",
    legendloc=:right,
)

## -- Initialize and solve strenths a la init function -- ##
# - body to body - #
A_bb = dt.doublet_panel_influence_matrix(panels.nodes, panels)

# - add internal panel stuff - #
LHS = dt.doublet_panel_influence_on_internal_panels(A_bb, panels, panels)

# apply Kutta condition
dt.body_lhs_kutta!(LHS, panels)

# right hand side from freestream
Vinf = 20.0
Vinfvec = [Vinf 0.0]
Vinfmat = repeat(Vinfvec, panels.npanels)
b_bf = dt.freestream_influence_vector(panels.normal, Vinfmat)

# - add internal panel stuff - #
RHS = dt.freestream_influence_on_internal_panels(b_bf, panels, Vinfvec)

# - From state initialziation - #
mub = dt.solve_body_strengths(LHS, RHS, prescribedpanels, panels.nbodies)

## -- Post Process -- ##
# - Body-induced Surface Velocity - #
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mub)

# - "Wake"-induced Surface Velocity - #
Vb_TE = dt.vfromTE(panels.controlpoint, panels.TEnodes, mub)

# - ∇μ/2 surface velocity - #
Vb_gradmu = dt.vfromgradmu(panels, mub)

#theoretically, Vtot dot nhat should be zero and Vtot dot that should be norm(Vtot)
Vtot = Vinfmat .+ Vb .+ Vb_TE .+ Vb_gradmu

# - Get magnitude and split - #
vs = dt.norm.(eachrow(Vtot))

# - surface pressure from steady cps - #
cp = dt.steady_cp(vs, Vinf, Vinf)

## -- Plot -- ##
xs = panels.controlpoint[:, 1]
plot!(pv, xs, vs; label="DuctTAPE")
plot!(pc, xs, cp; label="DuctTAPE")

savefig(pv, savepath * gc * "velocity-comp.pdf")
savefig(pc, savepath * gc * "pressure-comp.pdf")

## -- Check Internal Panel -- ##
# - Body-induced Surface Velocity - #
Vbit = dt.vfromdoubletpanels(panels.itcontrolpoint, panels.nodes, mub)

# - "Wake"-induced Surface Velocity - #
Vb_TEit = dt.vfromTE(panels.itcontrolpoint, panels.TEnodes, mub)

#theoretically, Vtotit should be zeros
Vtotit = [Vinf 0.0] .+ Vbit .+ Vb_TEit

println("Vtot on internal panel = ", Vtotit)

