#---------------------------------#
#           Initial Bits          #
#---------------------------------#
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE
using FLOWMath
using Statistics

datapath = project_dir * "/dev_debug_archive/dfdc_partial/"
savepath = datapath * "outputs/"

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/visualize_flowfield.jl")
include(project_dir * "/visualize/plots_default_new.jl")

function checkvdotn(V, n)
    vdn = dt.dot.(eachrow(V), eachrow(n))
    println("\tMax error V dot n = ", maximum(abs.(vdn)))
    println("\tRMS error V dot n = ", sqrt.(sum(vdn .^ 2) / length(vdn)))
    println("\tMedian error V dot n = ", median(abs.(vdn)))
    return vdn
end

#---------------------------------#
#        Set Up Comparion         #
#---------------------------------#

## -- OPTIONS -- ##
geomtype = "l"
comptype = "d"
plotexp = true
dtrp = false
edgeo = false
smoothgeo = true
closeTE = true
vte = true
gm = true
edkutta = true
simkutta = false

gc = geomtype * comptype

# - Initialize Plots - #
pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 30.0))
pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-1.25, 1.0), yflip=true)

include(project_dir * "/test/data/naca_662-015.jl")

include(project_dir * "/test/data/bodyofrevolutioncoords.jl")
if gc == "ld"
    # - get DFDC data - #
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
if closeTE
    duct_coordinates[1, :] .= duct_coordinates[end, :]
end
hub_coordinates = reverse([hubx hubr]; dims=1)

# - Genrate Panels - #
if comptype == "d"
    if edgeo
        include(project_dir * "/dev_debug_archive/new_panels/ed_duct_geom.jl")
    elseif smoothgeo
        include(datapath * "nacasmoothgeom.jl")
        revductcoords = reverse(smoothnormduct; dims=1)
        revductcoords[:, 2] .+= duct_coordinates[1, 2]
        duct_coordinates = dt.repanel_airfoil(revductcoords; normalize=false, N=100)
        # nanid = findall(x->isnan(x), myductcoords)[1][1]
        # myductcoordsraw = myductcoords[1:end .∉ nanid, :]

    end

    if dtrp
        # use coupled repanling for duct
        xwake, _, _, _ = dt.discretize_wake(
            duct_coordinates,
            hub_coordinates,
            0.5, #xrotor
            1.0, #wake length
            [20, 20, 40], #npanels
        )

        # - Repanel Bodies - #
        duct_coordinates, _ = dt.update_body_geometry(
            duct_coordinates,
            hub_coordinates,
            xwake,
            20, #nhub inlet
            20; #nduct inlet
            finterp=FLOWMath.akima,
        )
    end

    _, leid = findmin(duct_coordinates[:, 1])
    coordinates = [duct_coordinates]
    # prescribedpanels = [(leid, 0.0)]
    prescribedpanels = [(1, 0.0)]

    plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
    plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
elseif comptype == "h"
    coordinates = [hub_coordinates]
    prescribedpanels = [(length(hub_coordinates) - 1, 0.0)]
    # prescribedpanels = [(1, 0.0)]

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

# xtest = [2.0; 1.0; 0.0; 1.0; 2.0]
# rtest = [1.0; 0.8; 1.0; 1.2; 1.0]
# coordinates = [xtest rtest]

panels = dt.generate_panels(coordinates)
# panels.itcontrolpoint[1] = 0.5

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
LHS = A_bb

# - add internal panel stuff - #
# LHS = dt.doublet_panel_influence_on_internal_panels(A_bb, panels, panels)#, prescribedpanels)

# right hand side from freestream
Vinf = 30.0
Vinfvec = [Vinf 0.0]
Vinfmat = repeat(Vinfvec, panels.totpanel)
b_bf = dt.freestream_influence_vector(panels.normal, Vinfmat)
RHS = b_bf

# apply Kutta condition
if edkutta
    println("using Ed's kutta condition")
    dt.body_lhs_kutta!(LHS, panels)
end

if simkutta
    println("using simpler kutta condition")
    LHS = [A_bb; zeros(size(A_bb, 2))']
    LHS[end, 1] = LHS[end, end] = 1.0
    RHS = [b_bf; 0.0]
end

# - add internal panel stuff - #
# RHS = dt.freestream_influence_on_internal_panels(b_bf, panels, Vinfvec)

if closeTE
    # - From state initialziation - #
    # mub = dt.solve_body_strengths(LHS, RHS, prescribedpanels, nothing)
    # mub = dt.solve_body_strengths(LHS, RHS, prescribedpanels, panels.nbodies)

    LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

    mured = LHSlsq \ RHSlsq

    println("Inf norm error of Least Squares solver: ", maximum(LHSlsq * mured .- RHSlsq))

    # mured2mu!(mub, mured[1:(end - nbodies)], prescribedpanels)
    mub = dt.mured2mu(mured, prescribedpanels)

else
    if simkutta
        LHSlsq = LHS' * LHS
        RHSlsq = LHS' * RHS
        mub = LHSlsq \ RHSlsq
        mured = mub
        # mub=mulsq[1:end-1]
    else
        mub = LHS \ RHS
        println("Inf norm error of Direct solver (open TE): ", maximum(LHS * mub .- RHS))
    end
end
## -- Post Process -- ##
# - Body-induced Surface Velocity - #
Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mub)

#### ---- SANITY CHECK ---- ####
lm = LHS * mub
lms = LHSlsq * mured
println("LHSlsq*mublsq - RHSlsq max error: ", maximum(abs.(lms .- RHSlsq)))
println("LHS*mub - RHS max error: ", maximum(abs.(lm .- RHS)))

# - "Wake"-induced Surface Velocity - #
Vb_TE = dt.vfromTE(panels.controlpoint, panels.TEnodes, mub)

# - ∇μ/2 surface velocity - #
# Vb_gradmu = dt.vfromgradmu(panels, mub)
# Vb_gradmu = dt.gradmutry2a(panels, mub)
# Vb_gradmu = dt.gradmutry2b(panels, mub)
# Vb_gradmu = dt.vfromgradmutry3a(panels, mub)
Vb_gradmu = dt.vfromgradmutry3b(panels, mub)
println("check if gradmu is tangent")
checkvdotn(Vb_gradmu, panels.normal)

#theoretically, Vtot dot nhat should be zero and Vtot dot that should be norm(Vtot)
Vtot = copy(Vinfmat)
println("Vinf only")
vdn = checkvdotn(Vtot, panels.normal);
Vtot .+= Vb
println("After adding V_b")
vdn = checkvdotn(Vtot, panels.normal)

if vte
    Vtot .+= Vb_TE
    println("After adding V_TE")
    vdn = checkvdotn(Vtot, panels.normal)
end

if gm
    Vtot .+= Vb_gradmu
    println("After adding ∇μ")
    vdngm = checkvdotn(Vtot, panels.normal)
end

# - Get magnitude and split - #
vs = dt.norm.(eachrow(Vtot))

# - surface pressure from steady cps - #
cp = dt.steady_cp(vs, Vinf, Vinf)

## -- Plot -- ##
xs = panels.controlpoint[:, 1]
plot!(pv, xs, vs; label="DuctAPE")
plot!(pc, xs, cp; label="DuctAPE")

savefig(pv, savepath * gc * "velocity-comp.pdf")
savefig(pc, savepath * gc * "pressure-comp.pdf")

### -- Check Internal Panel -- ##
## - Body-induced Surface Velocity - #
#Vbit = dt.vfromdoubletpanels(panels.itcontrolpoint, panels.nodes, mub)

## - "Wake"-induced Surface Velocity - #
#Vb_TEit = dt.vfromTE(panels.itcontrolpoint, panels.TEnodes, mub)

##theoretically, Vtotit should be zeros
#Vtotit = [Vinf 0.0] .+ Vbit .+ Vb_TEit

#println("Vtot on internal panel = ", Vtotit)

# - Visualize Flow field - #
visualize_flowfield(
    Vinf; body_panels=panels, mub=mub, save_path=savepath, run_name="$gc-velocity-field"
)
visualize_surfaces(; body_panels=panels, run_name=gc)
