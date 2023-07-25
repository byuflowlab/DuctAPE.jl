#---------------------------------#
#           Initial Bits          #
#---------------------------------#
project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

using DuctTAPE
const dt = DuctTAPE
using FLOWMath
const fm = FLOWMath
using Statistics

datapath = project_dir * "/examples/dfdc_partial/"
savepath = datapath * "outputs/"

include(project_dir * "/visualize/visualize_geometry.jl")
include(project_dir * "/visualize/visualize_flowfield.jl")
include(project_dir * "/visualize/plots_default_new.jl")
# pyplot()

function checkvdotn(V, n)
    vdn = dt.dot.(eachrow(V), eachrow(n))
    println("\tMax error V dot n = ", maximum(abs.(vdn)))
    println("\tRMS error V dot n = ", sqrt.(sum(vdn .^ 2) / length(vdn)))
    println("\tMedian error V dot n = ", median(abs.(vdn)))
    return vdn
end

function load_smooth_duct(npanref, rtransform)

    # load smooth data
    include(datapath * "nacasmoothgeom.jl")
    # get it in the right order for DuctTAPE
    revductcoords = reverse(smoothnormduct; dims=1)
    # move the coordinates into postiion
    revductcoords[:, 2] .+= rtransform
    #overwrite duct coordinates
    duct_coordinates = dt.repanel_airfoil(revductcoords; normalize=false, N=npanref)

    #write for DFDC
    write_coordinates(
        reverse(duct_coordinates; dims=1);
        varname="lewis_duct_coordinates_smooth",
        savepath=savepath,
        filename="smoothed_duct_$(npanref)_panels.jl",
    )

    return duct_coordinates
end

function dtrepanel(duct_coordinates, hub_coordinates, dtpane, xrotor, wake_length=1.0)

    # use coupled repanling for duct
    xwake, _, _, _ = dt.discretize_wake(
        duct_coordinates,
        hub_coordinates,
        xrotor, #xrotor
        wake_length, #wake length
        dtpane[2:end],
        # [20, 20, 40], #npanels
    )
    # - Repanel Bodies - #
    duct_coordinates, hub_coordinates = dt.update_body_geometry(
        duct_coordinates,
        hub_coordinates,
        xwake,
        dtpane[1], #nhub inlet
        dtpane[1]; #nduct inlet
        finterp=FLOWMath.akima,
    )
    return duct_coordinates, hub_coordinates
end

include(project_dir * "/test/data/naca_662-015.jl")
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")

#---------------------------------#
#             Function            #
#---------------------------------#
"""
"""
function run_isolated_geometry(options, duct_coordinates_input, hub_coordinates_input)
    (; geomtype, comptype, geomsource, dtpane, Vinf, npanref) = options

    gc = geomtype * comptype
    gcg = geomtype * comptype * geomsource

    println("Initializing Plots and Plotting EXP Data")
    # - Initialize Plots - #
    # surface velocity
    pv = plot(; xlabel="x", ylabel=L"V_s", ylim=(0.0, 30.0))
    # surface pressure
    pc = plot(; xlabel="x", ylabel=L"c_p", ylim=(-1.25, 1.0), yflip=true)

    ## -- PLOT EXP DATA -- ##
    # - and select dfdc data file while you're at it - #
    if gc == "ld"
        # - get DFDC data - #
        if geomsource == "s"
            dfdcfile = "lewis_duct_smooth/"
        else
            dfdcfile = "lewis_duct/"
        end

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
        # - get DFDC data - #
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
        # - get DFDC data - #
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
        # - get DFDC data - #
        dfdcfile = "lewis_coupled/"
    else
        @error "no DFDC outputs exist for option $gc, need to run it first"
    end

    println("Loading and Plotting DFDC Data")
    ## -- GET DFDC DATA -- ##
    include(datapath * "dfdc_outs/" * dfdcfile * "DFDC_CPS.jl")

    # hub data
    hubx = dfdc_hub_cp[:, 1]
    hubr = dfdc_hub_cp[:, 2]
    hubvs = dfdc_hub_cp[:, end - 1]
    hubcp = dfdc_hub_cp[:, 4]

    # - duct data
    ductx = dfdc_duct_cp[:, 1]
    # println("ndfdc panels: ", length(ductx))
    ductr = dfdc_duct_cp[:, 2]
    ductvs = dfdc_duct_cp[:, end - 1]
    ductcp = dfdc_duct_cp[:, 4]

    ## -- PLOT DFDC DATA -- ##
    if comptype == "d"
        plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
        plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
    elseif comptype == "h"
        plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
        plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
    else
        plot!(pv, hubx, hubvs; linestyle=:dash, color=myred[2], label="DFDC Hub")
        plot!(pc, hubx, hubcp; linestyle=:dash, color=myred[2], label="DFDC Hub")
        plot!(pv, ductx, ductvs; linestyle=:dash, color=myblue[2], label="DFDC Duct")
        plot!(pc, ductx, ductcp; linestyle=:dash, color=myblue[2], label="DFDC Duct")
    end

    println("Setting up DuctTAPE Stuff")
    ## -- SET UP DUCTTAPE DATA -- ##

    if comptype == "d"
        # - set up duct data according to options
        if geomsource == "s"
            duct_coordinates = load_smooth_duct(npanref, duct_coordinates_input[1, 2])
        else
            duct_coordinates = duct_coordinates_input
        end
        hub_coordinates = hub_coordinates_input

        #if using ducttape repanling
        if dtpane != nothing
            duct_coordinates, hub_coordinates = dtrepanel(
                duct_coordinates, hub_coordinates, dtpane, 0.5
            )
        end

        # - set up prescribed panels and coordinates
        _, leid = findmin(duct_coordinates[:, 1])
        coordinates = [duct_coordinates]
        # prescribedpanels = [(leid, 0.0)]
        prescribedpanels = [(1, 0.0)]

    elseif comptype == "h"
        coordinates = [hub_coordinates]
        prescribedpanels = [(length(hub_coordinates) - 1, 0.0)]
        # prescribedpanels = [(1, 0.0)]
    else
        if geomsource == "s"
            duct_coordinates = load_smooth_duct(npanref, duct_coordinates_input[1, 2])
        else
            duct_coordinates = duct_coordinates_input
        end

        hub_coordinates = dt.repanel_revolution(
            hub_coordinates_input; N=ceil(Int, npanref / 2), normalize=false
        )

        #if using ducttape repanling
        if dtpane != nothing
            duct_coordinates, hub_coordinates = dtrepanel(
                duct_coordinates, hub_coordinates_input, dtpane, 0.5
            )
        end

        #TODO: need to fix why hub get's repaneled into having a negative last panel...
        if hub_coordinates[end,2] <= 0.0
            hub_coordinates = hub_coordinates[1:end-1,:]
        end

        println("min hubr: ",hub_coordinates[end,2])

        _, leid = findmin(duct_coordinates[:, 1])
        coordinates = [duct_coordinates, hub_coordinates]
        prescribedpanels = [(leid, 0.0); (size(duct_coordinates, 1), 0.0)]
    end

    # - Panel Geometry for DuctTAPE - #
    panels = dt.generate_panels(coordinates)

    println("Visualizing Paneling")
    # - Visualize paneling - #
    visualize_paneling(;
        body_panels=panels,
        coordinates=[duct_coordinates, hub_coordinates],
        controlpoints=true,
        nodes=true,
        normals=true,
        normal_scaling=0.1,
        savepath=savepath,
        filename=["$gcg-body-geometry.png"; "$gcg-body-geometry.pdf"],
        legendloc=:right,
    )

    println("Solving System")
    ## -- Initialize and solve strenths a la init function -- ##
    # - body to body - #
    A_bb = dt.doublet_panel_influence_matrix(panels.nodes, panels)
    LHS = A_bb

    # - add internal panel stuff - #
    # LHS = dt.doublet_panel_influence_on_internal_panels(A_bb, panels, panels)#, prescribedpanels)

    # right hand side from freestream
    Vinfvec = [Vinf 0.0]
    Vinfmat = repeat(Vinfvec, panels.npanels)
    b_bf = dt.freestream_influence_vector(panels.normal, Vinfmat)
    RHS = b_bf

    # apply Kutta condition
    dt.body_lhs_kutta!(LHS, panels)

    LHSlsq, RHSlsq = dt.prep_leastsquares(LHS, RHS, prescribedpanels)

    mured = LHSlsq \ RHSlsq
    mub = dt.mured2mu(mured, prescribedpanels)

    ## -- Post Process -- ##
    # - Body-induced Surface Velocity - #
    Vb = dt.vfromdoubletpanels(panels.controlpoint, panels.nodes, mub)

    # - "Wake"-induced Surface Velocity - #
    Vb_TE = dt.vfromTE(panels.controlpoint, panels.TEnodes, mub)

    # - ∇μ/2 surface velocity - #
    # Vb_gradmu = dt.vfromgradmu(panels, mub)
    # Vb_gradmu = dt.gradmutry2a(panels, mub)
    # Vb_gradmu = dt.gradmutry2b(panels, mub)
    # Vb_gradmu = dt.vfromgradmutry3a(panels, mub)
    Vb_gradmu = dt.vfromgradmutry3b(panels, mub)

    Vtot = copy(Vinfmat)
    Vtot .+= Vb
    Vtot .+= Vb_TE
    Vtot .+= Vb_gradmu

    # - Get magnitude and split - #
    vs = dt.norm.(eachrow(Vtot))

    # - surface pressure from steady cps - #
    cp = dt.steady_cp(vs, Vinf, Vinf)

    println("Plotting and Saving Outputs")
    ## -- Plot -- ##
    xs = panels.controlpoint[:, 1]
    plot!(pv, xs, vs; label="DuctTAPE")
    plot!(pc, xs, cp; label="DuctTAPE")

    savefig(pv, savepath * gcg * "-velocity-comp.pdf")
    savefig(pv, savepath * gcg * "-velocity-comp.png")
    savefig(pc, savepath * gcg * "-pressure-comp.pdf")
    savefig(pc, savepath * gcg * "-pressure-comp.png")

    println("Visualizing Flowfield")
    # - Visualize Flow field - #
    visualize_flowfield(
        Vinf; body_panels=panels, mub=mub, save_path=savepath, run_name="$gc-velocity-field"
    )
    visualize_surfaces(; body_panels=panels, run_name=gc)

    #TODO: return plot objects?
    return nothing
end

# # Duct: use ducttape repaneling with smooth input geometry
# options = (;
#     Vinf=20,
#     geomtype="l",
#     comptype="d",
#     geomsource="s",
#     dtpane=2*[40, 30, 20, 40],
#     npanref=200,
# )

# # Duct: dont use ducttape repanling, keep smooth geometry
# options = (;
#     Vinf=20, geomtype="l", comptype="d", geomsource="s", dtpane=nothing, npanref=400
# )

# # Body:  use smooth input for duct
# options = (;
#     Vinf=20, geomtype="l", comptype="b", geomsource="s", dtpane=nothing, npanref=400
# )

# Body: use ducttape repaneling
options = (;
    Vinf=20,
    geomtype="l",
    comptype="b",
    geomsource="s",
    dtpane=[40, 30, 20, 40],
    npanref=200,
)

run_isolated_geometry(options, lewis_duct_coordinates, lewis_hub_coordinates)
