project_dir = dirname(dirname(dirname(dirname(@__FILE__))))
if project_dir == ""
    project_dir = "."
end

using DuctAPE
const dt = DuctAPE

datapath = project_dir * "/test/data/dfdc_init_iter1/"
savepath = project_dir * "/test/manual_tests/dfdc_comparisons/"

include(project_dir * "/visualize/plots_default.jl")

# read in DFDC data
include(datapath * "reformat_dfdc_data.jl")
include(datapath * "ductape_parameters.jl")
include(datapath * "ductape_formatted_dfdc_geometry.jl")

# generate inputs
inputs = dt.precomputed_inputs(
    system_geometry,
    rotorstator_parameters, #vector of named tuples
    freestream,
    reference_parameters;
    debug=false,
)

gr()
#for w2b in 1:size(w2baic, 1)
#    # PLOT WAKE AIC'S
#    println("wake-body ", w2b, " of ", size(w2baic, 1))
#    pa = plot(
#        1:size(w2baic, 2),
#        inputs.A_bw[w2b, :];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="wake node influence",
#        # size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(w2baic, 2),
#        -w2baic[w2b, 1:end];
#        label="-DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        # size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        # size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, w2b]],
#        [inputs.body_vortex_panels.controlpoint[2, w2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:topleft)
#    savefig(p, savepath * "w2baic" * lpad(w2b, 3, "0") * ".png")
#end

#for r2b in 1:size(r2baic, 1)
#    #PLOT ROTOR-BODY AICS
#    println("rotor-body ", r2b, " of ", size(r2baic, 1))
#    pa = plot(
#        1:size(r2baic, 2),
#        inputs.A_br[r2b, :, 1];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="rotor node influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(r2baic, 2),
#        r2baic[r2b, 1:end];
#        label="+DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, r2b]],
#        [inputs.body_vortex_panels.controlpoint[2, r2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:topright)
#    savefig(p, savepath * "r2baic" * lpad(r2b, 3, "0") * ".png")
#end

#for b2r in 1:size(b2rvhat, 1)
#    #PLOT body-rotor V's
#    println("body-rotor-z ", b2r, " of ", size(b2rvhat, 1))
#    pa = plot(
#        1:size(b2rvhat, 2),
#        inputs.v_rb[1][b2r, :, 1]; #z unit velocity
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="body node z influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(b2rvhat, 2),
#        b2rvhat[b2r, 1:end, 1]; #z unit velocity
#        label="DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=mygray,
#        linestyle=:dot,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=mygray,
#        linestyle=:dot,
#    )
#    plot!(
#        pg,
#        inputs.rotor_source_panels.node[1][1, :],
#        inputs.rotor_source_panels.node[1][2, :];
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.rotor_source_panels[1].controlpoint[1, b2r]],
#        [inputs.rotor_source_panels[1].controlpoint[2, b2r]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250),legend=:topright)
#    savefig(p, savepath * "b2rvhatz" * lpad(b2r, 3, "0") * ".png")
#end

#for b2r in 1:size(b2rvhat, 1)
#    #PLOT body-rotor V's
#    println("body-rotor-r ", b2r, " of ", size(b2rvhat, 1))
#    pa = plot(
#        1:size(b2rvhat, 2),
#        inputs.v_rb[1][b2r, :, 2]; #r unit velocity
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="body node r influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(b2rvhat, 2),
#        b2rvhat[b2r, 1:end, 2]; #r unit velocity
#        label="DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=mygray,
#        linestyle=:dot,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=mygray,
#        linestyle=:dot,
#    )
#    plot!(
#        pg,
#        inputs.rotor_source_panels.node[1][1, :],
#        inputs.rotor_source_panels.node[1][2, :];
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.rotor_source_panels[1].controlpoint[1, b2r]],
#        [inputs.rotor_source_panels[1].controlpoint[2, b2r]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250),legend=:topright)
#    savefig(p, savepath * "b2rvhatr" * lpad(b2r, 3, "0") * ".png")
#end


#for b2b in 1:size(b2bvhat, 1)
#    #PLOT body-body V's
#    println("body-body-z ", b2b, " of ", size(b2bvhat, 1))
#    pa = plot(
#        1:size(b2bvhat, 2),
#        inputs.v_bb[b2b, :, 1];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="body node z influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(b2bvhat, 2),
#        -b2bvhat[b2b, 1:end, 1];
#        label="-DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, b2b]],
#        [inputs.body_vortex_panels.controlpoint[2, b2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
#    savefig(p, savepath * "b2bvhatz" * lpad(b2b, 3, "0") * ".png")
#end

#for b2b in 1:size(b2bvhat, 1)
#    #PLOT body-body V's
#    println("body-body-r ", b2b, " of ", size(b2bvhat, 1))
#    pa = plot(
#        1:size(b2bvhat, 2),
#        inputs.v_bb[b2b, :, 2];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="body node r influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(b2bvhat, 2),
#        -b2bvhat[b2b, 1:end, 2];
#        label="-DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, b2b]],
#        [inputs.body_vortex_panels.controlpoint[2, b2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
#    savefig(p, savepath * "b2bvhatr" * lpad(b2b, 3, "0") * ".png")
#end

#for r2b in 1:size(r2bvhat, 1)
#    #PLOT rotor-body V's
#    println("rotor-body-z ", r2b, " of ", size(r2bvhat, 1))
#    pa = plot(
#        1:size(r2bvhat, 2),
#        inputs.v_br[1][r2b, :, 1];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="rotor node z influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(r2bvhat, 2),
#        r2bvhat[r2b, 1:end, 1];
#        label="+DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
# pg,
# inputs.body_vortex_panels.node[1, 108:end],
# inputs.body_vortex_panels.node[2, 108:end];
# aspectratio=1,
# size=(300, 150),
# label="",
# color=:black,
# )
# plot!(
# pg,
# [inputs.body_vortex_panels.controlpoint[1, r2b]],
# [inputs.body_vortex_panels.controlpoint[2, r2b]];
# seriestype=:scatter,
# color=2,
# label="",
# )
# p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
# savefig(p, savepath * "r2bvhatz" * lpad(r2b, 3, "0") * ".png")
# end

#for r2b in 1:size(r2bvhat, 1)
#    #PLOT rotor-body V's
#    println("rotor-body-r ", r2b, " of ", size(r2bvhat, 1))
#    pa = plot(
#        1:size(r2bvhat, 2),
#        inputs.v_br[1][r2b, :, 2];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="rotor node r influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(r2bvhat, 2),
#        r2bvhat[r2b, 1:end, 2];
#        label="+DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
# plot!(
#     pg,
#     [inputs.body_vortex_panels.controlpoint[1, r2b]],
#     [inputs.body_vortex_panels.controlpoint[2, r2b]];
#     seriestype=:scatter,
#     color=2,
#     label="",
# )
# p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
# savefig(p, savepath * "r2bvhatr" * lpad(r2b, 3, "0") * ".png")
# end

#for w2b in 1:size(w2bvhat, 1)
#    #PLOT wake-body V's
#    println("wake-body-z ", w2b, " of ", size(w2bvhat, 1))
#    pa = plot(
#        1:size(w2bvhat, 2),
#        inputs.v_bw[w2b, :, 1];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="wake node z influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(w2bvhat, 2),
#        -w2bvhat[w2b, 1:end, 1];
#        label="-DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, w2b]],
#        [inputs.body_vortex_panels.controlpoint[2, w2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
#    savefig(p, savepath * "w2bvhatz" * lpad(w2b, 3, "0") * ".png")
#end

#for w2b in 1:size(w2bvhat, 1)
#    #PLOT wake-body V's
#    println("wake-body-r ", w2b, " of ", size(w2bvhat, 1))
#    pa = plot(
#        1:size(w2bvhat, 2),
#        inputs.v_bw[w2b, :, 2];
#        label="DuctAPE",
#        seriestype=:scatter,
#        markersize=2,
#        markerstrokewidth=0,
#        xlabel="node id",
#        ylabel="wake node r influence",
#        size=(300, 200),
#    )
#    plot!(
#        pa,
#        1:size(w2bvhat, 2),
#        -w2bvhat[w2b, 1:end, 2];
#        label="-DFDC",
#        seriestype=:scatter,
#        markersize=1.25,
#        markerstrokewidth=0,
#        ylim=(-0.05, 0.05),
#    )
#    pg = plot(
#        inputs.body_vortex_panels.node[1, 1:107],
#        inputs.body_vortex_panels.node[2, 1:107];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        inputs.body_vortex_panels.node[1, 108:end],
#        inputs.body_vortex_panels.node[2, 108:end];
#        aspectratio=1,
#        size=(300, 150),
#        label="",
#        color=:black,
#    )
#    plot!(
#        pg,
#        [inputs.body_vortex_panels.controlpoint[1, w2b]],
#        [inputs.body_vortex_panels.controlpoint[2, w2b]];
#        seriestype=:scatter,
#        color=2,
#        label="",
#    )
#    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
#    savefig(p, savepath * "w2bvhatr" * lpad(w2b, 3, "0") * ".png")
#end

for b2w in 1:size(b2wvhat, 1)
    #PLOT body-wake V's
    println("body-wake-z ", b2w, " of ", size(b2wvhat, 1))
    pa = plot(
        1:size(b2wvhat, 2),
        inputs.v_wb[b2w, :, 1];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="body node z influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(b2wvhat, 2),
        -b2wvhat[b2w, 1:end, 1];
        label="-DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
     for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
         plot!(
             inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
             inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
             label="",
             color=:black,
             linewidth=0.5,
         )
     end
     plot!(
         pg,
         [inputs.wake_vortex_panels.controlpoint[1, b2w]],
         [inputs.wake_vortex_panels.controlpoint[2, b2w]];
         seriestype=:scatter,
         color=2,
         label="",
     )
     p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
     savefig(p, savepath * "b2wvhatz" * lpad(b2w, 3, "0") * ".png")
 end

for b2w in 1:size(b2wvhat, 1)
    #PLOT body-wake V's
    println("body-wake-r ", b2w, " of ", size(b2wvhat, 1))
    pa = plot(
        1:size(b2wvhat, 2),
        inputs.v_wb[b2w, :, 2];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="body node r influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(b2wvhat, 2),
        -b2wvhat[b2w, 1:end, 2];
        label="-DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
     for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
         plot!(
             inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
             inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
             label="",
             color=:black,
             linewidth=0.5,
         )
     end
     plot!(
         pg,
         [inputs.wake_vortex_panels.controlpoint[1, b2w]],
         [inputs.wake_vortex_panels.controlpoint[2, b2w]];
         seriestype=:scatter,
         color=2,
         label="",
     )
     p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
     savefig(p, savepath * "b2wvhatr" * lpad(b2w, 3, "0") * ".png")
 end

#TODO: rotor to wake
for r2w in 1:size(r2wvhat, 1)
    #PLOT rotor-wake V's
    println("rotor-wake-z ", r2w, " of ", size(r2wvhat, 1))
    pa = plot(
        1:size(r2wvhat, 2),
        inputs.v_wr[1][r2w, :, 1];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="rotor node z influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(r2wvhat, 2),
        r2wvhat[r2w, 1:end, 1];
        label="DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )

    for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
        plot!(
            inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
            inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
            label="",
            color=:black,
            linewidth=0.5,
        )
    end
    plot!(
        pg,
        [inputs.wake_vortex_panels.controlpoint[1, r2w]],
        [inputs.wake_vortex_panels.controlpoint[2, r2w]];
        seriestype=:scatter,
        color=2,
        label="",
    )
    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
    savefig(p, savepath * "r2wvhatz" * lpad(r2w, 3, "0") * ".png")
end

for r2w in 1:size(r2wvhat, 1)
    #PLOT rotor-wake V's
    println("rotor-wake-r ", r2w, " of ", size(r2wvhat, 1))
    pa = plot(
        1:size(r2wvhat, 2),
        inputs.v_wr[1][r2w, :, 2];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="rotor node r influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(r2wvhat, 2),
        r2wvhat[r2w, 1:end, 2];
        label="DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )

    for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
        plot!(
            inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
            inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
            label="",
            color=:black,
            linewidth=0.5,
        )
    end
    plot!(
        pg,
        [inputs.wake_vortex_panels.controlpoint[1, r2w]],
        [inputs.wake_vortex_panels.controlpoint[2, r2w]];
        seriestype=:scatter,
        color=2,
        label="",
    )
    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
    savefig(p, savepath * "r2wvhatr" * lpad(r2w, 3, "0") * ".png")
end

#TODO: wake to wake
for w2w in 1:size(w2wvhat, 1)
    #PLOT wake-wake V's
    println("wake-wake-z ", w2w, " of ", size(w2wvhat, 1))
    pa = plot(
        1:size(w2wvhat, 2),
        inputs.v_ww[w2w, :, 1];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="wake node z influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(w2wvhat, 2),
        -w2wvhat[w2w, 1:end, 1];
        label="-DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )

    for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
        plot!(
            inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
            inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
            label="",
            color=:black,
            linewidth=0.5,
        )
    end
    plot!(
        pg,
        [inputs.wake_vortex_panels.controlpoint[1, w2w]],
        [inputs.wake_vortex_panels.controlpoint[2, w2w]];
        seriestype=:scatter,
        color=2,
        label="",
    )
    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
    savefig(p, savepath * "w2wvhatz" * lpad(w2w, 3, "0") * ".png")
end

for w2w in 1:size(w2wvhat, 1)
    #PLOT wake-wake V's
    println("wake-wake-r ", w2w, " of ", size(w2wvhat, 1))
    pa = plot(
        1:size(w2wvhat, 2),
        inputs.v_ww[w2w, :, 2];
        label="DuctAPE",
        seriestype=:scatter,
        markersize=2,
        markerstrokewidth=0,
        xlabel="node id",
        ylabel="wake node r influence",
        size=(300, 200),
    )
    plot!(
        pa,
        1:size(w2wvhat, 2),
        -w2wvhat[w2w, 1:end, 2];
        label="-DFDC",
        seriestype=:scatter,
        markersize=1.25,
        markerstrokewidth=0,
        ylim=(-0.05, 0.05),
    )
    pg = plot(
        inputs.body_vortex_panels.node[1, 1:107],
        inputs.body_vortex_panels.node[2, 1:107];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )
    plot!(
        pg,
        inputs.body_vortex_panels.node[1, 108:end],
        inputs.body_vortex_panels.node[2, 108:end];
        aspectratio=1,
        size=(300, 150),
        label="",
        color=mygray,
        linestyle=:dot,
    )

    for wid in eachcol(inputs.wake_vortex_panels.endnodeidxs)
        plot!(
            inputs.wake_vortex_panels.node[1, wid[1]:wid[2]],
            inputs.wake_vortex_panels.node[2, wid[1]:wid[2]];
            label="",
            color=:black,
            linewidth=0.5,
        )
    end
    plot!(
        pg,
        [inputs.wake_vortex_panels.controlpoint[1, w2w]],
        [inputs.wake_vortex_panels.controlpoint[2, w2w]];
        seriestype=:scatter,
        color=2,
        label="",
    )
    p = plot(pa, pg; margin=5mm, layout=(1, 2), size=(650, 250), legend=:bottomleft)
    savefig(p, savepath * "w2wvhatr" * lpad(w2w, 3, "0") * ".png")
end
