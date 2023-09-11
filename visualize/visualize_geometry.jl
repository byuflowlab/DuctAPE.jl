function visualize_paneling(;
    body_panels=nothing,
    rotor_panels=nothing,
    wake_panels=nothing,
    coordinates=nothing,
    controlpoints=true,
    nodes=true,
    wakeinterfaceid=[],
    prescribedpanels=nothing,
    TEnodes=true,
    normals=true,
    normal_scaling=0.1,
    savepath="",
    filename=["paneling.pdf"],
    legendloc=:best,
    zoom=false,
    limits=nothing,
    nodemarkersize=2,
    cpmarkersize=2,
)

    ## -- Initialize Plot -- ##
    # plot generated body_panels
    if limits != nothing
        p = plot(;
            ylim=limits.ylim,
            xlim=limits.xlim,
            aspectratio=1,
            xlabel="x",
            ylabel="r",
            legend=legendloc,
        )
    else
        p = plot(; aspectratio=1, xlabel="x", ylabel="r", legend=legendloc)
    end

    if zoom
        bcpsize = 6
        bnsize = 6
        rnsize = 5
        rcpsize = 4
        wcpsize = 3
        wnsize = 3
        tesize = 8
        psize = 3
        ifsize = 2

    else
        bcpsize = cpmarkersize
        bnsize = nodemarkersize
        rnsize = nodemarkersize
        rcpsize = cpmarkersize - 0.5
        wcpsize = cpmarkersize - 0.5
        wnsize = nodemarkersize - 0.5
        tesize = 3
        psize = 3
        ifsize = 1.5
    end

    #---------------------------------#
    #         Plot Body body_Panels        #
    #---------------------------------#
    if !isnothing(body_panels)

        # plot inputs
        if !isnothing(coordinates)
            for (ic, coords) in enumerate(coordinates)
                blab = ic == 1 ? "Input Body Coordinates" : ""
                plot!(
                    p,
                    coords[:, 1],
                    coords[:, 2];
                    color=mygray[1],
                    linewidth=cpmarkersize / 4,
                    label=blab,
                )
            end
        end

        # plot trailing edge (wake) nodes
        if TEnodes
            for i in 1:length(body_panels.TEnodes)
                lab = i == 1 ? "Body TE Nodes" : ""
                plot!(
                    p,
                    [body_panels.TEnodes[i].pos[1]],
                    [body_panels.TEnodes[i].pos[2]];
                    label=lab,
                    color=mygray[1],
                    seriestype=:scatter,
                    markersize=tesize,
                )
            end
        end

        #plot nodes
        if nodes
            for i in 1:length(body_panels.len)
                lab = i == 1 ? "Body Nodes" : ""
                plot!(
                    p,
                    body_panels.nodes[i, :, 1],
                    body_panels.nodes[i, :, 2];
                    label=lab,
                    color=myblue[2],
                    seriestype=:scatter,
                    markersize=bnsize,
                )
            end
        end

        #plot control points
        if controlpoints
            plot!(
                p,
                body_panels.controlpoint[:, 1],
                body_panels.controlpoint[:, 2];
                color=myblue[3],
                seriestype=:scatter,
                markershape=:rect,
                markersize=bcpsize,
                label="Body Control Points",
            )
        end

        if !isnothing(prescribedpanels)
            for ipp in 1:length(prescribedpanels)
                lab = ipp == 1 ? "Prescribed Panel(s)" : ""
                plot!(
                    p,
                    [body_panels.controlpoint[prescribedpanels[ipp][1], 1]],
                    [body_panels.controlpoint[prescribedpanels[ipp][1], 2]];
                    seriestype=:scatter,
                    color=mygreen[2],
                    markersize=psize,
                    markershape=:diamond,
                    label=lab,
                )
            end
        end

        #plot normal
        if normals
            for i in 1:length(body_panels.len)
                lab = i == 1 ? "Body Normals" : ""
                plot!(
                    p,
                    [0.0; normal_scaling * body_panels.normal[i, 1]] .+
                    body_panels.controlpoint[i, 1],
                    [0.0; normal_scaling * body_panels.normal[i, 2]] .+
                    body_panels.controlpoint[i, 2];
                    label=lab,
                    linewidth=cpmarkersize / 4,
                    color=myblue[1],
                )
            end
        end
    end

    #---------------------------------#
    #        Plot Rotor Panels        #
    #---------------------------------#
    if !isnothing(rotor_panels)
        for irotor in 1:length(rotor_panels)

            #plot control points
            if controlpoints
                plot!(
                    p,
                    rotor_panels[irotor].controlpoint[:, 1],
                    rotor_panels[irotor].controlpoint[:, 2];
                    color=mygray[irotor],
                    seriestype=:scatter,
                    markersize=rcpsize,
                    label="Rotor $irotor Nodes",
                )
            end

            #plot nodes
            if nodes
                for i in 1:length(rotor_panels[irotor].len)
                    lab = i == 1 ? "Rotor $irotor Panel Edges" : ""
                    plot!(
                        p,
                        rotor_panels[irotor].nodes[i, :, 1],
                        rotor_panels[irotor].nodes[i, :, 2];
                        label=lab,
                        color=mygray[irotor],
                        seriestype=:scatter,
                        markershape=:diamond,
                        markersize=rnsize,
                    )
                end
            end
        end
    end

    #---------------------------------#
    #         Plot Wake Panels        #
    #---------------------------------#
    if !isnothing(wake_panels)
        #plot nodes
        if nodes
            for i in 1:length(wake_panels.len)
                lab = i == 1 ? "Wake Panel Edges" : ""
                plot!(
                    p,
                    wake_panels.nodes[i, :, 1],
                    wake_panels.nodes[i, :, 2];
                    label=lab,
                    color=myred[2],
                    seriestype=:scatter,
                    markershape=:diamond,
                    markersize=wnsize,
                )
            end
        end

        #plot control points
        if controlpoints
            plot!(
                p,
                wake_panels.controlpoint[:, 1],
                wake_panels.controlpoint[:, 2];
                color=myred[3],
                seriestype=:scatter,
                markersize=wcpsize,
                label="Wake Nodes",
            )

            if !isempty(wakeinterfaceid)
                plot!(
                    p,
                    wake_panels.controlpoint[wakeinterfaceid, 1],
                    wake_panels.controlpoint[wakeinterfaceid, 2];
                    color=mygreen[2],
                    seriestype=:scatter,
                    markershape=:utriangle,
                    markersize=ifsize,
                    label="Body Interface Points",
                )
            end
        end
    end

    # save figure
    for fn in filename
        savefig(savepath * fn)
    end

    return nothing
end

function write_coordinates(
    coordinates; varname="coords", savepath="", filename="coordinates.jl"
)
    f = open(savepath * filename, "w")
    write(f, varname * " = [\n")

    for xr in eachrow(coordinates)
        write(f, "$(xr[1]) $(xr[2])\n")
    end

    write(f, "]")
    close(f)

    return nothing
end
