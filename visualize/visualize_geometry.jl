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
    filename="paneling.pdf",
    legendloc=:best,
)

    ## -- Initialize Plot -- ##
    # plot generated body_panels
    p = plot(; aspectratio=1, xlabel="x", ylabel="r", legend=legendloc)

    #---------------------------------#
    #         Plot Body body_Panels        #
    #---------------------------------#
    if !isnothing(body_panels)

        # plot inputs
        if !isnothing(coordinates)
            for coords in coordinates
                plot!(
                    p,
                    coords[:, 1],
                    coords[:, 2];
                    color=mygray[1],
                    linewidth=0.5,
                    label="body input coordinates",
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
                    markersize=3,
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
                    markersize=3,
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
                    markershape=:rect,
                    label="Rotor $irotor Control Points",
                )
            end

            #plot nodes
            if nodes
                for i in 1:length(rotor_panels[irotor].len)
                    lab = i == 1 ? "Rotor $irotor Nodes" : ""
                    plot!(
                        p,
                        rotor_panels[irotor].nodes[i, :, 1],
                        rotor_panels[irotor].nodes[i, :, 2];
                        label=lab,
                        color=mygray[irotor],
                        seriestype=:scatter,
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
                lab = i == 1 ? "Wake Nodes" : ""
                plot!(
                    p,
                    wake_panels.nodes[i, :, 1],
                    wake_panels.nodes[i, :, 2];
                    label=lab,
                    color=myred[2],
                    seriestype=:scatter,
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
                markershape=:rect,
                label="Wake Control Points",
            )

            if !isempty(wakeinterfaceid)
                plot!(
                    p,
                    wake_panels.controlpoint[wakeinterfaceid, 1],
                    wake_panels.controlpoint[wakeinterfaceid, 2];
                    color=mygreen[2],
                    seriestype=:scatter,
                    markershape=:utriangle,
                    markersize=1.5,
                    label="Body Interface Points",
                )
            end
        end
    end

    # save figure
    savefig(savepath * filename)

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
