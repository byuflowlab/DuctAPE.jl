function visualize_paneling(panels; coordinates=nothing, control_points=true, nodes=true, TEnodes=true, normals=true, normal_scaling = 0.1, savepath="", filename="paneling.pdf")

    # plot generated panels
    plot(aspectratio=1, xlabel="x", ylabel="r")

    # plot inputs
    if !isnothing(coordinates)
        plot!(coordinates[:,1],coordinates[:,2], color=gray[1], linewidth=0.5, label="input coordinates")
    end

    #plot control points
    if control_points
        plot!(panels.control_point[:,1], panels.control_point[:,2], color=red[2], seriestype=:scatter, markershape=:rect, label="Control Points")
    end

    #plot nodes
    if nodes
        for i in 1:length(panels.length)
            lab= i==1 ? "Nodes" : ""
            plot!(panels.nodes[i,:,1],panels.nodes[i,:,2], label=lab, color=blue[2], seriestype=:scatter)
        end
    end

    # plot trailing edge (wake) nodes
    if TEnodes
        for i in 1:length(panels.TEnodes[:,1,1])
            lab = i == 1 ? "TE Nodes" : ""
            plot!(panels.TEnodes[i,:,1], panels.TEnodes[i,:,2], label=lab, color = red[3], seriestype=:scatter, markersize=3)
        end
    end

    #plot normal
    if normals
        for i in 1:length(panels.length)
            lab= i==1 ? "Normals" : ""
        plot!([0.0;normal_scaling*panels.normal[i,1]].+panels.control_point[i,1],[0.0;normal_scaling*panels.normal[i,2]].+panels.control_point[i,2], label=lab, color=blue[3])
        end
    end

    # save figure
    savefig(savepath*filename)

    return nothing

end
