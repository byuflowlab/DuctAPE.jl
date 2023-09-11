function plot_rotor_geometry(inputs, savepath; num=false, offset=0.0075)
    nwp = inputs.wake_vortex_panels.npanels
    wcp = inputs.wake_vortex_panels.controlpoint
    nrp = inputs.rotor_source_panels[1].npanels
    rcp = inputs.rotor_source_panels[1].controlpoint

    xt = unique(wcp[:, 1])
    rt = unique(wcp[:, 2])
    plot(;
        aspcectratio=1,
        xlabel="x",
        ylabel="r",
        ylim=(0.0, 0.15),
        xlim=(0.0, 0.25),
        xticks=xt,
        yticks=rt,
    )
    plot!(wcp[:, 1], wcp[:, 2]; seriestype=:scatter, label="wake cps", color=myblue[2], markersize=2)

    plot!(rcp[:, 1], rcp[:, 2]; seriestype=:scatter, label="rotor cps", color=myred[2], markersize=2)

    if num
        for i in 1:nwp
            annotate!(wcp[i, 1] + offset, wcp[i, 2] + offset, text("$i", color=myblue[2]))
        end
        for i in 1:nrp
            annotate!(rcp[i, 1] + offset, rcp[i, 2] + offset, text("$i", color=myred[2]))
        end
    end

    savefig(savepath * "rotor-wake-geometry.pdf")
    return nothing
end
