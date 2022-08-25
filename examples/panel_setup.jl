include("geometry_setup.jl")

function setup_panels(; plotpanels=false)

    # get initial geometry and wake grid
    ductgeometry, ductsplines, rotors, wakegrid = initialize_geometry()

    # get paneling of various objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = DuctTAPE.generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    # plot panels
    if plotpanels
        plot(; xlabel="x", ylabel="r", aspectratio=:equal, legend=true, label="")

        # wall panels:
        for i in 1:length(wall_panels.panel_edges_x)
            plot!(
                [wall_panels.panel_edges_x[i][1]; wall_panels.panel_edges_x[i][2]],
                [wall_panels.panel_edges_r[i][1]; wall_panels.panel_edges_r[i][2]];
                color=1,
                linewidth=0.5,
                markershape=:diamond,
                markersize=1,
                label="",
            )
        end

        scatter!(
            getindex.(wall_panels.panel_centers, 1),
            getindex.(wall_panels.panel_centers, 2);
            color=1,
            markersize=1,
            markershape=:circle,
            label="wall panel centers",
        )

        #hub panels:
        for i in 1:length(hub_panels.panel_edges_x)
            plot!(
                [hub_panels.panel_edges_x[i][1]; hub_panels.panel_edges_x[i][2]],
                [hub_panels.panel_edges_r[i][1]; hub_panels.panel_edges_r[i][2]];
                markersize=1,
                markershape=:diamond,
                color=2,
                linewidth=0.5,
                label="",
            )
        end

        scatter!(
            getindex.(hub_panels.panel_centers, 1),
            getindex.(hub_panels.panel_centers, 2);
            markersize=1,
            markershape=:circle,
            color=2,
            label="hub panel centers",
        )

        # println("wpc: ", wake_panels.panel_centers)
        #vortex sheet panels
        for i in 1:length(wake_panels.panel_centers)
            plot!(
                [wake_panels.panel_edges_x[i][1]; wake_panels.panel_edges_x[i][2]],
                [wake_panels.panel_edges_r[i][1]; wake_panels.panel_edges_r[i][2]];
                markersize=1,
                markershape=:diamond,
                color=3,
                linewidth=0.5,
                label="",
            )
        end

        scatter!(
            getindex.(wake_panels.panel_centers, 1),
            getindex.(wake_panels.panel_centers, 2);
            markersize=1,
            markershape=:circle,
            color=3,
            label="vortex sheet panel centers",
        )

        #rotor source panels:
        for i in 1:length(rotor_source_panels.panel_centers)
            plot!(
                [
                    rotor_source_panels.panel_edges_x[i][1]
                    rotor_source_panels.panel_edges_x[i][2]
                ],
                [
                    rotor_source_panels.panel_edges_r[i][1]
                    rotor_source_panels.panel_edges_r[i][2]
                ];
                markersize=1,
                markershape=:diamond,
                color=4,
                linewidth=0.5,
                label="",
            )
        end

        scatter!(
            getindex.(rotor_source_panels.panel_centers, 1),
            getindex.(rotor_source_panels.panel_centers, 2);
            markersize=1,
            markershape=:circle,
            color=4,
            label="rotor source panel centers",
        )

        savefig("examples/test_panels.pdf")
    end

    return nothing
end
