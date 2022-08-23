include("geometry_setup.jl")

function setup_panels()

    # Get initial geometry and wake grid
    ductgeometry, ductsplines, rotors, wakegrid = initialize_geometry(; hubscale=0.75)

    # Get paneling of various objects
    wall_panels, hub_panels, wake_panels, rotor_source_panels = DuctTAPE.generate_paneling(
        ductgeometry, ductsplines, rotors, wakegrid
    )

    # PLOT PANELS

    figure(2, figsize=(6,4))
    clf()

    # println("wallp:")
    # display(wake_panels.panel_centers)
    # display(wake_panels.panel_edges_x)
    # display(wake_panels.panel_edges_r)

    # wall panels:
    for i in 1:length(wall_panels.panel_edges_x)
        plot(wall_panels.panel_edges_x[i], wall_panels.panel_edges_r[i], "+-C0")
    end

    plot(
        getindex.(wall_panels.panel_centers, 1),
        getindex.(wall_panels.panel_centers, 2),
        "oC0";
        label="wall panels",
    )

    #hub panels:
    for i in 1:length(hub_panels.panel_edges_x)
        plot(hub_panels.panel_edges_x[i], hub_panels.panel_edges_r[i], "+-C1")
    end
    plot(
        getindex.(hub_panels.panel_centers, 1),
        getindex.(hub_panels.panel_centers, 2),
        "oC1";
        label="hub panels",
    )

    # println("wpc: ", wake_panels.panel_centers)
    #vortex sheet panels
    for i in 1:length(wake_panels.panel_centers)
        plot(wake_panels.panel_edges_x[i], wake_panels.panel_edges_r[i], "+-C2")
    end
    plot(
        getindex.(wake_panels.panel_centers, 1),
        getindex.(wake_panels.panel_centers, 2),
        "oC2";
        label="vortex sheet panels",
    )

    #rotor source panels:
    for i in 1:length(rotor_source_panels.panel_centers)
        plot(
            rotor_source_panels.panel_edges_x[i],
            rotor_source_panels.panel_edges_r[i],
            "+-C3",
        )
    end
    plot(
        getindex.(rotor_source_panels.panel_centers, 1),
        getindex.(rotor_source_panels.panel_centers, 2),
        "oC3";
        label="rotor source panels",
    )

    axis("equal")
    legend()

    savefig("examples/test_panels.png"; bbox_inches="tight")

    return nothing
end

setup_panels()

