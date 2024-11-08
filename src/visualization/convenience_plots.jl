#---------------------------------#
#         DISPATCH TYPES          #
#---------------------------------#
struct staticPlots end
struct animatedPlots end

"""
    generate_plots(
        ::staticPlots, (or ::animatedPlots)
        Plots,
        ins,
        outs;
        save_path="",
        static_file_type=".pdf",
        suffix=nothing,
        plot_geometry=true,
        plot_pressure=false,
        plot_velocity=false,
        plot_boundary_layer=false,
        plot_streamlines=false,
        verbose=false,
        kwargs...,
    )

Generate standard suite of plots or animations from input and output objects.

# Arguments:
- `::staticPlots (or ::animatedPlots)` :
- `Plots::` : the Plots package namespace
- `ins::NamedTuple` : returned inputs from `analyze` function
- `outs::Vector{NamedTuple}` : retured outputs from `analyze` function

# Keyword Arguments:
- `save_path=""` : custom save path
- `static_file_type=".pdf"` : file type for static files (must be compatible with the desired backend)
- `suffix=nothing` : custom suffixes, if unused plots files will be numbered starting from 1.
- `plot_geometry=true` : flag to generate geometry plot
- `plot_panels=false` : flag to include markers indicating panel edges in geometry plot
- `plot_pressure=false` : flag to generate surface pressures plot
- `plot_velocity=false` : flag to generate surface velocities plot
- `plot_boundary_layer=false` : flag to generate boundary layer plot
- `plot_streamlines=false` : flag to generate streamlines plot
- `verbose=false` : print verbose statements
- `kwargs...` : arguments passed into the plot functions (Plots keyword arguments/defaults to be used in every plot)
"""
function generate_plots(
    ::staticPlots,
    Plots,
    ins,
    outs;
    save_path="",
    static_file_type=".pdf",
    suffix=nothing,
    plot_geometry=true,
    plot_panels=false,
    plot_pressure=false,
    plot_velocity=false,
    plot_boundary_layer=false,
    plot_streamlines=false,
    verbose=false,
    kwargs...,
)
    if isnothing(suffix)
        suffix = ["_$(o)" for o in eachindex(outs)]
        lpad.(suffix, ndigits(length(outs)) + 1, "0")
    end

    # Extract useful items
    (; body_vortex_panels, wake_vortex_panels, rotor_source_panels) = ins.panels

    # - Geometry - #
    if plot_geometry
        verbose && println("Plotting Geometry")
        Plots.plot(
            plotGeometry(),
            body_vortex_panels,
            rotor_source_panels,
            wake_vortex_panels;
            plot_panels=plot_panels,
            kwargs...,
        )
        Plots.savefig(save_path * "geometry" * static_file_type)
    end

    # - Surface Pressure - #
    if plot_pressure
        verbose && println("Plotting Surface Pressure")
        for (i, out) in zip(suffix, outs)

            # underlay geometry first
            plt = Plots.plot(
                underlayGeometry(), body_vortex_panels, rotor_source_panels; kwargs...
            )

            # then Plots.plot pressure distributions
            Plots.plot!(
                plotCP(),
                body_vortex_panels,
                out.bodies,
                rotor_source_panels;
                subplot=2,
                inset=(1, Plots.bbox(0, 0, 1, 1)),
                kwargs...,
            )

            Plots.savefig(save_path * "surface_pressure" * i * static_file_type)
        end
    end

    # - Tangential Velocity - #
    if plot_velocity
        verbose && println("Plotting Surface Velocity")
        for (i, out) in zip(suffix, outs)
            # underlay geometry first
            plt = Plots.plot(
                underlayGeometry(), body_vortex_panels, rotor_source_panels; kwargs...
            )

            # then Plots.plot velocity distributions
            Plots.plot!(
                plotVtan(),
                body_vortex_panels,
                out.bodies,
                out.reference_values.Vref[],
                rotor_source_panels;
                subplot=2,
                inset=(1, Plots.bbox(0, 0, 1, 1)),
                kwargs...,
            )
            Plots.savefig(save_path * "surface_velocity" * i * static_file_type)
        end
    end

    # - Boundary Layer Stuff - #
    if plot_boundary_layer
        verbose && println("Plotting Boundary Layer")
        for (i, out) in zip(suffix, outs)
            # Plots.plot momentum thicknesses on top
            plt = Plots.plot(
                plotDuctGeometry(), body_vortex_panels; color=6, linewidth=0.5, kwargs...
            )

            Plots.plot!(
                plotMomentum(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                scale_thickness=5.0,
                kwargs...,
            )

            # Plots.plot stagnation point on top
            Plots.plot!(
                plotStagnation(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                markersize=4,
                color=3,
                markerstrokecolor=3,
                kwargs...,
            )
            Plots.savefig(save_path * "boundary_layer" * i * static_file_type)
        end
    end

    # - Streamlines - #
    if plot_streamlines
        verbose && println("Plotting Streamlines")
        for (i, out) in zip(suffix, outs)
            plt = Plots.plot(plotBodyGeometry(), body_vortex_panels; kwargs...)

            Plots.plot!(
                plotStreamlines(),
                body_vortex_panels,
                out.bodies.panel_strengths,
                wake_vortex_panels,
                out.wake.panel_strengths,
                rotor_source_panels,
                out.rotors.panel_strengths,
                out.reference_values.Vinf[];
                starting_radial_points=range(0.01, 0.7; length=50),
                axial_range=[-0.15, 0.5],
                step_limit=75,
                nominal_step_size=1e-2,
                integration_options=IntegrationOptions(),
                stag_tol=0.0,
                kwargs...,
            )

            Plots.plot!(
                plotStagnation(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                markersize=4,
                color=2,
                markerstrokecolor=2,
                kwargs...,
            )
            Plots.savefig(save_path * "streamlines" * i * static_file_type)
        end
    end

    return nothing
end

function generate_plots(
    ::animatedPlots,
    Plots,
    ins,
    outs;
    save_path="",
    static_file_type=".pdf",
    plot_geometry=true,
    plot_panels=false,
    plot_pressure=false,
    plot_velocity=false,
    plot_boundary_layer=false,
    plot_streamlines=false,
    verbose=false,
    kwargs...,
)

    # Extract useful items
    (; body_vortex_panels, wake_vortex_panels, rotor_source_panels) = ins.panels

    # - Geometry - #
    if plot_geometry
        verbose && println("Plotting Geometry")
        Plots.plot(
            plotGeometry(),
            body_vortex_panels,
            rotor_source_panels,
            wake_vortex_panels;
            plot_panels=plot_panels,
            kwargs...,
        )
        Plots.savefig(save_path * "geometry" * static_file_type)
    end

    # - Surface Pressure - #
    if plot_pressure
        verbose && println("Animating Surface Pressure")
        anim = Plots.Animation()
        for out in outs

            # underlay geometry first
            plt = Plots.plot(
                underlayGeometry(), body_vortex_panels, rotor_source_panels; kwargs...
            )

            # then Plots.plot pressure distributions
            Plots.plot!(
                plotCP(),
                body_vortex_panels,
                out.bodies,
                rotor_source_panels;
                subplot=2,
                inset=(1, Plots.bbox(0, 0, 1, 1)),
                (; kwargs..., background_color=:white)...,
            )

            Plots.frame(anim, plt)
        end
        Plots.gif(anim, save_path * "surface_pressure.gif")
    end

    # - Tangential Velocity - #
    if plot_velocity
        verbose && println("Animating Surface Velocity")
        anim = Plots.Animation()
        for out in outs
            # underlay geometry first
            plt = Plots.plot(
                underlayGeometry(), body_vortex_panels, rotor_source_panels; kwargs...
            )

            # then Plots.plot velocity distributions
            Plots.plot!(
                plotVtan(),
                body_vortex_panels,
                out.bodies,
                out.reference_values.Vref[],
                rotor_source_panels;
                subplot=2,
                inset=(1, Plots.bbox(0, 0, 1, 1)),
                (; kwargs..., background_color=:white)...,
            )
            Plots.frame(anim, plt)
        end
        Plots.gif(anim, save_path * "surface_velocity.gif")
    end

    # - Boundary Layer Stuff - #
    if plot_boundary_layer
        verbose && println("Animating Boundary Layer")
        anim = Plots.Animation()
        for out in outs
            # Plots.plot momentum thicknesses on top
            plt = Plots.plot(
                plotDuctGeometry(), body_vortex_panels; color=6, linewidth=0.5, kwargs...
            )

            Plots.plot!(
                plotMomentum(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                scale_thickness=5.0,
                kwargs...,
            )

            # Plots.plot stagnation point on top
            Plots.plot!(
                plotStagnation(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                color=3,
                markerstrokecolor=3,
                (; kwargs..., background_color=:white)...,
            )
            Plots.frame(anim, plt)
        end
        Plots.gif(anim, save_path * "boundary_layer.gif")
    end

    # - Streamlines - #
    if plot_streamlines
        verbose && println("Animating Streamlines")
        anim = Plots.Animation()
        for out in outs
            plt = Plots.plot(plotBodyGeometry(), body_vortex_panels; kwargs...)

            Plots.plot!(
                plotStreamlines(),
                body_vortex_panels,
                out.bodies.panel_strengths,
                wake_vortex_panels,
                out.wake.panel_strengths,
                rotor_source_panels,
                out.rotors.panel_strengths,
                out.reference_values.Vinf[];
                starting_radial_points=range(0.01, 0.7; length=50),
                axial_range=[-0.15, 0.5],
                step_limit=75,
                nominal_step_size=1e-2,
                integration_options=IntegrationOptions(),
                stag_tol=0.0,
                kwargs...,
            )

            Plots.plot!(
                plotStagnation(),
                out.bodies.boundary_layers,
                body_vortex_panels;
                markersize=2,
                color=2,
                markerstrokecolor=2,
                (; kwargs..., background_color=:white)...,
            )
            Plots.frame(anim, plt)
        end
        Plots.gif(anim, save_path * "streamlines.gif")
    end

    return nothing
end
