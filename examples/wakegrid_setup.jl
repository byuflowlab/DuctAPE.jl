"""
"""
function setup_wakegrid(ductgeometry, ductsplines, rotors; plotgrid=false)
    #set up grid options
    num_radial_stations = length(rnondim1)
    grid_options = DuctTAPE.defineGridOptions(num_radial_stations)

    #generate starting grid points
    # xg, rg, nx, nr = DuctTAPE.generate_grid_points(
    #     ductgeometry, ductsplines, rotors, grid_options
    # )

    # relax grid points
    # grid = relax_grid(xg, rg, nx, nr)

    # Generate initialized grid points (generates and relaxes grid points in one function)
    grid = DuctTAPE.initialize_grid(ductgeometry, ductsplines, rotors, grid_options)

    if plotgrid
        #extract grid points
        xg = grid.x_grid_points
        rg = grid.r_grid_points
        nx = grid.nx
        nr = grid.nr

        # PLOT GEOMETRY
        plot(; aspectratio=:equal, xlabel="x", ylabel="r", legend=true, label="")

        plot!(
            ductgeometry.wallinnerxcoordinates,
            ductgeometry.wallinnerrcoordinates;
            color=1,
            linewidth=2,
            label="inner duct wall",
        )
        plot!(
            ductgeometry.wallouterxcoordinates,
            ductgeometry.wallouterrcoordinates;
            color=1,
            linestyle=:dot,
            linewidth=2,
            label="outer duct wall",
        )

        plot!(
            ductgeometry.hubxcoordinates,
            ductgeometry.hubrcoordinates;
            color=2,
            linewidth=2,
            label="hub",
        )

        blade1 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotors[1])
        blade2 = DuctTAPE.initialize_blade_dimensions(ductgeometry, ductsplines, rotors[2])

        plot!(
            rotors[1].xlocation .* ductgeometry.chord * ones(nr),
            blade1.rdim;
            color=4,
            linewidth=2,
            label="rotors",
        )
        plot!(
            rotors[2].xlocation .* ductgeometry.chord * ones(nr),
            blade2.rdim;
            color=4,
            linewidth=2,
            label="",
        )

        plot!(xg, rg; color=3, label="", linewidth=0.5)
        plot!(xg', rg'; color=3, label="", linewidth=0.5)
        plot!(
            xg[:, 1] .* NaN, rg[:, 1] .* NaN; color=3, label="Initial Grid", linewidth=0.5
        )

        savefig("./examples/test_grid.pdf")
    end

    return grid
end

