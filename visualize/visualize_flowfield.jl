using GeometricTools
const gt = GeometricTools
using DuctAPE
const dt = DuctAPE
using FLOWMath
const fm = FLOWMath

function visualize_flowfield(
    Vinf;
    body_panels=nothing,
    rotor_panels=nothing,
    wake_panels=nothing,
    mub=nothing,
    sigr=nothing,
    gamw=nothing,
    Gamr=nothing,
    num_blades=[1],
    rotorzloc=[0.0],
    Pmax=nothing,
    Pmin=nothing,
    verbose=false,
    run_name="velocity_field",
    save_path="",
    cellsizescale=0.005,
)

    ## -- SET UP -- ##
    # Bounds of grid
    if Pmax == nothing
        if body_panels == nothing
            Pmax = [
                1.1 * maximum(wake_panels.controlpoint[:, 1]),
                1.25 *
                maximum([rotor_panel.controlpoint[end, 2] for rotor_panel in rotor_panels]),
                0,
            ]
        elseif wake_panels == nothing
            Pmax = [
                1.25 * (maximum(body_panels.controlpoint[:, 1])),
                1.25 * (maximum(body_panels.controlpoint[:, 2])),
                0,
            ]
        else
            Pmax = [
                1.25 * (maximum(wake_panels.controlpoint[:, 1])),
                1.25 * (maximum(body_panels.controlpoint[:, 2])),
                0,
            ]
        end
    end
    if Pmin == nothing
        if body_panels == nothing
            Pmin = [-0.1 * maximum(wake_panels.controlpoint[:, 1]), eps(), 0]
        else
            Pmin = [-0.25 * (maximum(body_panels.controlpoint[:, 1])), eps(), 0]
        end
    end

    # Grid discretization
    dr = cellsizescale * (Pmax[2] - Pmin[2])                            # Cell size
    dx = dr

    NDIVS = ceil.(Int, (Pmax .- Pmin) ./ [dx, dr, 1]) # Divisions in each dimension

    # Generate grid
    @time fieldgrid = gt.Grid(Pmin, Pmax, NDIVS) # Grid

    if verbose
        println("\t"^(1) * "Grid size:\t\t$(NDIVS)")
    end
    if verbose
        println("\t"^(1) * "Number of nodes :\t$(fieldgrid.totnode)")
    end

    # Targets where to probe the velocity
    targets = fieldgrid.nodes'
    ntargets = size(targets, 1)

    ## -- ADD VELOCITIES -- ##

    # - Freestream velocities - #
    Uinf = repeat([Vinf 0.0 0.0], ntargets)
    Usinf = [[U[1], U[2], 0] for U in eachrow(Uinf)]
    gt.add_field(fieldgrid, "Usinf", "vector", Usinf, "node")

    Utot = copy(Uinf)

    if body_panels != nothing

        # - Body-induced velocities - #
        Ubody = dt.vfromdoubletpanels(targets[:, 1:2], body_panels.nodes, mub)

        # - Bodywake-induced velocities - #
        Ubodywake = dt.vfromTE(targets[:, 1:2], body_panels.TEnodes, mub)

        # body+bodywake
        Ubbw = Ubody .+ Ubodywake

        Utot[:, 1:2] .+= Ubbw

        Usbody = [[U[1], U[2], 0] for U in eachrow(Ubody)]
        Usbodywake = [[U[1], U[2], 0] for U in eachrow(Ubodywake)]
        Usbbw = [[U[1], U[2], 0] for U in eachrow(Ubbw)]
        gt.add_field(fieldgrid, "Ubody", "vector", Usbody, "node")
        gt.add_field(fieldgrid, "Ubodywake", "vector", Usbodywake, "node")
        gt.add_field(fieldgrid, "Ubodyandbodywake", "vector", Usbbw, "node")
    end

    if rotor_panels != nothing
        # - Rotor panel-induced velocities - #
        Urotor = similar(Ubody) .= 0.0
        for (ir, rp) in enumerate(rotor_panels)
            dt.vfromsourcepanels!(
                Urotor, targets[:, 1:2], rp.controlpoint, rp.len, sigr[ir]
            )
        end

        Utot[:, 1:2] .+= Urotor
        Usrotor = [[U[1], U[2], 0] for U in eachrow(Urotor)]
        gt.add_field(fieldgrid, "Urotor", "vector", Usrotor, "node")
    end

    if Gamr != nothing && wake_panels != nothing
        # - Rotor induced tangential velocity - #
        # TODO: only within the wake, only defined on wake lines.
        # TODO: need to get the x spacing of the targets,
        # then need to get the wake values at those points,
        # then need to interpolate the wake values to the target radial locations

        # Initialize
        Utheta = zeros(size(targets, 1))

        # reshape the wake panel control points into the wake sheets
        nsheets = size(Gamr, 1) + 1
        nwakex = Int(wake_panels.totpanel / nsheets)
        wakecpx = reshape(wake_panels.controlpoint[:, 1], (nwakex, nsheets))'
        wakecpr = reshape(wake_panels.controlpoint[:, 2], (nwakex, nsheets))'

        # get the target xpositions between the wake start and end points.
        tx1 = findfirst(x -> x >= wakecpx[1, 1], targets[:, 1])
        tx2 = findfirst(x -> x >= wakecpx[1, end], targets[:, 1]) - 1
        tgridx = repeat(targets[tx1:tx2, 1]; inner=(1, nsheets))'

        # spline the wake sheets
        tgridr = similar(tgridx) .= 0.0
        for i in 1:nsheets
            tgridr[i, :] = FLOWMath.linear(wakecpx[i, :], wakecpr[i, :], tgridx[i, :])
        end

        # Third, get the vtheta values at the interpolated points
        vthetagrid = similar(tgridx) .= 0.0

        Gamma_tilde_be = dt.calculate_net_circulation(Gamr, num_blades)
        Gamma_tilde_ws = [
            Gamma_tilde_be[1]
            (Gamma_tilde_be[2:end] .+ Gamma_tilde_be[1:(end - 1)]) / 2
            Gamma_tilde_be[end]
        ]
        for ix in 1:size(tgridx, 2)
            gtid = findlast(x -> x <= tgridx[1, ix], rotorzloc)
            vthetagrid[:, ix] = dt.calculate_vtheta(Gamma_tilde_ws[:, gtid], tgridr[:, ix])
        end

        # Fourth, loop through targets, if they are within the wake, interpolate vtheta values.
        for (it, tar) in enumerate(eachrow(targets))
            if tgridx[1, 1] <= tar[1] <= tgridx[1, end]
                # println("target: ", tar, "inside xrange")
                #we're inside the wake x range
                # find which xstation we're at
                xid = searchsortedfirst(tgridx[1, :], tar[1])
                if tgridr[1, xid] <= tar[2] <= tgridr[end, xid]
                    # println("target: ", tar, "inside rrange")
                    # we're inside the wake radially
                    # get the tangential velocity
                    Utheta[it] = FLOWMath.akima(tgridr[:, xid], vthetagrid[:, xid], tar[2])
                end
            end
        end

        # put Utheta into Utot
        Utot[:, 3] .+= Utheta

        Ustheta = [[0, 0, U[1]] for U in eachrow(Utheta)]
        gt.add_field(fieldgrid, "Utheta", "vector", Ustheta, "node")
    end

    if wake_panels != nothing
        # - Wake-induced velocities - #
        Uwake = dt.vfromvortexpanels(
            targets[:, 1:2], wake_panels.controlpoint, wake_panels.len, gamw
        )

        Utot[:, 1:2] .+= Uwake
        Uswake = [[U[1], U[2], 0] for U in eachrow(Uwake)]
        gt.add_field(fieldgrid, "Uwake", "vector", Uswake, "node")
    end

    # Put the Velocities together in the right formats
    Ustot = [[U[1], U[2], U[3]] for U in eachrow(Utot)]

    # - Save VTKs - #
    # Save fields
    gt.add_field(fieldgrid, "Utot", "vector", Ustot, "node")

    # Output fluid domain
    @time vtks = gt.save(fieldgrid, run_name; path=save_path, format="vtk")

    return nothing
end

function visualize_surfaces(;
    body_panels=nothing,
    rotor_panels=nothing,
    wake_panels=nothing,
    verbose=false,
    run_name="",
    save_path="",
)
    if body_panels != nothing

        # - Body Geometry - #
        cps = [[p[1], p[2], 0] for p in eachrow(body_panels.controlpoint)]

        point_data = [
            Dict(
                "field_name" => "normal",
                "field_type" => "vector",
                "field_data" => [[n[1], n[2], 0] for n in eachrow(body_panels.normal)],
            ),
            Dict(
                "field_name" => "tangent",
                "field_type" => "vector",
                "field_data" => [[t[1], t[2], 0] for t in eachrow(body_panels.tangent)],
            ),
        ]

        gt.generateVTK(
            run_name * "body_controlpoints", cps; point_data=point_data, path=save_path
        )
    end

    if wake_panels != nothing

        # - Geometry - #
        cps = [[p[1], p[2], 0] for p in eachrow(wake_panels.controlpoint)]

        point_data = [
            Dict(
                "field_name" => "normal",
                "field_type" => "vector",
                "field_data" => [[n[1], n[2], 0] for n in eachrow(wake_panels.normal)],
            ),
            Dict(
                "field_name" => "tangent",
                "field_type" => "vector",
                "field_data" => [[t[1], t[2], 0] for t in eachrow(wake_panels.tangent)],
            ),
        ]

        gt.generateVTK(
            run_name * "wake_controlpoints", cps; point_data=point_data, path=save_path
        )
    end
    if rotor_panels != nothing
        for (ir, rp) in enumerate(rotor_panels)
            # - Geometry - #
            cps = [[p[1], p[2], 0] for p in eachrow(rp.controlpoint)]

            point_data = [
                Dict(
                    "field_name" => "normal",
                    "field_type" => "vector",
                    "field_data" => [[n[1], n[2], 0] for n in eachrow(rp.normal)],
                ),
                Dict(
                    "field_name" => "tangent",
                    "field_type" => "vector",
                    "field_data" => [[t[1], t[2], 0] for t in eachrow(rp.tangent)],
                ),
            ]

            gt.generateVTK(
                run_name * "rotor$(ir)_controlpoints",
                cps;
                point_data=point_data,
                path=save_path,
            )
        end
    end

    return nothing
end
