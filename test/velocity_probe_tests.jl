@testset "Velocity Probe" begin

    # define probe locations
    probe_poses = [
        -1.0 1.0 # out in front
        0.25 0.01 # out below
        0.375 2.5 # out above
        1.75 1.0 # out behind
        0.0 1.0 # on rotor mid
        0.0 1.5 # on rotor top edge
        0.0 0.5 # on rotor bottom edge
        0.25 1.5 # on wake top edge on control point
        0.375 1.5 # on wake top edge not on control point
        0.375 0.5 # on wake bottom edge
        0.75 0.75 # on wake back edge
        0.25 1.0 # aligned with first wake cp
        0.625 0.75 # arbitrary postition in wake, close to rotor
        1.25 0.75 # arbitrary postition in wake, not on anything
    ]

    np = length(probe_poses[:, 1])

    # define required body panel geometry
    duct_coordinates = [1.0 2.0; 0.0 2.0]
    body_doublet_panels = dt.generate_panels([duct_coordinates])
    # duct_coordinates = [1.0 2.0; 0.5 1.75; 0.0 2.0]
    # hub_coordinates = [0.0 0.0; 0.5 0.25; 1.0 0.125]
    # body_doublet_panels = dt.generate_panels([duct_coordinates, hub_coordinates])

    # define required rotor geometry
    nrotor = 2
    rotor_parameters = [(; xrotor=0.0); (; xrotor=0.5)]
    rrotor = [0.5; 1.0; 1.5]
    rotor_source_panels = [
        dt.generate_rotor_panels(rotor_parameters[i].xrotor, rrotor) for i in 1:nrotor
    ]

    # define required blade element stuff
    blade_elements = [(; B=1, xrotor=0.0); (; B=1, xrotor=0.5)]

    # define required wake geometry
    xwake = [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0; 1.5 1.5 1.5]
    rwake = [0.5 1.0 1.5; 0.5 1.0 1.5; 0.5 1.0 1.5; 0.5 1.0 1.5]
    wake_vortex_panels = dt.generate_wake_panels(xwake, rwake)

    # sanity plot for manual check
    plot()
    plot!(probe_poses[:, 1], probe_poses[:, 2]; seriestype=:scatter, label="")
    plot!(
        duct_coordinates[:, 1],
        duct_coordinates[:, 2];
        seriestype=:scatter,
        label="",
        color=:red,
    )
    # plot!(hub_coordinates[:,1], hub_coordinates[:,2], seriestype=:scatter, label="",color=:red)
    plot!(xwake, rwake; seriestype=:scatter, label="", color=:gray)
    plot!(
        [0.0; 0.5; 0.0; 0.5],
        [0.5; 0.5; 1.5; 1.5];
        seriestype=:scatter,
        label="",
        color=:black,
    )

    # put geometry into inputs along with Vinf
    inputs = (;
        body_doublet_panels,
        num_body_panels=body_doublet_panels.npanels,
        rotor_source_panels,
        rotor_panel_centers=reduce(
            hcat, [rotor_source_panels[i].controlpoint[:, 2] for i in 1:nrotor]
        ),
        wake_vortex_panels,
        num_wake_x_panels=2,
        blade_elements,
    )

    # prescribe strengths for the panels (body, wake, circulation, rotor)
    states = ones(15)
    mub, gamw, Gamr, sigr = dt.extract_state_variables(states, inputs)
    Gamr[2, :] .= 2.0

    # run velocity probing function
    vx, vr, vt = dt.probe_velocity_field(probe_poses, inputs, states)

    # Get manual x and r values
    # body induced velocities
    vb = zeros(np, 2)
    dt.vfromdoubletpanels!(vb, probe_poses, body_doublet_panels.nodes, mub)

    # rotor induced velocities
    vr1 = zeros(np, 2)
    dt.vfromsourcepanels!(
        vr1,
        probe_poses,
        rotor_source_panels[1].controlpoint,
        rotor_source_panels[1].len,
        sigr[:, 1],
    )

    vr2 = zeros(np, 2)
    dt.vfromsourcepanels!(
        vr2,
        probe_poses,
        rotor_source_panels[2].controlpoint,
        rotor_source_panels[2].len,
        sigr[:, 2],
    )

    # wake induced velocities
    vw = zeros(np, 2)
    dt.vfromvortexpanels!(
        vw, probe_poses, wake_vortex_panels.controlpoint, wake_vortex_panels.len, gamw
    )

    vxmanual = vb[:, 1] .+ vr1[:, 1] .+ vr2[:, 1] .+ vw[:, 1]
    vrmanual = vb[:, 2] .+ vr1[:, 2] .+ vr2[:, 2] .+ vw[:, 2]

    @test isapprox(vx, vxmanual)
    @test isapprox(vr, vrmanual)

    # get manual tangential values

    vtmanual = [
        0.0 # out in front
        0.0 # out below
        0.0 # out above
        0.0 # out behind
        0.75 / (2 * pi) # on rotor mid
        1.0 / (2 * pi * 1.5) # on rotor top edge
        0.5 / (2 * pi * 0.5) # on rotor bottom edge
        2.0 / (2 * pi * 1.5) # on wake top edge on control point
        2.0 / (2 * pi * 1.5) # on wake top edge not on control point
        1.0 / (2 * pi * 0.5) # on wake bottom edge
        2.5 / (2 * pi * probe_poses[11, 2]) # on wake back edge
        1.5 / (2 * pi * probe_poses[12, 2]) # aligned with first wake cp
        2.34375 / (2 * pi * probe_poses[13, 2]) # arbitrary postition in wake, close to rotor
        2.5 / (2 * pi * probe_poses[13, 2]) # arbitrary postition in wake, not on anything
    ]

    @test vt == vtmanual
end
