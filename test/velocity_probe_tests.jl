@testset "Velocity Probe" begin

    # define probe locations
    probe_poses = [0.25 1.5; 0.75 1.0]

    # define required body panel geometry
    duct_coordinates = [1.0 2.0; 0.0 2.0]
    body_doublet_panels = dt.generate_panels([duct_coordinates])
    # duct_coordinates = [1.0 2.0; 0.5 1.75; 0.0 2.0]
    # hub_coordinates = [0.0 0.0; 0.5 0.25; 1.0 0.125]
    # body_doublet_panels = dt.generate_panels([duct_coordinates, hub_coordinates])

    # define required rotor geometry
    nrotor = 2
    rotor_parameters = [(; xrotor=0.0); (; xrotor=0.5)]
    rrotor = [0.5; 1.5]
    rotor_source_panels = [
        dt.generate_rotor_panels(rotor_parameters[i].xrotor, rrotor) for i in 1:nrotor
    ]

    # define required blade element stuff
    blade_elements = [(; B=1); (; B=1)]

    # define required wake geometry
    xwake = [0.0 0.0; 0.5 0.5; 1.0 1.0]
    rwake = [0.5 1.5; 0.5 1.5; 0.5 1.5]
    wake_vortex_panels = dt.generate_wake_panels(xwake, rwake)

    # # sanity plot for manual check
    # plot()
    # plot!(probe_poses[:,1], probe_poses[:,2], seriestype=:scatter,label="")
    # plot!(duct_coordinates[:,1], duct_coordinates[:,2], seriestype=:scatter, label="",color=:red)
    # plot!(hub_coordinates[:,1], hub_coordinates[:,2], seriestype=:scatter, label="",color=:red)
    # plot!(xwake,rwake, seriestype=:scatter, label="",color=:gray)
    # plot!(xrotor*ones(2),rrotor, seriestype=:scatter, label="",color=:black)

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
    states = ones(12)

    # run velocity probing function
    vx, vr, vt = dt.probe_velocity_field(probe_poses, inputs, states)

    # Get manual x and r values
    # body induced velocities
    vb = zeros(2, 2)
    dt.vfromdoubletpanels!(vb, probe_poses, body_doublet_panels.nodes, [1.0])

    # rotor induced velocities
    vr1 = zeros(2, 2)
    dt.vfromsourcepanels!(
        vr1,
        probe_poses,
        rotor_source_panels[1].controlpoint,
        rotor_source_panels[1].len,
        [1.0],
    )

    vr2 = zeros(2, 2)
    dt.vfromsourcepanels!(
        vr2,
        probe_poses,
        rotor_source_panels[2].controlpoint,
        rotor_source_panels[2].len,
        [1.0],
    )

    # wake induced velocities
    vw = zeros(2, 2)
    dt.vfromvortexpanels!(
        vw, probe_poses, wake_vortex_panels.controlpoint, wake_vortex_panels.len, ones(4)
    )

    vxmanual = vb[:, 1] .+ vr1[:, 1] .+ vr2[:, 1] .+ vw[:, 1]
    vrmanual = vb[:, 2] .+ vr1[:, 2] .+ vr2[:, 2] .+ vw[:, 2]

    @test isapprox(vx, vxmanual)
    @test isapprox(vr, vrmanual)

    # get manual tangential values

    vtmanual = zeros(2,2)



    @test vt == vtmanual
end
