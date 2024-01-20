"""
"""
function reformat_controlpoint_velocities(
    vels, pidx, nhub_interface_panels, nduct_interface_panels
)

    # - hub is element 1 - #
    hub_vel = vels[pidx[1, 1]:pidx[1, 2], :]
    # - duct is element 2 - #
    duct_vel = vels[pidx[2, 1]:pidx[2, 2], :]
    # - rotor is element 3 - #
    rotor_vel = vels[pidx[3, 1]:pidx[3, 2], :]
    # - hubwake is element 4 - #
    hubwake_vel = vels[pidx[4, 1]:pidx[4, 2], :]
    # - interior wakes are element 5-13 - #
    interior_wake_vel = vels[pidx[5, 1]:pidx[13, 2], :]
    # - ductwake is element 14 - #
    ductwake_vel = vels[pidx[14, 1]:pidx[14, 2], :]

    # - assemble full wake velocities - #
    full_wake_vel = [
        reverse(hub_vel[1:nhub_interface_panels, :]; dims=1)
        hubwake_vel
        interior_wake_vel
        duct_vel[(end - nduct_interface_panels + 1):end, :]
        ductwake_vel[1:end, :]
    ]

    return (;
        duct_vel,
        hub_vel,
        ductwake_vel,
        hubwake_vel,
        interior_wake_vel,
        full_wake_vel,
        rotor_vel,
    )
end

"""
"""
function reformat_circulation(BGAM, B)
    return BGAM ./ B
end

"""
"""
function reformat_extended_blade_elements(vals, B)

    # VREL(1,IR,NR),VTBG,WTB,WWB,PHIB,XI,ALF,BETAR(IR,NR),REY,SECSIG,SECTAGR,CLR(IR,NR),CDR(IR,NR),BGX(IR)
    return (
        Wz_rotor=vals[:, 1],
        vtheta_self=vals[:, 2],
        Wtheta_rotor=vals[:, 3],
        Wmag_rotor=vals[:, 4],
        phi=vals[:, 5],
        xi=vals[:, 6],
        alpha=vals[:, 7],
        twist=vals[:, 8],
        reynolds=vals[:, 9],
        solidity=vals[:, 10],
        stagger=vals[:, 11],
        cl=vals[:, 12],
        cd=vals[:, 13],
        Gamr_est=vals[:, end] ./ B,
    )
end

"""
"""
function reformat_gamw(GAMTH, nidx, nhub_interface_nodes, nduct_interface_nodes)

    # Format the Wake Strengths
    hub_gamw = GAMTH[nidx[1, 1]:nidx[1, 2]]
    duct_gamw = GAMTH[nidx[2, 1]:nidx[2, 2]]
    hubwake_gamw = GAMTH[nidx[4, 1]:nidx[4, 2]]
    interior_gamw = GAMTH[nidx[5, 1]:nidx[13, 2]]
    ductwake_gamw = GAMTH[nidx[14, 1]:nidx[14, 2]]

    return [
        reverse(hub_gamw[1:nhub_interface_nodes]) #don't repeat the TE node
        hubwake_gamw
        interior_gamw
        duct_gamw[(end - nduct_interface_nodes + 1):end]
        ductwake_gamw[1:end]#don't repeat the TE node
    ]
end

"""
"""
function reformat_wake_panel_vels(VMAV, nhub_interface_panels, nduct_interface_panels)
    # fill in first and last wake sheet with zeros?
    # these are the average velocities at the wake panels
    return [
        zeros(nhub_interface_panels) #start with zeros on hub wake interface, don't repeat TE node
        VMAV[1:(end - pidranges[end])] #include hub wake and interior wakes
        zeros(nduct_interface_panels) # start with zeros in duct wake interface
        VMAV[(end - pidranges[end] + 1):end] # finish with duct wake, don't repeate TE node
    ]
end

"""
"""
function reformat_wmavg(VMAVG, nhub_interface_nodes, nduct_interface_nodes)
    # fill in first and last wake sheet with zeros?
    # these are the average velocities at the wake nodes
    return [
        zeros(nhub_interface_nodes) #start with zeros on hub wake interface, don't repeat TE node
        VMAVG[1:(end - nidranges[end] - 1)] #include hub wake and interior wakes
        zeros(nduct_interface_nodes) # start with zeros in duct wake interface
        VMAVG[(end - nidranges[end]):end] # finish with duct wake, don't repeate TE node
    ]
end
