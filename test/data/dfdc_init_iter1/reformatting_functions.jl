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
        reverse(hub_gamw[2:(nhub_interface_nodes + 1)])
        hubwake_gamw
        interior_gamw
        duct_gamw[(end - nduct_interface_nodes):(end - 1)]#don't repeat the TE node
        ductwake_gamw[1:end]
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

"""
"""
function reformat_body_cp(cp, nidx)
    duct_cp = cp[nidx[2, 1]:nidx[2, 2]]
    hub_cp = cp[nidx[1, 1]:nidx[1, 2]]

    return duct_cp, hub_cp
end

"""
"""
function reformat_rhs(aic, vec, pidx)
    rhs = zeros(length(vec))
    for k in 1:size(aic, 1)
        for j in 1:size(aic, 2)
            rhs[k] += aic[k, j] * vec[j]
        end
    end

    # rhs = aic * vec

    hub = reverse(rhs[pidx[1, 1]:(pidx[1, 2])])
    duct = reverse(rhs[(pidx[2, 1] + 2):(pidx[2, 2] + 2)])

    return [
        duct
        hub
        rhs[pidx[2, 2] + 4] #duct interior controlpoint
        rhs[pidx[2, 2] + 3] #kutta
        rhs[(pidx[2, 1])] #hub prescribed LE
        rhs[pidx[2, 1] + 1] #hub interior controlpoint
    ]
end

function reformat_rhs(rhs, pidx)
    hub = reverse(rhs[pidx[1, 1]:(pidx[1, 2])])
    duct = reverse(rhs[(pidx[2, 1] + 2):(pidx[2, 2] + 2)])

    return [
        duct
        hub
        rhs[pidx[2, 2] + 4] #duct interior controlpoint
        rhs[pidx[2, 2] + 3] #kutta
        rhs[(pidx[2, 1])] #hub prescribed LE
        rhs[pidx[2, 1] + 1] #hub interior controlpoint
    ]
end

function reformat_sol(sol, nidx)
    hub = reverse(sol[nidx[1, 1]:(nidx[1, 2])])
    duct = reverse(sol[(nidx[2, 1]):(nidx[2, 2])])

    return [
        duct
        hub
        sol[nidx[2, 2] + 1] #duct interior controlpoint
        sol[nidx[2, 2] + 2] #hub interior controlpoint
    ]
end

"""
"""
function reformat_wake_aic(
    aic, nidx, pidx, nidranges, pidranges, nhub_interface_nodes, nduct_interface_nodes
)

    # define column ranges
    # NOTE: there are 2 extra nodes. It seems that the TE node is calculated both from the TE body panel as well as the first wake panel.
    # they should be able to be simply added together.
    hwi_cols = nidx[1, 1]:(nhub_interface_nodes + 1)
    restofhub = (nhub_interface_nodes + 2):nidx[1, 2]
    restofduct = nidx[2, 1]:(nidx[2, 2] - nduct_interface_nodes - 1)
    dwi_cols = (nidx[2, 2] - nduct_interface_nodes):nidx[2, 2]
    rotor = nidx[3, 1]:nidx[3, 2]
    hw_cols = nidx[4, 1]:nidx[4, 2]
    iw_cols = nidx[5, 1]:nidx[end - 1, 2]
    dw_cols = nidx[end, 1]:nidx[end, 2]

    # define row ranges
    h_rows = reverse(pidx[1, 1]:pidx[1, 2]) # hub panels
    h_pn = pidx[2, 1]:pidx[2, 1] #prescribed node
    h_itcp = (pidx[2, 1] + 1):(pidx[2, 1] + 1) # interior ctrpt
    d_rows = reverse((pidx[2, 1] + 2):(pidx[2, 2] + 2)) # duct panels
    kutta = (pidx[2, 2] + 3):(pidx[2, 2] + 3) #kutta
    d_itcp = (pidx[2, 2] + 4):(pidx[2, 2] + 4) # interior ctrpt

    # put things in DuctAPE order
    # initialize matrix first, then fill in
    mat = zeros(167, 561)

    # define dt row and column ranges
    duct_r = 1:pidranges[2]
    hub_r = (pidranges[2] + 1):(pidranges[2] + pidranges[1])

    hubinterface_c = 1:(nhub_interface_nodes + 1)
    hubwake_c = (nhub_interface_nodes + 1):(nhub_interface_nodes + nidranges[4])
    interior_c = (hubwake_c[end] + 1):(hubwake_c[end] + nidranges[5] * 9)
    ductinterface_c = (interior_c[end] + 1):(interior_c[end] + nduct_interface_nodes + 1)
    ductwake_c = ductinterface_c[end]:(ductinterface_c[end] + nidranges[end] - 1)

    #hub wake interface on duct
    mat[duct_r, hubinterface_c] = reverse(aic[d_rows, hwi_cols]; dims=2)
    #hub wake on duct
    mat[duct_r, hubwake_c] += aic[d_rows, hw_cols]
    #interior wakes on duct
    mat[duct_r, interior_c] = aic[d_rows, iw_cols]
    # duct wake interface on duct
    mat[duct_r, ductinterface_c] = aic[d_rows, dwi_cols]
    # duct wake on duct
    mat[duct_r, ductwake_c] = aic[d_rows, dw_cols]

    #hub wake interface on hub
    mat[hub_r, hubinterface_c] = reverse(aic[h_rows, hwi_cols]; dims=2)
    #hub wake on hub
    mat[hub_r, hubwake_c] += aic[h_rows, hw_cols]
    #interior wakes on hub
    mat[hub_r, interior_c] = aic[h_rows, iw_cols]
    # duct wake interface on hub
    mat[hub_r, ductinterface_c] = aic[h_rows, dwi_cols]
    # duct wake on hub
    mat[hub_r, ductwake_c] = aic[h_rows, dw_cols]

    # interior point influence
    pcp_mat = zeros(2, 561)

    #hub wake interface on duct interior point
    pcp_mat[1:1, hubinterface_c] = reverse(aic[d_itcp, hwi_cols]; dims=2)
    #hub wake on duct interior point
    pcp_mat[1:1, hubwake_c] += aic[d_itcp, hw_cols]
    #interior wakes on duct interior point
    pcp_mat[1:1, interior_c] = aic[d_itcp, iw_cols]
    # duct wake interface on duct interior point
    pcp_mat[1:1, ductinterface_c] = aic[d_itcp, dwi_cols]
    # duct wake on duct interior point
    pcp_mat[1:1, ductwake_c] = aic[d_itcp, dw_cols]

    #hub wake interface on hub interior point
    pcp_mat[2:2, hubinterface_c] = reverse(aic[h_itcp, hwi_cols]; dims=2)
    #hub wake on hub interior point
    pcp_mat[2:2, hubwake_c] += aic[h_itcp, hw_cols]
    #interior wakes on hub interior point
    pcp_mat[2:2, interior_c] = aic[h_itcp, iw_cols]
    # duct wake interface on hub interior point
    pcp_mat[2:2, ductinterface_c] = aic[h_itcp, dwi_cols]
    # duct wake on hub interior point
    pcp_mat[2:2, ductwake_c] = aic[h_itcp, dw_cols]

    return mat, pcp_mat
end

"""
"""
function reformat_rotor_aic(aic, nidx, pidx, nidranges, pidranges)

    # define column ranges
    rotor = nidx[3, 1]:nidx[3, 2]

    # define row ranges
    h_rows = reverse(pidx[1, 1]:pidx[1, 2]) # hub panels
    h_pn = pidx[2, 1]:pidx[2, 1] #prescribed node
    h_itcp = (pidx[2, 1] + 1):(pidx[2, 1] + 1) # interior ctrpt
    d_rows = reverse((pidx[2, 1] + 2):(pidx[2, 2] + 2)) # duct panels
    kutta = (pidx[2, 2] + 3):(pidx[2, 2] + 3) #kutta
    d_itcp = (pidx[2, 2] + 4):(pidx[2, 2] + 4) # interior ctrpt

    mat = zeros(167, 11)

    # define dt row and column ranges
    duct_r = 1:pidranges[2]
    hub_r = (pidranges[2] + 1):(pidranges[2] + pidranges[1])

    #hub wake interface on duct
    mat[duct_r, :] = aic[d_rows, rotor]
    mat[hub_r, :] = aic[h_rows, rotor]

    # interior point influence
    pcp_mat = zeros(2, 11)

    #hub wake interface on duct interior point
    pcp_mat[1:1, :] = aic[d_itcp, rotor]
    pcp_mat[2:2, :] = aic[h_itcp, rotor]

    return mat, pcp_mat
end

"""
"""
function reformat_body_vhat(
    vz, vr, nidx, pidx, pidranges, nhub_interface_panels, nduct_interface_panels
)
    # body to rotor
    b2rvhat = zeros(pidranges[3], nidx[2, 2], 2)
    rotorrows = pidx[3, 1]:pidx[3, 2]
    bodycols = nidx[1, 1]:nidx[2, 2]
    bodyrows = pidx[1, 1]:pidx[2, 2]

    b2rvhat[:, :, 1] = vz[rotorrows, bodycols]
    b2rvhat[:, :, 2] = vr[rotorrows, bodycols]

    # body to wake
    b2wvhat = zeros(
        sum(pidranges[4:end]) + nhub_interface_panels + nduct_interface_panels,
        nidx[2, 2],
        2,
    )

    hwirows = reverse(1:nhub_interface_panels)
    hwrows = pidx[4, 1]:pidx[4, 2]
    wrows = pidx[5, 1]:pidx[13, 2]
    dwirows = (pidx[2, 2] - nduct_interface_panels + 1):pidx[2, 2]
    dwrows = pidx[end, 1]:pidx[end, 2]

    b2wvhat[:, :, 1] = [
        vz[hwirows, bodycols]
        vz[hwrows, bodycols]
        vz[wrows, bodycols]
        vz[dwirows, bodycols]
        vz[dwrows, bodycols]
    ]
    b2wvhat[:, :, 2] = [
        vr[hwirows, bodycols]
        vr[hwrows, bodycols]
        vr[wrows, bodycols]
        vr[dwirows, bodycols]
        vr[dwrows, bodycols]
    ]

    # body to body
    b2bvhat = zeros(sum(pidranges[1:2]), nidx[2, 2], 2)
    b2bvhat[:, :, 1] = vz[bodyrows, bodycols]
    b2bvhat[:, :, 2] = vr[bodyrows, bodycols]

    return reverse(reverse(b2bvhat; dims=2); dims=1),
    reverse(b2rvhat; dims=2),
    reverse(b2wvhat; dims=2)
end

"""
"""
function reformat_rotor_vhat(
    vz, vr, nidx, nidranges, pidx, pidranges, nhub_interface_panels, nduct_interface_panels
)
    # rotor body
    r2bvhat = zeros(pidx[2, 2], nidranges[3], 2) #panels, nodes, z/r
    bodyrows = pidx[1, 1]:pidx[2, 2]
    rotorcols = nidx[3, 1]:nidx[3, 2]
    rotorrows = pidx[3, 1]:pidx[3, 2]

    r2bvhat[:, :, 1] = vz[bodyrows, rotorcols]
    r2bvhat[:, :, 2] = vr[bodyrows, rotorcols]

    # rotor to wake
    r2wvhat = zeros(
        sum(pidranges[4:end]) + nhub_interface_panels + nduct_interface_panels,
        nidranges[3],
        2,
    )

    hwirows = reverse(1:nhub_interface_panels)
    hwrows = pidx[4, 1]:pidx[4, 2]
    wrows = pidx[5, 1]:pidx[13, 2]
    dwirows = (pidx[2, 2] - nduct_interface_panels + 1):pidx[2, 2]
    dwrows = pidx[end, 1]:pidx[end, 2]

    r2wvhat[:, :, 1] = [
        vz[hwirows, rotorcols]
        vz[hwrows, rotorcols]
        vz[wrows, rotorcols]
        vz[dwirows, rotorcols]
        vz[dwrows, rotorcols]
    ]
    r2wvhat[:, :, 2] = [
        vr[hwirows, rotorcols]
        vr[hwrows, rotorcols]
        vr[wrows, rotorcols]
        vr[dwirows, rotorcols]
        vr[dwrows, rotorcols]
    ]

    # rotor to rotor
    r2rvhat = zeros(pidranges[3], nidranges[3], 2)
    r2rvhat[:, :, 1] = vz[rotorrows, rotorcols]
    r2rvhat[:, :, 2] = vr[rotorrows, rotorcols]

    return reverse(r2bvhat; dims=1), r2rvhat, r2wvhat
end

"""
"""
function reformat_wake_vhat(
    vz, vr, nidx, nidranges, pidx, pidranges, nhub_interface_panels, nduct_interface_panels
)

    # - define row (panels) ranges - #
    # body
    bodyrows = pidx[1, 1]:pidx[2, 2]
    #rotor
    rotorrows = pidx[3, 1]:pidx[3, 2]
    #wake
    hwirows = reverse(1:nhub_interface_panels)
    hwrows = pidx[4, 1]:pidx[4, 2]
    wrows = pidx[5, 1]:pidx[13, 2]
    dwirows = (pidx[2, 2] - nduct_interface_panels + 1):pidx[2, 2]
    dwrows = pidx[end, 1]:pidx[end, 2]

    # - Define wake column (node) ranges - #
    hwicols = reverse(nidx[1, 1]:(nhub_interface_nodes))
    restofhub = (nhub_interface_nodes + 1):nidx[1, 2]
    restofduct = nidx[2, 1]:(nidx[2, 2] - nduct_interface_nodes)
    dwicols = (nidx[2, 2] - nduct_interface_nodes + 1):(nidx[2, 2])
    rotor = nidx[3, 1]:nidx[3, 2]
    hwcols = nidx[4, 1]:nidx[4, 2]
    wcols = nidx[5, 1]:nidx[end - 1, 2]
    dwcols = nidx[end, 1]:nidx[end, 2]
    nwn = sum(
        [length(hwicols); length(hwcols); length(wcols); length(dwicols); length(dwcols)]
    )

    # - wake to body - #
    w2bvhat = zeros(pidx[2, 2], nwn, 2) #panels, nodes, z/r

    w2bvhat[:, :, 1] = [
        vz[bodyrows, hwicols] vz[bodyrows, hwcols] vz[bodyrows, wcols] vz[bodyrows, dwicols] vz[bodyrows, dwcols]
    ]
    w2bvhat[:, :, 2] = [
        vr[bodyrows, hwicols] vr[bodyrows, hwcols] vr[bodyrows, wcols] vr[bodyrows, dwicols] vr[bodyrows, dwcols]
    ]

    # wake to wake
    w2wvhat = zeros(
        sum(pidranges[4:end]) + nhub_interface_panels + nduct_interface_panels, nwn, 2
    )

    hwirows = reverse(1:nhub_interface_panels)
    hwrows = pidx[4, 1]:pidx[4, 2]
    wrows = pidx[5, 1]:pidx[13, 2]
    dwirows = (pidx[2, 2] - nduct_interface_panels + 1):pidx[2, 2]
    dwrows = pidx[end, 1]:pidx[end, 2]

    w2wvhat[:, :, 1] = [
        vz[hwirows, hwicols] vz[hwirows, hwcols] vz[hwirows, wcols] vz[hwirows, dwicols] vz[hwirows, dwcols]
        vz[hwrows, hwicols] vz[hwrows, hwcols] vz[hwrows, wcols] vz[hwrows, dwicols] vz[hwrows, dwcols]
        vz[wrows, hwicols] vz[wrows, hwcols] vz[wrows, wcols] vz[wrows, dwicols] vz[wrows, dwcols]
        vz[dwirows, hwicols] vz[dwirows, hwcols] vz[dwirows, wcols] vz[dwirows, dwicols] vz[dwirows, dwcols]
        vz[dwrows, hwicols] vz[dwrows, hwcols] vz[dwrows, wcols] vz[dwrows, dwicols] vz[dwrows, dwcols]
    ]
    w2wvhat[:, :, 2] = [
        vr[hwirows, hwicols] vr[hwirows, hwcols] vr[hwirows, wcols] vr[hwirows, dwicols] vr[hwirows, dwcols]
        vr[hwrows, hwicols] vr[hwrows, hwcols] vr[hwrows, wcols] vr[hwrows, dwicols] vr[hwrows, dwcols]
        vr[wrows, hwicols] vr[wrows, hwcols] vr[wrows, wcols] vr[wrows, dwicols] vr[wrows, dwcols]
        vr[dwirows, hwicols] vr[dwirows, hwcols] vr[dwirows, wcols] vr[dwirows, dwicols] vr[dwirows, dwcols]
        vr[dwrows, hwicols] vr[dwrows, hwcols] vr[dwrows, wcols] vr[dwrows, dwicols] vr[dwrows, dwcols]
    ]

    # wake to rotor
    w2rvhat = zeros(pidranges[3], nwn, 2)
    w2rvhat[:, :, 1] = [
        vz[rotorrows, hwicols] vz[rotorrows, hwcols] vz[rotorrows, wcols] vz[rotorrows, dwicols] vz[rotorrows, dwcols]
    ]
    w2rvhat[:, :, 2] = [
        vr[rotorrows, hwicols] vr[rotorrows, hwcols] vr[rotorrows, wcols] vr[rotorrows, dwicols] vr[rotorrows, dwcols]
    ]

    return reverse(w2bvhat; dims=1), w2rvhat, w2wvhat
end

function reformat_rotor_velocities(vels)
    return (; r=vels[:, 1], Wz=vels[:, 2], Wr=vels[:, 3], Wm=vels[:, 4], Wtheta=vels[:, 5])
end

function reformat_blade_elements(be)
    return (;
        Wz=be[:,1], # VREL(1, IR, NR),
        Wtheta=be[:,2], # WTB,
        Wmag=be[:,3], # WWB,
        phi=be[:,4], # PHIB,
        xi=be[:,5], # XI,
        alpha=be[:,6], # ALF,
        twist=be[:,7], # BETAR(IR, NR),
        reynolds=be[:,8], # REY,
        solidity=be[:,9], # SECSIG,
        stagger=be[:,10], # SECTAGR,
        cl=be[:,11], # CLR(IR, NR),
        cd=be[:,12], # CDR(IR, NR),
        Gamr_est=be[:,13], # BGX(IR),
    )
end
