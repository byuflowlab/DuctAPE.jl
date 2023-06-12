#=
function to grab all the dfdc extraction file data and put things into a format that is easily accessible
=#

project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

datapath = project_dir * "/examples/dfdc_testing/"

function get_dfdc()

    #####################################
    ## -- Get Panel Indices -- ##
    #####################################
    # element, cp1, cp2, node1, node2
    include(datapath * "DFDC_ELEMENT_INDICES.jl")
    hub_cp_range = idv2r(dfdc_indices[1, 2:3])
    hub_node_range = idv2r(dfdc_indices[1, 4:5])
    duct_cp_range = idv2r(dfdc_indices[2, 2:3])
    duct_node_range = idv2r(dfdc_indices[2, 4:5])
    rotor_cp_range = idv2r(dfdc_indices[3, 2:3])
    rotor_node_range = idv2r(dfdc_indices[3, 4:5])
    hubwake_cp_range = idv2r(dfdc_indices[4, 2:3])
    hubwake_node_range = idv2r(dfdc_indices[4, 4:5])
    ductwake_cp_range = idv2r(dfdc_indices[14, 2:3])
    ductwake_node_range = idv2r(dfdc_indices[14, 4:5])
    wake_cp_range = idv2r(dfdc_indices[5:13, 2:3])
    wake_node_range = idv2r(dfdc_indices[5:13, 4:5])

    #####################################
    ## -- Extract Panel Strengths -- ##
    #####################################
    #columns: node index, gamw, gamb, sigr
    include(datapath * "DFDC_PANEL_STATES.jl")

    # - Get indices on body of non-zero gamma theta - #
    # these are the wake-body interface panels
    # first get non-zero indices for whole column
    nzid = findall(x -> x != 0.0, dfdc_states[:, 2])
    # next get the indices that fall in the hub range
    hwiid = findall(x -> x in hub_node_range, nzid)
    #
    hubwake_interface_node_range = nzid[hwiid]
    # next get the indices that fall in the duct range
    dwiid = findall(x -> x in duct_node_range, nzid)
    #
    ductwake_interface_node_range = nzid[dwiid]

    #TODO: try getting cp ranges from these
    hubwake_interface_cp_range = hub_cp_range[1]:(length(hubwake_interface_node_range) )
    ductwake_interface_cp_range =
        (duct_cp_range[end] - length(ductwake_interface_node_range) +1):duct_cp_range[end]

    # - Body Panel Strengths - #
    gamb_duct = dfdc_states[duct_node_range, 3]
    gamb_hub = dfdc_states[hub_node_range, 3]

    # - Body Wake Strengths - #
    gamw_duct = dfdc_states[ductwake_node_range, 2]
    gamw_hub = dfdc_states[hubwake_node_range, 2]

    # - Body-Wake Interface Strengths - #
    gamw_ductinterface = dfdc_states[ductwake_interface_node_range, 2]
    gamw_hubinterface = dfdc_states[hubwake_interface_node_range, 2]

    # - Rotor Source Strengths - #
    sigr = dfdc_states[rotor_node_range, 4]

    # - Wake Strengths - #
    gamw = zeros(
        length(wake_node_range), wake_node_range[1][end] - wake_node_range[1][1] + 1
    )
    for i in 1:length(wake_node_range)
        gamw[i, :] = dfdc_states[wake_node_range[i], 2]
    end

    #####################################
    ## -- Extract Node Locations -- ##
    #####################################
    # column 2: x location
    # column 3: r location
    include(datapath * "DFDC_ELEMENT_EDGES.jl")
    hub_node_xr = elem1[:, 2:3]
    duct_node_xr = elem2[:, 2:3]
    rotor_node_xr = elem3[:, 2:3]
    hubwake_node_xr = elem4[:, 2:3]
    ductwake_node_xr = elem14[:, 2:3]
    hubinterface_node_xr = elem1[1:length(hubwake_interface_node_range), 2:3]
    ductinterface_node_xr = elem2[
        (end - length(ductwake_interface_node_range) + 1):end, 2:3
    ]
    wake_node_x =
        hcat(
            elem5[:, 2],
            elem6[:, 2],
            elem7[:, 2],
            elem8[:, 2],
            elem9[:, 2],
            elem10[:, 2],
            elem11[:, 2],
            elem12[:, 2],
            elem13[:, 2],
        )'
    wake_node_r =
        hcat(
            elem5[:, 3],
            elem6[:, 3],
            elem7[:, 3],
            elem8[:, 3],
            elem9[:, 3],
            elem10[:, 3],
            elem11[:, 3],
            elem12[:, 3],
            elem13[:, 3],
        )'

    #####################################
    ## -- Extract Control Point Locations -- ##
    #####################################
    # column 2: x location
    # column 3: r location
    # column 6: total Vm
    include(datapath * "DFDC_ELEMENT_CENTERS.jl")
    hub_ctrlpt_xr = elem1[:, 2:3]
    duct_ctrlpt_xr = elem2[:, 2:3]
    rotor_ctrlpt_xr = elem3[:, 2:3]
    hubwake_ctrlpt_xr = elem4[:, 2:3]
    ductwake_ctrlpt_xr = elem14[:, 2:3]
    hubinterface_ctrlpt_xr = elem1[1:length(hubwake_interface_cp_range), 2:3]
    ductinterface_ctrlpt_xr = elem2[
        (end - length(ductwake_interface_cp_range) + 1):end, 2:3
    ]
    wake_ctrlpt_x =
        hcat(
            elem5[:, 2],
            elem6[:, 2],
            elem7[:, 2],
            elem8[:, 2],
            elem9[:, 2],
            elem10[:, 2],
            elem11[:, 2],
            elem12[:, 2],
            elem13[:, 2],
        )'
    wake_ctrlpt_r =
        hcat(
            elem5[:, 3],
            elem6[:, 3],
            elem7[:, 3],
            elem8[:, 3],
            elem9[:, 3],
            elem10[:, 3],
            elem11[:, 3],
            elem12[:, 3],
            elem13[:, 3],
        )'

    # - This file also contains Vm's - #
    Vm_hub = elem1[:, 6]
    Vm_duct = elem2[:, 6]
    Vm_rotor = elem3[:, 6]
    Vm_hubwake = elem4[:, 6]
    Vm_ductwake = elem14[:, 6]
    Vm_wake =
        hcat(
            elem5[:, 6],
            elem6[:, 6],
            elem7[:, 6],
            elem8[:, 6],
            elem9[:, 6],
            elem10[:, 6],
            elem11[:, 6],
            elem12[:, 6],
            elem13[:, 6],
        )'

    #####################################
    ## -- Exract Induced Velocities -- ##
    #####################################

    # 1: index, get from cp ranges
    # 2: vx from body
    # 3: vr from body
    # 4: vx from rotor
    # 5: vr from rotor
    # 6: vx from wake
    # 7: vr from wake
    # 8: vx TOTAL (does not include Vinf)
    # 9: vr TOTAL (does not include Omega*r)
    include(datapath * "DFDC_VELOCITY_BREAKDOWN.jl")

    # body on body
    vx_bb_duct = dfdc_velocities[duct_cp_range, 2]
    vr_bb_duct = dfdc_velocities[duct_cp_range, 3]
    vx_bb_hub = dfdc_velocities[hub_cp_range, 2]
    vr_bb_hub = dfdc_velocities[hub_cp_range, 3]
    # rotor on body
    vx_br_duct = dfdc_velocities[duct_cp_range, 4]
    vr_br_duct = dfdc_velocities[duct_cp_range, 5]
    vx_br_hub = dfdc_velocities[hub_cp_range, 4]
    vr_br_hub = dfdc_velocities[hub_cp_range, 5]
    # wake on body
    vx_bw_duct = dfdc_velocities[duct_cp_range, 6]
    vr_bw_duct = dfdc_velocities[duct_cp_range, 7]
    vx_bw_hub = dfdc_velocities[hub_cp_range, 6]
    vr_bw_hub = dfdc_velocities[hub_cp_range, 7]
    # total on body
    vx_duct = dfdc_velocities[duct_cp_range, 8]
    vr_duct = dfdc_velocities[duct_cp_range, 9]
    vx_hub = dfdc_velocities[hub_cp_range, 8]
    vr_hub = dfdc_velocities[hub_cp_range, 9]

    # body on rotor
    vx_rb = dfdc_velocities[rotor_cp_range, 2]
    vr_rb = dfdc_velocities[rotor_cp_range, 3]
    # rotor on rotor
    vx_rr = dfdc_velocities[rotor_cp_range, 4]
    vr_rr = dfdc_velocities[rotor_cp_range, 5]
    # wake on rotor
    vx_rw = dfdc_velocities[rotor_cp_range, 6]
    vr_rw = dfdc_velocities[rotor_cp_range, 7]
    # total on rotor
    vx_rotor = dfdc_velocities[rotor_cp_range, 8]
    vr_rotor = dfdc_velocities[rotor_cp_range, 9]

    #TODO: get interface and body wake values
    # wake on hub interface
    vx_hiw = dfdc_velocities[hubwake_interface_cp_range, 6]
    vr_hiw = dfdc_velocities[hubwake_interface_cp_range, 7]
    # wake on duct interface
    vx_diw = dfdc_velocities[ductwake_interface_cp_range, 6]
    vr_diw = dfdc_velocities[ductwake_interface_cp_range, 7]
    # wake on hub wake
    vx_hww = dfdc_velocities[hubwake_cp_range, 6]
    vr_hww = dfdc_velocities[hubwake_cp_range, 7]
    # wake on duct wake
    vx_dww = dfdc_velocities[ductwake_cp_range, 6]
    vr_dww = dfdc_velocities[ductwake_cp_range, 7]

    # body on wake
    vx_wb = zeros(length(wake_cp_range), wake_cp_range[1][end] - wake_cp_range[1][1] + 1)
    vr_wb = similar(vx_wb) .= 0.0
    # rotor on wake
    vx_wr = similar(vx_wb) .= 0.0
    vr_wr = similar(vx_wb) .= 0.0
    # wake on wake
    vx_ww = similar(vx_wb) .= 0.0
    vr_ww = similar(vx_wb) .= 0.0
    # total on wake
    vx_wake = similar(vx_wb) .= 0.0
    vr_wake = similar(vx_wb) .= 0.0

    for i in 1:length(wake_node_range)
        # body on wake
        vx_wb[i, :] = dfdc_velocities[wake_cp_range[i], 2]
        vr_wb[i, :] = dfdc_velocities[wake_cp_range[i], 3]
        # rotor on wake
        vx_wr[i, :] = dfdc_velocities[wake_cp_range[i], 4]
        vr_wr[i, :] = dfdc_velocities[wake_cp_range[i], 5]
        # wake on wake
        vx_ww[i, :] = dfdc_velocities[wake_cp_range[i], 6]
        vr_ww[i, :] = dfdc_velocities[wake_cp_range[i], 7]
        # total on wake
        vx_wake[i, :] = dfdc_velocities[wake_cp_range[i], 8]
        vr_wake[i, :] = dfdc_velocities[wake_cp_range[i], 9]
    end

    ## -- Extract Blade Element Data -- ##
    include(datapath * "ALL_THE_STUFF.jl")
    # dfdc_blade_elements
    #c,r,aoa,twist,PHI,W,Re,solidity,STAGR,CL,CD
    include(datapath * "DFDC_BLADE_ELEMENT_INFO.jl")

    dfdc = (;
        # - Indices
        hub_cp_range,
        hub_node_range,
        duct_cp_range,
        duct_node_range,
        rotor_cp_range,
        rotor_node_range,
        hubwake_cp_range,
        hubwake_node_range,
        ductwake_cp_range,
        ductwake_node_range,
        hubwake_interface_node_range,
        hubwake_interface_cp_range,
        ductwake_interface_node_range,
        ductwake_interface_cp_range,
        wake_cp_range,
        wake_node_range,
        # - Locations
        duct_node_x=duct_node_xr[:, 1],
        duct_node_r=duct_node_xr[:, 2],
        duct_ctrlpt_x=duct_ctrlpt_xr[:, 1],
        duct_ctrlpt_r=duct_ctrlpt_xr[:, 2],
        hub_node_x=hub_node_xr[:, 1],
        hub_node_r=hub_node_xr[:, 2],
        hub_ctrlpt_x=hub_ctrlpt_xr[:, 1],
        hub_ctrlpt_r=hub_ctrlpt_xr[:, 2],
        ductwake_node_x=ductwake_node_xr[:, 1],
        ductwake_node_r=ductwake_node_xr[:, 2],
        ductwake_ctrlpt_x=ductwake_ctrlpt_xr[:, 1],
        ductwake_ctrlpt_r=ductwake_ctrlpt_xr[:, 2],
        hubwake_node_x=hubwake_node_xr[:, 1],
        hubwake_node_r=hubwake_node_xr[:, 2],
        hubwake_ctrlpt_x=hubwake_ctrlpt_xr[:, 1],
        hubwake_ctrlpt_r=hubwake_ctrlpt_xr[:, 2],
        ductinterface_node_x=ductinterface_node_xr[:, 1],
        ductinterface_node_r=ductinterface_node_xr[:, 2],
        ductinterface_ctrlpt_x=ductinterface_ctrlpt_xr[:, 1],
        ductinterface_ctrlpt_r=ductinterface_ctrlpt_xr[:, 2],
        hubinterface_node_x=hubinterface_node_xr[:, 1],
        hubinterface_node_r=hubinterface_node_xr[:, 2],
        hubinterface_ctrlpt_x=hubinterface_ctrlpt_xr[:, 1],
        hubinterface_ctrlpt_r=hubinterface_ctrlpt_xr[:, 2],
        rotor_node_x=rotor_node_xr[:, 1],
        rotor_node_r=rotor_node_xr[:, 2],
        rotor_ctrlpt_x=rotor_ctrlpt_xr[:, 1],
        rotor_ctrlpt_r=rotor_ctrlpt_xr[:, 2],
        wake_node_x,
        wake_node_r,
        wake_ctrlpt_x,
        wake_ctrlpt_r,
        # - States
        gamb_duct,
        gamb_hub,
        gamw_duct, #duct wake
        gamw_hub, #hub wake
        gamw_ductinterface, #gamma theta on duct-wake interface
        gamw_hubinterface, #gamma theta on hub-wake interface
        sigr,
        gamw,
        Gamr=dfdcGamr,
        # - Induced Velocities
        # body on body
        vx_bb_duct,
        vr_bb_duct,
        vx_bb_hub,
        vr_bb_hub,
        # rotor on body
        vx_br_duct,
        vr_br_duct,
        vx_br_hub,
        vr_br_hub,
        # wake on body
        vx_bw_duct,
        vr_bw_duct,
        vx_bw_hub,
        vr_bw_hub,
        # total on body
        vx_duct,
        vr_duct,
        vx_hub,
        vr_hub,
        # body on rotor
        vx_rb,
        vr_rb,
        # rotor on rotor
        vx_rr,
        vr_rr,
        # wake on rotor
        vx_rw,
        vr_rw,
        # total on rotor
        vx_rotor,
        vr_rotor,
        # body on wake
        vx_wb,
        vr_wb,
        # rotor on wake
        vx_wr,
        vr_wr,
        # wake on wake
        vx_ww,
        vr_ww,
        # total on wake
        vx_wake,
        vr_wake,
        # wake on hub interface
        vx_hiw,
        vr_hiw,
        # wake on duct interface
        vx_diw,
        vr_diw,
        # wake on hub wake
        vx_hww,
        vr_hww,
        # wake on duct wake
        vx_dww,
        vr_dww,
        # - Meridional Velocities
        Vm_hub,
        Vm_duct,
        Vm_rotor,
        Vm_hubwake,
        Vm_ductwake,
        Vm_wake,
        # - Blade Element Data
        inflow_angle=dfdc_blade_elements[:, 5],
        alpha=dfdc_blade_elements[:, 3],
        clift=dfdc_blade_elements[:, 10],
        cdrag=dfdc_blade_elements[:, 11],
        twist=dfdc_blade_elements[:, 4],
        chord=dfdc_blade_elements[:, 1],
        solidity=dfdc_blade_elements[:, 8],
        stagger=dfdc_blade_elements[:, 9],
        deltaH,
        deltaS,
        Wmag_rotor=dfdc_blade_elements[:, 6],
    )

    return dfdc
end

"""
index vector to range
"""
function idv2r(vec)
    if typeof(vec) <: AbstractVector
        return vec[1]:vec[2]
    elseif typeof(vec) <: AbstractMatrix
        ranges = []
        for i in 1:length(vec[:, 1])
            push!(ranges, vec[i, 1]:vec[i, 2])
        end
        return ranges
    end

    return nothing
end
