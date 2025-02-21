function withdraw_prepost_body_container_cache(vec, dims)
    zpts = (;
        centerbody_zpts=reshape(
            @view(vec[dims.zpts.centerbody_zpts.index]), dims.zpts.centerbody_zpts.shape
        ),
        casing_zpts=reshape(
            @view(vec[dims.zpts.casing_zpts.index]), dims.zpts.casing_zpts.shape
        ),
        nacelle_zpts=reshape(
            @view(vec[dims.zpts.nacelle_zpts.index]), dims.zpts.nacelle_zpts.shape
        ),
    )

    vtan_tuple = (;
        Vtot_in=reshape(
            @view(vec[dims.vtan_tuple.Vtot_in.index]), dims.vtan_tuple.Vtot_in.shape
        ),
        Vtot_out=reshape(
            @view(vec[dims.vtan_tuple.Vtot_out.index]), dims.vtan_tuple.Vtot_out.shape
        ),
        Vtan_in=reshape(
            @view(vec[dims.vtan_tuple.Vtan_in.index]), dims.vtan_tuple.Vtan_in.shape
        ),
        Vtan_out=reshape(
            @view(vec[dims.vtan_tuple.Vtan_out.index]), dims.vtan_tuple.Vtan_out.shape
        ),
        Vtot_prejump=reshape(
            @view(vec[dims.vtan_tuple.Vtot_prejump.index]),
            dims.vtan_tuple.Vtot_prejump.shape,
        ),
        vtot_body=reshape(
            @view(vec[dims.vtan_tuple.vtot_body.index]), dims.vtan_tuple.vtot_body.shape
        ),
        duct_jump=reshape(
            @view(vec[dims.vtan_tuple.duct_jump.index]), dims.vtan_tuple.duct_jump.shape
        ),
        centerbody_jump=reshape(
            @view(vec[dims.vtan_tuple.centerbody_jump.index]),
            dims.vtan_tuple.centerbody_jump.shape,
        ),
        body_jump_term=reshape(
            @view(vec[dims.vtan_tuple.body_jump_term.index]),
            dims.vtan_tuple.body_jump_term.shape,
        ),
        vtot_jump=reshape(
            @view(vec[dims.vtan_tuple.vtot_jump.index]), dims.vtan_tuple.vtot_jump.shape
        ),
        vtot_wake=reshape(
            @view(vec[dims.vtan_tuple.vtot_wake.index]), dims.vtan_tuple.vtot_wake.shape
        ),
        vtot_rotors=reshape(
            @view(vec[dims.vtan_tuple.vtot_rotors.index]), dims.vtan_tuple.vtot_rotors.shape
        ),
        # Splits:
        vtan_casing_in=reshape(
            @view(vec[dims.vtan_tuple.vtan_casing_in.index]),
            dims.vtan_tuple.vtan_casing_in.shape,
        ),
        vtan_casing_out=reshape(
            @view(vec[dims.vtan_tuple.vtan_casing_out.index]),
            dims.vtan_tuple.vtan_casing_out.shape,
        ),
        vtan_nacelle_in=reshape(
            @view(vec[dims.vtan_tuple.vtan_nacelle_in.index]),
            dims.vtan_tuple.vtan_nacelle_in.shape,
        ),
        vtan_nacelle_out=reshape(
            @view(vec[dims.vtan_tuple.vtan_nacelle_out.index]),
            dims.vtan_tuple.vtan_nacelle_out.shape,
        ),
        vtan_centerbody_in=reshape(
            @view(vec[dims.vtan_tuple.vtan_centerbody_in.index]),
            dims.vtan_tuple.vtan_centerbody_in.shape,
        ),
        vtan_centerbody_out=reshape(
            @view(vec[dims.vtan_tuple.vtan_centerbody_out.index]),
            dims.vtan_tuple.vtan_centerbody_out.shape,
        ),
    )

    cp_tuple = (;
        cp_in=reshape(@view(vec[dims.cp_tuple.cp_in.index]), dims.cp_tuple.cp_in.shape),
        cp_out=reshape(@view(vec[dims.cp_tuple.cp_out.index]), dims.cp_tuple.cp_out.shape),
        cp_casing_in=reshape(
            @view(vec[dims.cp_tuple.cp_casing_in.index]), dims.cp_tuple.cp_casing_in.shape
        ),
        cp_casing_out=reshape(
            @view(vec[dims.cp_tuple.cp_casing_out.index]), dims.cp_tuple.cp_casing_out.shape
        ),
        cp_nacelle_in=reshape(
            @view(vec[dims.cp_tuple.cp_nacelle_in.index]), dims.cp_tuple.cp_nacelle_in.shape
        ),
        cp_nacelle_out=reshape(
            @view(vec[dims.cp_tuple.cp_nacelle_out.index]),
            dims.cp_tuple.cp_nacelle_out.shape,
        ),
        cp_centerbody_in=reshape(
            @view(vec[dims.cp_tuple.cp_centerbody_in.index]),
            dims.cp_tuple.cp_centerbody_in.shape,
        ),
        cp_centerbody_out=reshape(
            @view(vec[dims.cp_tuple.cp_centerbody_out.index]),
            dims.cp_tuple.cp_centerbody_out.shape,
        ),
    )

    body_inviscid_thrust = reshape(@view(vec[dims.body_inviscid_thrust.index]), dims.body_inviscid_thrust.shape)

    body_viscous_drag = reshape(@view(vec[dims.body_viscous_drag.index]), dims.body_viscous_drag.shape)

    body_thrust = reshape(@view(vec[dims.body_thrust.index]), dims.body_thrust.shape)

    body_force_coefficient = reshape(
        @view(vec[dims.body_force_coefficient.index]), dims.body_force_coefficient.shape
    )

    return zpts, vtan_tuple, cp_tuple, body_inviscid_thrust, body_viscous_drag,body_thrust, body_force_coefficient
end

"""
    withdraw_prepost_container_cache(vec, dims)

Reshape the prepost cache vector using the saved dimensions tuple.

# Arguments
- `vec::Vector{Float}` : vector cache of pre- and post-processing intermediate containers.
- `dims::NamedTuple` : Named tuple containing the indices and shape of the various items stored in the cache vector.

# Returns
- `prepost_container_caching::NamedTuple` : Named tuple containing reshaped views of sections of the cache vector.
"""
function withdraw_prepost_container_cache(vec, dims)
    wake_grid = reshape(@view(vec[dims.wake_grid.index]), dims.wake_grid.shape)
    rp_duct_coordinates = reshape(
        @view(vec[dims.rp_duct_coordinates.index]), dims.rp_duct_coordinates.shape
    )
    rp_centerbody_coordinates = reshape(
        @view(vec[dims.rp_centerbody_coordinates.index]),
        dims.rp_centerbody_coordinates.shape,
    )
    rotor_indices_in_wake = reshape(
        @view(vec[dims.rotor_indices_in_wake.index]), dims.rotor_indices_in_wake.shape
    )
    AICn = reshape(@view(vec[dims.AICn.index]), dims.AICn.shape)
    AICpcp = reshape(@view(vec[dims.AICpcp.index]), dims.AICpcp.shape)
    vdnb = reshape(@view(vec[dims.vdnb.index]), dims.vdnb.shape)
    vdnpcp = reshape(@view(vec[dims.vdnpcp.index]), dims.vdnpcp.shape)

    ivb = (;
        v_bb=reshape(@view(vec[dims.ivb.v_bb.index]), dims.ivb.v_bb.shape),
        v_br=reshape(@view(vec[dims.ivb.v_br.index]), dims.ivb.v_br.shape),
        v_bw=reshape(@view(vec[dims.ivb.v_bw.index]), dims.ivb.v_bw.shape),
    )

    panels = (;
        body_vortex_panels=(;
            nnode=reshape(
                @view(vec[dims.panels.body_vortex_panels.nnode.index]),
                dims.panels.body_vortex_panels.nnode.shape,
            ),
            npanel=reshape(
                @view(vec[dims.panels.body_vortex_panels.npanel.index]),
                dims.panels.body_vortex_panels.npanel.shape,
            ),
            nbodies=reshape(
                @view(vec[dims.panels.body_vortex_panels.nbodies.index]),
                dims.panels.body_vortex_panels.nbodies.shape,
            ),
            totpanel=reshape(
                @view(vec[dims.panels.body_vortex_panels.totpanel.index]),
                dims.panels.body_vortex_panels.totpanel.shape,
            ),
            totnode=reshape(
                @view(vec[dims.panels.body_vortex_panels.totnode.index]),
                dims.panels.body_vortex_panels.totnode.shape,
            ),
            node=reshape(
                @view(vec[dims.panels.body_vortex_panels.node.index]),
                dims.panels.body_vortex_panels.node.shape,
            ),
            controlpoint=reshape(
                @view(vec[dims.panels.body_vortex_panels.controlpoint.index]),
                dims.panels.body_vortex_panels.controlpoint.shape,
            ),
            normal=reshape(
                @view(vec[dims.panels.body_vortex_panels.normal.index]),
                dims.panels.body_vortex_panels.normal.shape,
            ),
            tangent=reshape(
                @view(vec[dims.panels.body_vortex_panels.tangent.index]),
                dims.panels.body_vortex_panels.tangent.shape,
            ),
            nodemap=reshape(
                @view(vec[dims.panels.body_vortex_panels.nodemap.index]),
                dims.panels.body_vortex_panels.nodemap.shape,
            ),
            influence_length=reshape(
                @view(vec[dims.panels.body_vortex_panels.influence_length.index]),
                dims.panels.body_vortex_panels.influence_length.shape,
            ),
            endnodes=reshape(
                @view(vec[dims.panels.body_vortex_panels.endnodes.index]),
                dims.panels.body_vortex_panels.endnodes.shape,
            ),
            tenode=reshape(
                @view(vec[dims.panels.body_vortex_panels.tenode.index]),
                dims.panels.body_vortex_panels.tenode.shape,
            ),
            itcontrolpoint=reshape(
                @view(vec[dims.panels.body_vortex_panels.itcontrolpoint.index]),
                dims.panels.body_vortex_panels.itcontrolpoint.shape,
            ),
            itnormal=reshape(
                @view(vec[dims.panels.body_vortex_panels.itnormal.index]),
                dims.panels.body_vortex_panels.itnormal.shape,
            ),
            ittangent=reshape(
                @view(vec[dims.panels.body_vortex_panels.ittangent.index]),
                dims.panels.body_vortex_panels.ittangent.shape,
            ),
            tenormal=reshape(
                @view(vec[dims.panels.body_vortex_panels.tenormal.index]),
                dims.panels.body_vortex_panels.tenormal.shape,
            ),
            tendotn=reshape(
                @view(vec[dims.panels.body_vortex_panels.tendotn.index]),
                dims.panels.body_vortex_panels.tendotn.shape,
            ),
            tencrossn=reshape(
                @view(vec[dims.panels.body_vortex_panels.tencrossn.index]),
                dims.panels.body_vortex_panels.tencrossn.shape,
            ),
            endnodeidxs=reshape(
                @view(vec[dims.panels.body_vortex_panels.endnodeidxs.index]),
                dims.panels.body_vortex_panels.endnodeidxs.shape,
            ),
            endpanelidxs=reshape(
                @view(vec[dims.panels.body_vortex_panels.endpanelidxs.index]),
                dims.panels.body_vortex_panels.endpanelidxs.shape,
            ),
            teadjnodeidxs=reshape(
                @view(vec[dims.panels.body_vortex_panels.teadjnodeidxs.index]),
                dims.panels.body_vortex_panels.teadjnodeidxs.shape,
            ),
            teinfluence_length=reshape(
                @view(vec[dims.panels.body_vortex_panels.teinfluence_length.index]),
                dims.panels.body_vortex_panels.teinfluence_length.shape,
            ),
            prescribednodeidxs=reshape(
                @view(vec[dims.panels.body_vortex_panels.prescribednodeidxs.index]),
                dims.panels.body_vortex_panels.prescribednodeidxs.shape,
            ),
        ),
        rotor_source_panels=(;
            nnode=reshape(
                @view(vec[dims.panels.rotor_source_panels.nnode.index]),
                dims.panels.rotor_source_panels.nnode.shape,
            ),
            npanel=reshape(
                @view(vec[dims.panels.rotor_source_panels.npanel.index]),
                dims.panels.rotor_source_panels.npanel.shape,
            ),
            nbodies=reshape(
                @view(vec[dims.panels.rotor_source_panels.nbodies.index]),
                dims.panels.rotor_source_panels.nbodies.shape,
            ),
            totpanel=reshape(
                @view(vec[dims.panels.rotor_source_panels.totpanel.index]),
                dims.panels.rotor_source_panels.totpanel.shape,
            ),
            totnode=reshape(
                @view(vec[dims.panels.rotor_source_panels.totnode.index]),
                dims.panels.rotor_source_panels.totnode.shape,
            ),
            node=reshape(
                @view(vec[dims.panels.rotor_source_panels.node.index]),
                dims.panels.rotor_source_panels.node.shape,
            ),
            controlpoint=reshape(
                @view(vec[dims.panels.rotor_source_panels.controlpoint.index]),
                dims.panels.rotor_source_panels.controlpoint.shape,
            ),
            normal=reshape(
                @view(vec[dims.panels.rotor_source_panels.normal.index]),
                dims.panels.rotor_source_panels.normal.shape,
            ),
            tangent=reshape(
                @view(vec[dims.panels.rotor_source_panels.tangent.index]),
                dims.panels.rotor_source_panels.tangent.shape,
            ),
            nodemap=reshape(
                @view(vec[dims.panels.rotor_source_panels.nodemap.index]),
                dims.panels.rotor_source_panels.nodemap.shape,
            ),
            influence_length=reshape(
                @view(vec[dims.panels.rotor_source_panels.influence_length.index]),
                dims.panels.rotor_source_panels.influence_length.shape,
            ),
            endnodes=reshape(
                @view(vec[dims.panels.rotor_source_panels.endnodes.index]),
                dims.panels.rotor_source_panels.endnodes.shape,
            ),
            tenode=reshape(
                @view(vec[dims.panels.rotor_source_panels.tenode.index]),
                dims.panels.rotor_source_panels.tenode.shape,
            ),
            itcontrolpoint=reshape(
                @view(vec[dims.panels.rotor_source_panels.itcontrolpoint.index]),
                dims.panels.rotor_source_panels.itcontrolpoint.shape,
            ),
            itnormal=reshape(
                @view(vec[dims.panels.rotor_source_panels.itnormal.index]),
                dims.panels.rotor_source_panels.itnormal.shape,
            ),
            ittangent=reshape(
                @view(vec[dims.panels.rotor_source_panels.ittangent.index]),
                dims.panels.rotor_source_panels.ittangent.shape,
            ),
            tenormal=reshape(
                @view(vec[dims.panels.rotor_source_panels.tenormal.index]),
                dims.panels.rotor_source_panels.tenormal.shape,
            ),
            tendotn=reshape(
                @view(vec[dims.panels.rotor_source_panels.tendotn.index]),
                dims.panels.rotor_source_panels.tendotn.shape,
            ),
            tencrossn=reshape(
                @view(vec[dims.panels.rotor_source_panels.tencrossn.index]),
                dims.panels.rotor_source_panels.tencrossn.shape,
            ),
            endnodeidxs=reshape(
                @view(vec[dims.panels.rotor_source_panels.endnodeidxs.index]),
                dims.panels.rotor_source_panels.endnodeidxs.shape,
            ),
            endpanelidxs=reshape(
                @view(vec[dims.panels.rotor_source_panels.endpanelidxs.index]),
                dims.panels.rotor_source_panels.endpanelidxs.shape,
            ),
            teadjnodeidxs=reshape(
                @view(vec[dims.panels.rotor_source_panels.teadjnodeidxs.index]),
                dims.panels.rotor_source_panels.teadjnodeidxs.shape,
            ),
            teinfluence_length=reshape(
                @view(vec[dims.panels.rotor_source_panels.teinfluence_length.index]),
                dims.panels.rotor_source_panels.teinfluence_length.shape,
            ),
            prescribednodeidxs=reshape(
                @view(vec[dims.panels.rotor_source_panels.prescribednodeidxs.index]),
                dims.panels.rotor_source_panels.prescribednodeidxs.shape,
            ),
        ),
        wake_vortex_panels=(;
            nnode=reshape(
                @view(vec[dims.panels.wake_vortex_panels.nnode.index]),
                dims.panels.wake_vortex_panels.nnode.shape,
            ),
            npanel=reshape(
                @view(vec[dims.panels.wake_vortex_panels.npanel.index]),
                dims.panels.wake_vortex_panels.npanel.shape,
            ),
            nbodies=reshape(
                @view(vec[dims.panels.wake_vortex_panels.nbodies.index]),
                dims.panels.wake_vortex_panels.nbodies.shape,
            ),
            totpanel=reshape(
                @view(vec[dims.panels.wake_vortex_panels.totpanel.index]),
                dims.panels.wake_vortex_panels.totpanel.shape,
            ),
            totnode=reshape(
                @view(vec[dims.panels.wake_vortex_panels.totnode.index]),
                dims.panels.wake_vortex_panels.totnode.shape,
            ),
            node=reshape(
                @view(vec[dims.panels.wake_vortex_panels.node.index]),
                dims.panels.wake_vortex_panels.node.shape,
            ),
            controlpoint=reshape(
                @view(vec[dims.panels.wake_vortex_panels.controlpoint.index]),
                dims.panels.wake_vortex_panels.controlpoint.shape,
            ),
            normal=reshape(
                @view(vec[dims.panels.wake_vortex_panels.normal.index]),
                dims.panels.wake_vortex_panels.normal.shape,
            ),
            tangent=reshape(
                @view(vec[dims.panels.wake_vortex_panels.tangent.index]),
                dims.panels.wake_vortex_panels.tangent.shape,
            ),
            nodemap=reshape(
                @view(vec[dims.panels.wake_vortex_panels.nodemap.index]),
                dims.panels.wake_vortex_panels.nodemap.shape,
            ),
            influence_length=reshape(
                @view(vec[dims.panels.wake_vortex_panels.influence_length.index]),
                dims.panels.wake_vortex_panels.influence_length.shape,
            ),
            endnodes=reshape(
                @view(vec[dims.panels.wake_vortex_panels.endnodes.index]),
                dims.panels.wake_vortex_panels.endnodes.shape,
            ),
            tenode=reshape(
                @view(vec[dims.panels.wake_vortex_panels.tenode.index]),
                dims.panels.wake_vortex_panels.tenode.shape,
            ),
            itcontrolpoint=reshape(
                @view(vec[dims.panels.wake_vortex_panels.itcontrolpoint.index]),
                dims.panels.wake_vortex_panels.itcontrolpoint.shape,
            ),
            itnormal=reshape(
                @view(vec[dims.panels.wake_vortex_panels.itnormal.index]),
                dims.panels.wake_vortex_panels.itnormal.shape,
            ),
            ittangent=reshape(
                @view(vec[dims.panels.wake_vortex_panels.ittangent.index]),
                dims.panels.wake_vortex_panels.ittangent.shape,
            ),
            tenormal=reshape(
                @view(vec[dims.panels.wake_vortex_panels.tenormal.index]),
                dims.panels.wake_vortex_panels.tenormal.shape,
            ),
            tendotn=reshape(
                @view(vec[dims.panels.wake_vortex_panels.tendotn.index]),
                dims.panels.wake_vortex_panels.tendotn.shape,
            ),
            tencrossn=reshape(
                @view(vec[dims.panels.wake_vortex_panels.tencrossn.index]),
                dims.panels.wake_vortex_panels.tencrossn.shape,
            ),
            endnodeidxs=reshape(
                @view(vec[dims.panels.wake_vortex_panels.endnodeidxs.index]),
                dims.panels.wake_vortex_panels.endnodeidxs.shape,
            ),
            endpanelidxs=reshape(
                @view(vec[dims.panels.wake_vortex_panels.endpanelidxs.index]),
                dims.panels.wake_vortex_panels.endpanelidxs.shape,
            ),
            teadjnodeidxs=reshape(
                @view(vec[dims.panels.wake_vortex_panels.teadjnodeidxs.index]),
                dims.panels.wake_vortex_panels.teadjnodeidxs.shape,
            ),
            teinfluence_length=reshape(
                @view(vec[dims.panels.wake_vortex_panels.teinfluence_length.index]),
                dims.panels.wake_vortex_panels.teinfluence_length.shape,
            ),
            prescribednodeidxs=reshape(
                @view(vec[dims.panels.wake_vortex_panels.prescribednodeidxs.index]),
                dims.panels.wake_vortex_panels.prescribednodeidxs.shape,
            ),
        ),
    )

    # - ROTOR POST CACHE - #

    rotor_inviscid_thrust = reshape(
        @view(vec[dims.rotor_inviscid_thrust.index]), dims.rotor_inviscid_thrust.shape
    )
    rotor_inviscid_thrust_dist = reshape(
        @view(vec[dims.rotor_inviscid_thrust_dist.index]),
        dims.rotor_inviscid_thrust_dist.shape,
    )
    rotor_viscous_thrust = reshape(
        @view(vec[dims.rotor_viscous_thrust.index]), dims.rotor_viscous_thrust.shape
    )
    rotor_viscous_thrust_dist = reshape(
        @view(vec[dims.rotor_viscous_thrust_dist.index]),
        dims.rotor_viscous_thrust_dist.shape,
    )
    rotor_thrust = reshape(@view(vec[dims.rotor_thrust.index]), dims.rotor_thrust.shape)
    rotor_inviscid_torque = reshape(
        @view(vec[dims.rotor_inviscid_torque.index]), dims.rotor_inviscid_torque.shape
    )
    rotor_inviscid_torque_dist = reshape(
        @view(vec[dims.rotor_inviscid_torque_dist.index]),
        dims.rotor_inviscid_torque_dist.shape,
    )
    rotor_viscous_torque = reshape(
        @view(vec[dims.rotor_viscous_torque.index]), dims.rotor_viscous_torque.shape
    )
    rotor_viscous_torque_dist = reshape(
        @view(vec[dims.rotor_viscous_torque_dist.index]),
        dims.rotor_viscous_torque_dist.shape,
    )
    rotor_torque = reshape(@view(vec[dims.rotor_torque.index]), dims.rotor_torque.shape)
    rotor_inviscid_power = reshape(
        @view(vec[dims.rotor_inviscid_power.index]), dims.rotor_inviscid_power.shape
    )
    rotor_inviscid_power_dist = reshape(
        @view(vec[dims.rotor_inviscid_power_dist.index]),
        dims.rotor_inviscid_power_dist.shape,
    )
    rotor_viscous_power = reshape(
        @view(vec[dims.rotor_viscous_power.index]), dims.rotor_viscous_power.shape
    )
    rotor_viscous_power_dist = reshape(
        @view(vec[dims.rotor_viscous_power_dist.index]), dims.rotor_viscous_power_dist.shape
    )
    rotor_power = reshape(@view(vec[dims.rotor_power.index]), dims.rotor_power.shape)
    rotor_CT = reshape(@view(vec[dims.rotor_CT.index]), dims.rotor_CT.shape)
    rotor_CQ = reshape(@view(vec[dims.rotor_CQ.index]), dims.rotor_CQ.shape)
    rotor_CP = reshape(@view(vec[dims.rotor_CP.index]), dims.rotor_CP.shape)
    rotor_efficiency = reshape(
        @view(vec[dims.rotor_efficiency.index]), dims.rotor_efficiency.shape
    )
    induced_efficiency = reshape(
        @view(vec[dims.induced_efficiency.index]), dims.induced_efficiency.shape
    )
    blade_normal_force_per_unit_span = reshape(
        @view(vec[dims.blade_normal_force_per_unit_span.index]),
        dims.blade_normal_force_per_unit_span.shape,
    )
    blade_tangential_force_per_unit_span = reshape(
        @view(vec[dims.blade_tangential_force_per_unit_span.index]),
        dims.blade_tangential_force_per_unit_span.shape,
    )
    blade_loading_intermediate_containers = (;
        cn=reshape(
            @view(vec[dims.blade_loading_intermediate_containers.cn.index]),
            dims.blade_loading_intermediate_containers.cn.shape,
        ),
        ct=reshape(
            @view(vec[dims.blade_loading_intermediate_containers.ct.index]),
            dims.blade_loading_intermediate_containers.ct.shape,
        ),
        cphi=reshape(
            @view(vec[dims.blade_loading_intermediate_containers.cphi.index]),
            dims.blade_loading_intermediate_containers.cphi.shape,
        ),
        sphi=reshape(
            @view(vec[dims.blade_loading_intermediate_containers.sphi.index]),
            dims.blade_loading_intermediate_containers.sphi.shape,
        ),
    )

    # - BODY POST CACHE - #
    zpts, vtan_tuple, cp_tuple, body_inviscid_thrust, body_viscous_drag, body_thrust, body_force_coefficient = withdraw_prepost_body_container_cache(
        vec, dims
    )

    # - TOTALS POST CACHE - #
    total_thrust = reshape(@view(vec[dims.total_thrust.index]), dims.total_thrust.shape)
    total_torque = reshape(@view(vec[dims.total_torque.index]), dims.total_torque.shape)
    total_power = reshape(@view(vec[dims.total_power.index]), dims.total_power.shape)
    total_CT = reshape(@view(vec[dims.total_CT.index]), dims.total_CT.shape)
    total_CQ = reshape(@view(vec[dims.total_CQ.index]), dims.total_CQ.shape)
    total_CP = reshape(@view(vec[dims.total_CP.index]), dims.total_CP.shape)
    total_efficiency = reshape(
        @view(vec[dims.total_efficiency.index]), dims.total_efficiency.shape
    )
    ideal_efficiency = reshape(
        @view(vec[dims.ideal_efficiency.index]), dims.ideal_efficiency.shape
    )

    return (;
        # pre
        rp_duct_coordinates,
        rp_centerbody_coordinates,
        wake_grid,
        rotor_indices_in_wake,
        AICn,
        AICpcp,
        vdnb,
        vdnpcp,
        panels,
        ivb,
        #rotor post
        rotor_inviscid_thrust,
        rotor_inviscid_thrust_dist,
        rotor_viscous_thrust,
        rotor_viscous_thrust_dist,
        rotor_thrust,
        rotor_inviscid_torque,
        rotor_inviscid_torque_dist,
        rotor_viscous_torque,
        rotor_viscous_torque_dist,
        rotor_torque,
        rotor_inviscid_power,
        rotor_inviscid_power_dist,
        rotor_viscous_power,
        rotor_viscous_power_dist,
        rotor_power,
        rotor_CT,
        rotor_CQ,
        rotor_CP,
        rotor_efficiency,
        induced_efficiency,
        blade_normal_force_per_unit_span,
        blade_tangential_force_per_unit_span,
        blade_loading_intermediate_containers,
        # Body Post
        zpts,
        vtan_tuple,
        cp_tuple,
        body_inviscid_thrust,
        body_viscous_drag,
        body_thrust,
        body_force_coefficient,
        # Totals Post
        total_thrust,
        total_torque,
        total_power,
        total_CT,
        total_CQ,
        total_CP,
        total_efficiency,
        ideal_efficiency,
    )
end

"""
    withdraw_solve_parameter_cache(solver_options::SolverOptionsType, vec, dims)

Reshape the solve parameter cache vector using the saved dimensions tuple.

# Arguments
- `solver_options::SolverOptionsType` : Solver options type for dispatch.
- `vec::Vector{Float}` : vector cache of pre- and post-processing intermediate containers.
- `dims::NamedTuple` : Named tuple containing the indices and shape of the various items stored in the cache vector.

# Returns
- `solve_parameter_caching::NamedTuple` : Named tuple containing reshaped views of sections of the cache vector.
"""
function withdraw_solve_parameter_cache(solver_options::TS, vec, dims) where {TS<:InternalSolverOptions}

    # - Initial Guesses - #
    Gamr = reshape(@view(vec[dims.Gamr.index]), dims.Gamr.shape)
    sigr = reshape(@view(vec[dims.sigr.index]), dims.sigr.shape)
    gamw = reshape(@view(vec[dims.gamw.index]), dims.gamw.shape)

    # - Operating Point - #
    operating_point = (;
        Vinf=reshape(
            @view(vec[dims.operating_point.Vinf.index]), dims.operating_point.Vinf.shape
        ),
        Minf=reshape(
            @view(vec[dims.operating_point.Minf.index]), dims.operating_point.Minf.shape
        ),
        rhoinf=reshape(
            @view(vec[dims.operating_point.rhoinf.index]), dims.operating_point.rhoinf.shape
        ),
        muinf=reshape(
            @view(vec[dims.operating_point.muinf.index]), dims.operating_point.muinf.shape
        ),
        asound=reshape(
            @view(vec[dims.operating_point.asound.index]), dims.operating_point.asound.shape
        ),
        Ptot=reshape(
            @view(vec[dims.operating_point.Ptot.index]), dims.operating_point.Ptot.shape
        ),
        Ttot=reshape(
            @view(vec[dims.operating_point.Ttot.index]), dims.operating_point.Ttot.shape
        ),
        Omega=reshape(
            @view(vec[dims.operating_point.Omega.index]), dims.operating_point.Omega.shape
        ),
    )

    # - induced velocities on rotor - #
    ivr = (;
        v_rb=reshape(@view(vec[dims.ivr.v_rb.index]), dims.ivr.v_rb.shape),
        v_rr=reshape(@view(vec[dims.ivr.v_rr.index]), dims.ivr.v_rr.shape),
        v_rw=reshape(@view(vec[dims.ivr.v_rw.index]), dims.ivr.v_rw.shape),
    )

    # - induced velocities on wake - #
    ivw = (;
        v_wb=reshape(@view(vec[dims.ivw.v_wb.index]), dims.ivw.v_wb.shape),
        v_wr=reshape(@view(vec[dims.ivw.v_wr.index]), dims.ivw.v_wr.shape),
        v_ww=reshape(@view(vec[dims.ivw.v_ww.index]), dims.ivw.v_ww.shape),
    )

    # - linear system - #
    linsys = (;
        A_bb=reshape(@view(vec[dims.linsys.A_bb.index]), dims.linsys.A_bb.shape),
        # A_bb_LU=LinearAlgebra.LU(
        #     reshape(@view(vec[dims.linsys.A_bb_LU.index]), dims.linsys.A_bb_LU.shape),
        #     [i for i in 1:dims.linsys.A_bb_LU.shape[1]],
        #     0,
        # ),
        b_bf=reshape(@view(vec[dims.linsys.b_bf.index]), dims.linsys.b_bf.shape),
        A_bw=reshape(@view(vec[dims.linsys.A_bw.index]), dims.linsys.A_bw.shape),
        A_pw=reshape(@view(vec[dims.linsys.A_pw.index]), dims.linsys.A_pw.shape),
        A_br=reshape(@view(vec[dims.linsys.A_br.index]), dims.linsys.A_br.shape),
        A_pr=reshape(@view(vec[dims.linsys.A_pr.index]), dims.linsys.A_pr.shape),
        # lu_decomp_flag=reshape(@view(vec[dims.linsys.lu_decomp_flag.index]), dims.linsys.lu_decomp_flag.shape),
    )

    # - blade element geometry - #
    blade_elements = (;
        B=reshape(@view(vec[dims.blade_elements.B.index]), dims.blade_elements.B.shape),
        Rtip=reshape(
            @view(vec[dims.blade_elements.Rtip.index]), dims.blade_elements.Rtip.shape
        ),
        Rhub=reshape(
            @view(vec[dims.blade_elements.Rhub.index]), dims.blade_elements.Rhub.shape
        ),
        fliplift=reshape(
            @view(vec[dims.blade_elements.fliplift.index]),
            dims.blade_elements.fliplift.shape,
        ),
        chords=reshape(
            @view(vec[dims.blade_elements.chords.index]), dims.blade_elements.chords.shape
        ),
        twists=reshape(
            @view(vec[dims.blade_elements.twists.index]), dims.blade_elements.twists.shape
        ),
        stagger=reshape(
            @view(vec[dims.blade_elements.stagger.index]), dims.blade_elements.stagger.shape
        ),
        solidity=reshape(
            @view(vec[dims.blade_elements.solidity.index]),
            dims.blade_elements.solidity.shape,
        ),
        rotor_panel_centers=reshape(
            @view(vec[dims.blade_elements.rotor_panel_centers.index]),
            dims.blade_elements.rotor_panel_centers.shape,
        ),
        inner_fraction=reshape(
            @view(vec[dims.blade_elements.inner_fraction.index]),
            dims.blade_elements.inner_fraction.shape,
        ),
    )

    wakeK = reshape(@view(vec[dims.wakeK.index]), dims.wakeK.shape)

    # # - index maps - #
    # idmaps = (;
    #     wake_node_ids_along_casing_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_casing_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_casing_wake_interface.shape,
    #     ),
    #     wake_node_ids_along_centerbody_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_centerbody_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_centerbody_wake_interface.shape,
    #     ),
    #     id_of_first_casing_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.shape,
    #     ),
    #     id_of_first_centerbody_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.shape,
    #     ),
    #     rotorwakenodeid=reshape(
    #         @view(vec[dims.idmaps.rotorwakenodeid.index]), dims.idmaps.rotorwakenodeid.shape
    #     ),
    #     wake_nodemap=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    #     wake_endnodeidxs=reshape(
    #         @view(vec[dims.idmaps.wake_endnodeidxs.index]),
    #         dims.idmaps.wake_endnodeidxs.shape,
    #     ),
    #     rotor_indices_in_wake=reshape(
    #         @view(vec[dims.idmaps.rotor_indices_in_wake.index]),
    #         dims.idmaps.rotor_indices_in_wake.shape,
    #     ),
    #     body_totnodes=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    # )

    return (; Gamr, sigr, gamw, operating_point, ivr, ivw, linsys, blade_elements, wakeK)
end

function withdraw_solve_parameter_cache(
    solver_options::TS, vec, dims
   ) where {TS<:ExternalSolverOptions}

    # - Initial Guesses - #
    vz_rotor = reshape(@view(vec[dims.vz_rotor.index]), dims.vz_rotor.shape)
    vtheta_rotor = reshape(@view(vec[dims.vtheta_rotor.index]), dims.vtheta_rotor.shape)
    Cm_wake = reshape(@view(vec[dims.Cm_wake.index]), dims.Cm_wake.shape)

    # - Operating Point - #
    operating_point = (;
        Vinf=reshape(
            @view(vec[dims.operating_point.Vinf.index]), dims.operating_point.Vinf.shape
        ),
        Minf=reshape(
            @view(vec[dims.operating_point.Minf.index]), dims.operating_point.Minf.shape
        ),
        rhoinf=reshape(
            @view(vec[dims.operating_point.rhoinf.index]), dims.operating_point.rhoinf.shape
        ),
        muinf=reshape(
            @view(vec[dims.operating_point.muinf.index]), dims.operating_point.muinf.shape
        ),
        asound=reshape(
            @view(vec[dims.operating_point.asound.index]), dims.operating_point.asound.shape
        ),
        Ptot=reshape(
            @view(vec[dims.operating_point.Ptot.index]), dims.operating_point.Ptot.shape
        ),
        Ttot=reshape(
            @view(vec[dims.operating_point.Ttot.index]), dims.operating_point.Ttot.shape
        ),
        Omega=reshape(
            @view(vec[dims.operating_point.Omega.index]), dims.operating_point.Omega.shape
        ),
    )

    # - induced velocities on rotor - #
    ivr = (;
        v_rb=reshape(@view(vec[dims.ivr.v_rb.index]), dims.ivr.v_rb.shape),
        v_rr=reshape(@view(vec[dims.ivr.v_rr.index]), dims.ivr.v_rr.shape),
        v_rw=reshape(@view(vec[dims.ivr.v_rw.index]), dims.ivr.v_rw.shape),
    )

    # - induced velocities on wake - #
    ivw = (;
        v_wb=reshape(@view(vec[dims.ivw.v_wb.index]), dims.ivw.v_wb.shape),
        v_wr=reshape(@view(vec[dims.ivw.v_wr.index]), dims.ivw.v_wr.shape),
        v_ww=reshape(@view(vec[dims.ivw.v_ww.index]), dims.ivw.v_ww.shape),
    )

    # - linear system - #
    linsys = (;
        A_bb=reshape(@view(vec[dims.linsys.A_bb.index]), dims.linsys.A_bb.shape),
        b_bf=reshape(@view(vec[dims.linsys.b_bf.index]), dims.linsys.b_bf.shape),
        A_bw=reshape(@view(vec[dims.linsys.A_bw.index]), dims.linsys.A_bw.shape),
        A_pw=reshape(@view(vec[dims.linsys.A_pw.index]), dims.linsys.A_pw.shape),
        A_br=reshape(@view(vec[dims.linsys.A_br.index]), dims.linsys.A_br.shape),
        A_pr=reshape(@view(vec[dims.linsys.A_pr.index]), dims.linsys.A_pr.shape),
    )

    # - blade element geometry - #
    blade_elements = (;
        B=reshape(@view(vec[dims.blade_elements.B.index]), dims.blade_elements.B.shape),
        Rtip=reshape(
            @view(vec[dims.blade_elements.Rtip.index]), dims.blade_elements.Rtip.shape
        ),
        Rhub=reshape(
            @view(vec[dims.blade_elements.Rhub.index]), dims.blade_elements.Rhub.shape
        ),
        fliplift=reshape(
            @view(vec[dims.blade_elements.fliplift.index]),
            dims.blade_elements.fliplift.shape,
        ),
        chords=reshape(
            @view(vec[dims.blade_elements.chords.index]), dims.blade_elements.chords.shape
        ),
        twists=reshape(
            @view(vec[dims.blade_elements.twists.index]), dims.blade_elements.twists.shape
        ),
        stagger=reshape(
            @view(vec[dims.blade_elements.stagger.index]), dims.blade_elements.stagger.shape
        ),
        solidity=reshape(
            @view(vec[dims.blade_elements.solidity.index]),
            dims.blade_elements.solidity.shape,
        ),
        rotor_panel_centers=reshape(
            @view(vec[dims.blade_elements.rotor_panel_centers.index]),
            dims.blade_elements.rotor_panel_centers.shape,
        ),
        inner_fraction=reshape(
            @view(vec[dims.blade_elements.inner_fraction.index]),
            dims.blade_elements.inner_fraction.shape,
        ),
    )

    wakeK = reshape(@view(vec[dims.wakeK.index]), dims.wakeK.shape)

    # # - index maps - #
    # idmaps = (;
    #     wake_node_ids_along_casing_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_casing_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_casing_wake_interface.shape,
    #     ),
    #     wake_node_ids_along_centerbody_wake_interface=reshape(
    #         @view(vec[dims.idmaps.wake_node_ids_along_centerbody_wake_interface.index]),
    #         dims.idmaps.wake_node_ids_along_centerbody_wake_interface.shape,
    #     ),
    #     id_of_first_casing_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_casing_panel_aft_of_each_rotor.shape,
    #     ),
    #     id_of_first_centerbody_panel_aft_of_each_rotor=reshape(
    #         @view(vec[dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.index]),
    #         dims.idmaps.id_of_first_centerbody_panel_aft_of_each_rotor.shape,
    #     ),
    #     rotorwakenodeid=reshape(
    #         @view(vec[dims.idmaps.rotorwakenodeid.index]), dims.idmaps.rotorwakenodeid.shape
    #     ),
    #     wake_nodemap=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    #     wake_endnodeidxs=reshape(
    #         @view(vec[dims.idmaps.wake_endnodeidxs.index]),
    #         dims.idmaps.wake_endnodeidxs.shape,
    #     ),
    #     rotor_indices_in_wake=reshape(
    #         @view(vec[dims.idmaps.rotor_indices_in_wake.index]),
    #         dims.idmaps.rotor_indices_in_wake.shape,
    #     ),
    #     body_totnodes=reshape(
    #         @view(vec[dims.idmaps.wake_nodemap.index]), dims.idmaps.wake_nodemap.shape
    #     ),
    # )

    return (;
        vz_rotor,
        vtheta_rotor,
        Cm_wake,
        operating_point,
        ivr,
        ivw,
        linsys,
        blade_elements,
        wakeK,
    )
end

function withdraw_solve_parameter_cache(solver_options::SIAMFANLEOptions, vec, dims)
    tuple = withdraw_solve_parameter_cache(NonlinearSolveOptions(), vec, dims)

    resid_cache_vec = reshape(
        @view(vec[dims.resid_cache_vec.index]), dims.resid_cache_vec.shape
    )
    krylov_cache_vec = reshape(
        @view(vec[dims.krylov_cache_vec.index]), dims.krylov_cache_vec.shape
    )
    jvp_cache_vec = reshape(@view(vec[dims.jvp_cache_vec.index]), dims.jvp_cache_vec.shape)

    return (; tuple..., resid_cache_vec, krylov_cache_vec, jvp_cache_vec)
end

"""
    withdraw_solve_container_cache(solver_options::SolverOptionsType, vec, dims)

Reshape the intermediate solve container cache vector using the saved dimensions tuple.

# Arguments
- `solver_options::SolverOptionsType` : Solver options type for dispatch.
- `vec::Vector{Float}` : vector cache of pre- and post-processing intermediate containers.
- `dims::NamedTuple` : Named tuple containing the indices and shape of the various items stored in the cache vector.

# Returns
- `solve_container_caching::NamedTuple` : Named tuple containing reshaped views of sections of the cache vector.
"""
function withdraw_solve_container_cache(solver_options::TS, vec, dims) where {TS<:InternalSolverOptions}
    return (;
        # Strengths
        gamb=reshape(@view(vec[dims.gamb.index]), dims.gamb.shape),
        rhs=reshape(@view(vec[dims.rhs.index]), dims.rhs.shape),

        # Blade Element Values
        Cz_rotor=reshape(@view(vec[dims.Cz_rotor.index]), dims.Cz_rotor.shape),
        Ctheta_rotor=reshape(@view(vec[dims.Ctheta_rotor.index]), dims.Ctheta_rotor.shape),
        Cmag_rotor=reshape(@view(vec[dims.Cmag_rotor.index]), dims.Cmag_rotor.shape),
        cl=reshape(@view(vec[dims.cl.index]), dims.cl.shape),
        cd=reshape(@view(vec[dims.cd.index]), dims.cd.shape),
        beta1=reshape(@view(vec[dims.beta1.index]), dims.beta1.shape),
        alpha=reshape(@view(vec[dims.alpha.index]), dims.alpha.shape),
        reynolds=reshape(@view(vec[dims.reynolds.index]), dims.reynolds.shape),
        mach=reshape(@view(vec[dims.mach.index]), dims.mach.shape),
        vz_rotor=reshape(@view(vec[dims.vz_rotor.index]), dims.vz_rotor.shape),
        vtheta_rotor=reshape(@view(vec[dims.vtheta_rotor.index]), dims.vtheta_rotor.shape),

        # Wake Velocities
        vz_wake=reshape(@view(vec[dims.vz_wake.index]), dims.vz_wake.shape),
        vr_wake=reshape(@view(vec[dims.vr_wake.index]), dims.vr_wake.shape),
        Cm_wake=reshape(@view(vec[dims.Cm_wake.index]), dims.Cm_wake.shape),
        Cm_avg=reshape(@view(vec[dims.Cm_avg.index]), dims.Cm_avg.shape),
        Gamma_tilde=reshape(@view(vec[dims.Gamma_tilde.index]), dims.Gamma_tilde.shape),
        H_tilde=reshape(@view(vec[dims.H_tilde.index]), dims.H_tilde.shape),
        deltaGamma2=reshape(@view(vec[dims.deltaGamma2.index]), dims.deltaGamma2.shape),
        deltaH=reshape(@view(vec[dims.deltaH.index]), dims.deltaH.shape),

        # State estimates
        Gamr_est=reshape(@view(vec[dims.Gamr_est.index]), dims.Gamr_est.shape),
        sigr_est=reshape(@view(vec[dims.sigr_est.index]), dims.sigr_est.shape),
        gamw_est=reshape(@view(vec[dims.gamw_est.index]), dims.gamw_est.shape),

        # Convergence items
        deltaG=reshape(@view(vec[dims.deltaG.index]), dims.deltaG.shape),
        deltaG_prev=reshape(@view(vec[dims.deltaG_prev.index]), dims.deltaG_prev.shape),
        deltag=reshape(@view(vec[dims.deltag.index]), dims.deltag.shape),
        deltag_prev=reshape(@view(vec[dims.deltag_prev.index]), dims.deltag_prev.shape),
        deltas=reshape(@view(vec[dims.deltas.index]), dims.deltas.shape),
        maxBGamr=reshape(@view(vec[dims.maxBGamr.index]), dims.maxBGamr.shape),
        maxdeltaBGamr=reshape(
            @view(vec[dims.maxdeltaBGamr.index]), dims.maxdeltaBGamr.shape
        ),
        maxdeltagamw=reshape(@view(vec[dims.maxdeltagamw.index]), dims.maxdeltagamw.shape),
    )
end

function withdraw_solve_container_cache(
    solver_options::TS, vec, dims
   ) where {TS<:ExternalSolverOptions}
    return (;
        # Strengths
        gamb=reshape(@view(vec[dims.gamb.index]), dims.gamb.shape),
        rhs=reshape(@view(vec[dims.rhs.index]), dims.rhs.shape),
        Gamr=reshape(@view(vec[dims.Gamr.index]), dims.Gamr.shape),
        sigr=reshape(@view(vec[dims.sigr.index]), dims.sigr.shape),
        gamw=reshape(@view(vec[dims.gamw.index]), dims.gamw.shape),

        # Blade Element Values
        Cz_rotor=reshape(@view(vec[dims.Cz_rotor.index]), dims.Cz_rotor.shape),
        Ctheta_rotor=reshape(@view(vec[dims.Ctheta_rotor.index]), dims.Ctheta_rotor.shape),
        Cmag_rotor=reshape(@view(vec[dims.Cmag_rotor.index]), dims.Cmag_rotor.shape),
        cl=reshape(@view(vec[dims.cl.index]), dims.cl.shape),
        cd=reshape(@view(vec[dims.cd.index]), dims.cd.shape),
        beta1=reshape(@view(vec[dims.beta1.index]), dims.beta1.shape),
        alpha=reshape(@view(vec[dims.alpha.index]), dims.alpha.shape),
        reynolds=reshape(@view(vec[dims.reynolds.index]), dims.reynolds.shape),
        mach=reshape(@view(vec[dims.mach.index]), dims.mach.shape),

        # Wake Velocities
        vz_wake=reshape(@view(vec[dims.vz_wake.index]), dims.vz_wake.shape),
        vr_wake=reshape(@view(vec[dims.vr_wake.index]), dims.vr_wake.shape),
        Cm_avg=reshape(@view(vec[dims.Cm_avg.index]), dims.Cm_avg.shape),
        Gamma_tilde=reshape(@view(vec[dims.Gamma_tilde.index]), dims.Gamma_tilde.shape),
        H_tilde=reshape(@view(vec[dims.H_tilde.index]), dims.H_tilde.shape),
        deltaGamma2=reshape(@view(vec[dims.deltaGamma2.index]), dims.deltaGamma2.shape),
        deltaH=reshape(@view(vec[dims.deltaH.index]), dims.deltaH.shape),

        # State estimates
        vz_est=reshape(@view(vec[dims.vz_est.index]), dims.vz_est.shape),
        vtheta_est=reshape(@view(vec[dims.vtheta_est.index]), dims.vtheta_est.shape),
        Cm_est=reshape(@view(vec[dims.Cm_est.index]), dims.Cm_est.shape),
    )
end
