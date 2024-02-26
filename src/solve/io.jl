"""
"""
function vectorize!(propulsor)

    # - Extract propulsor - #
    (;
        duct_coordinates, # Matrix
        centerbody_coordinates, # Matrix
        rotorstator_parameters, # Vector of NamedTuples of a bunch of stuff...
        paneling_constants, # NamedTuple of numbers and vectors of numbers
        freestream, # NamedTuple of numbers
        reference_parameters, # NamedTuple of numbers
    ) = propulsor

    # Rename rotorstator_parameters for convenience
    rsp = rotorstator_parameters

    # - Assemble (almost) everything into one big vector - #
    # Note: we use reduce for the vectors and namedtuples that could contain vectors
    xv = [
        duct_coordinates[:]
        centerbody_coordinates[:]
        reduce(vcat, rsp.B)
        reduce(vcat, rsp.Omega)
        reduce(vcat, rsp.r)
        reduce(vcat, rsp.chords)
        reduce(vcat, rsp.twists)
        reduce(vcat, rsp.Rtip)
        reduce(vcat, rsp.Rhub)
        reduce(vcat, rsp.tip_gap)
        reduce(vcat, rsp.rotor_z_loc)
        #TODO: need to remove nwake_sheets from rotor parameters, it's already in the paneling constants
        reduce(vcat, [paneling_constants...])
        reduce(vcat, [freestream...])
        reduce(vcat, [reference_parameters...])
    ]

    # - not quite ready for everything to be a differentiated variable so some things won't go into the vector - #
    pv = (; rsp.fliplift, rsp.airfoils)

    return xv, pv
end

"""
"""
function extract_state_vars(vars, dims)

    # - Separate out - #
    Vz_rotor = @views reshape(vars[dims.Vz.index], dims.Vz.shape)
    Vtheta_rotor = @views reshape(vars[dims.Vtheta.index], dims.Vtheta.shape)
    Cm_wake = @views reshape(vars[dims.Cm.index], dims.Cm.shape)

    return Vz_rotor, Vtheta_rotor, Cm_wake
end
