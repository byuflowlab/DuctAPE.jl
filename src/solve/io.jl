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
    return [
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
        reduce(vcat, rsp.fliplift) # TODO: change this from a bool to 0.0 or 1.0 and update the coefficient calculation function for it.
        #TODO: need to remove nwake_sheets from rotor parameters, it's already in the paneling constants
        reduce(vcat, [paneling_constants...])
        reduce(vcat, [freestream...])
        reduce(vcat, [reference_parameters...])
        c4b.vectorize_airfoils(rsp.airfoils) #TODO add official test for this function
    ]
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
