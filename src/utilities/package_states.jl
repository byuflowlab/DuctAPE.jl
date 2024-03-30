"""
"""
function vectorize_velocity_states(vz_rotor, vtheta_rotor, Cm_wake)
    total_length = 0

    # vz_rotor Dims
    s = size(vz_rotor)
    lvz = lfs(s)
    Vz = (; index=(total_length + 1):(total_length + lvz), shape=s)
    total_length += lvz

    # vtheta_rotor Dims
    s = size(vtheta_rotor)
    lvt = lfs(s)
    Vtheta = (; index=(total_length + 1):(total_length + lvt), shape=s)
    total_length += lvt

    # Cm_wake Dims
    s = size(Cm_wake)
    l = lfs(s)
    Cm = (; index=(total_length + 1):(total_length + l), shape=s)

    return [reshape(vz_rotor, lvz); reshape(vtheta_rotor, lvt); Cm_wake],
    (; vz_rotor=Vz, vtheta_rotor=Vtheta, Cm_wake=Cm)
end

"""
"""
function vectorize_strength_states(Gamr, sigr, gamw)
    total_length = 0

    # Gamr Dims
    s = size(Gamr)
    lG = lfs(s)
    Γr = (; index=(total_length + 1):(total_length + lG), shape=s)
    total_length += lG

    # sigr Dims
    s = size(sigr)
    ls = lfs(s)
    σr = (; index=(total_length + 1):(total_length + ls), shape=s)
    total_length += ls

    # gamw Dims
    s = size(gamw)
    l = lfs(s)
    γw = (; index=(total_length + 1):(total_length + l), shape=s)

    return [reshape(Gamr, lG); reshape(sigr, ls); gamw], (; Gamr=Γr, sigr=σr, gamw=γw)
end

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
function extract_state_variables(solve_options::SolverOptions, vars, dims)

    # - Separate out - #
    vz_rotor = @views reshape(vars[dims.vz_rotor.index], dims.vz_rotor.shape)
    vtheta_rotor = @views reshape(vars[dims.vtheta_rotor.index], dims.vtheta_rotor.shape)
    Cm_wake = @views reshape(vars[dims.Cm_wake.index], dims.Cm_wake.shape)

    return vz_rotor, vtheta_rotor, Cm_wake
end

"""
"""
function extract_state_variables(solve_options::CSORSolverOptions, vars, dims)

    # - Separate out - #
    Gamr = @views reshape(vars[dims.Gamr.index], dims.Gamr.shape)
    sigr = @views reshape(vars[dims.sigr.index], dims.sigr.shape)
    gamw = @views reshape(vars[dims.gamw.index], dims.gamw.shape)

    return Gamr, sigr, gamw
end
