"""
    allocate_body_panel_containers!(total_length, problem_dimensions::ProblemDimensions)

A helper function is assembling the prepost_container_cache.

# Arguments
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Returns
- `body_vortex_panels::NamedTuple` : A named tuple containing the dimensions needed to reshape the cache with regards to the body vortex panel object
"""
function allocate_body_panel_container!(total_length, problem_dimensions::ProblemDimensions)
    (;
        ndn,    # number of duct nodes
        ncbn,   # number of center_body nodes
    ) = problem_dimensions

    nn = [ndn, ncbn]
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container!(total_length, nn, np, tn, tp, nb)
end

"""
    allocate_rotor_panel_containers!(total_length, problem_dimensions::ProblemDimensions)

A helper function is assembling the prepost_container_cache.

# Arguments
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Returns
- `rotor_source_panels::NamedTuple` : A named containing the dimensions needed to reshape the cache with regards to the rotor source panel object
"""
function allocate_rotor_panel_container!(
    total_length, problem_dimensions::ProblemDimensions
)
    (;
        nrotor,     # number of rotors
        nws,    # number of wake sheets (also rotor nodes)
    ) = problem_dimensions

    nn = nws * ones(Int, nrotor)
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container!(total_length, nn, np, tn, tp, nb)
end

"""
    allocate_wake_panel_containers!(total_length, problem_dimensions::ProblemDimensions)

A helper function is assembling the prepost_container_cache.

# Arguments
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Returns
- `wake_vortex_panels::NamedTuple` : A named containing the dimensions needed to reshape the cache with regards to the wake vortex panel object
"""
function allocate_wake_panel_container!(total_length, problem_dimensions::ProblemDimensions)
    (;
        nws,    # number of wake sheets (also rotor nodes)
        nwsn,   # number of nodes in each wake sheet
    ) = problem_dimensions

    nn = nwsn * ones(Int, nws)
    # number of panels to generate for each body
    np = nn .- 1
    # number of bodies
    nb = length(np)
    # total number of panels in system
    tp = sum(np)
    # total number of nodes in system
    tn = tp + nb

    return allocate_panel_container!(total_length, nn, np, tn, tp, nb)
end

"""
    allocate_panel_container!(total_length, nn, np, tn, tp, nb)

A helper function is assembling the prepost_container_cache.

# Arguments
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.
- `nn::Int` : number of nodes in each body, rotor, or wake sheet
- `np::Int` : number of panels in each body, rotor, or wake sheet
- `tn::Int` : number of total nodes among the bodies, rotors, or wake sheets
- `tp::Int` : number of total panels among the bodies, rotors, or wake sheets
- `nb::Int` : number of bodies, rotors, or wake sheets

# Returns
- `panel::NamedTuple` : A named containing the dimensions needed to reshape the cache with regards to an arbitrary panel set
"""
function allocate_panel_container!(total_length, nn, np, tn, tp, nb)
    s = size(nn)
    l = lfs(s)
    nnode = cache_dims!(total_length, l, s)

    # number of panels to generate for each body
    npanel = cache_dims!(total_length, l, s)

    s = (1,)
    l = lfs(s)
    nbodies = cache_dims!(total_length, l, s)

    s = (1,)
    l = lfs(s)
    totpanel = cache_dims!(total_length, l, s)

    s = (1,)
    l = lfs(s)
    totnode = cache_dims!(total_length, l, s)

    s = (2, tn)
    l = lfs(s)
    node = cache_dims!(total_length, l, s)

    s = (2, tp)
    l = lfs(s)
    controlpoint = cache_dims!(total_length, l, s)

    normal = cache_dims!(total_length, l, s)

    tangent = cache_dims!(total_length, l, s)

    nodemap = cache_dims!(total_length, l, s)

    s = (tp)
    l = lfs(s)
    influence_length = cache_dims!(total_length, l, s)

    s = (nb, 2, 2)
    l = lfs(s)
    endnodes = cache_dims!(total_length, l, s)

    tenode = cache_dims!(total_length, l, s)

    s = (2, nb)
    l = lfs(s)

    itcontrolpoint = cache_dims!(total_length, l, s)

    itnormal = cache_dims!(total_length, l, s)

    ittangent = cache_dims!(total_length, l, s)

    tenormal = cache_dims!(total_length, l, s)

    tendotn = cache_dims!(total_length, l, s)

    tencrossn = cache_dims!(total_length, l, s)

    endnodeidxs = cache_dims!(total_length, l, s)

    endpanelidxs = cache_dims!(total_length, l, s)

    teadjnodeidxs = cache_dims!(total_length, l, s)

    s = (nb,)
    l = lfs(s)
    teinfluence_length = cache_dims!(total_length, l, s)

    prescribednodeidxs = cache_dims!(total_length, l, s)

    return (;
        nnode,
        npanel,
        nbodies,
        totpanel,
        totnode,
        node,
        controlpoint,
        normal,
        tangent,
        nodemap,
        influence_length,
        endnodes,
        tenode,
        itcontrolpoint,
        itnormal,
        ittangent,
        tenormal,
        tendotn,
        tencrossn,
        endnodeidxs,
        endpanelidxs,
        teadjnodeidxs,
        teinfluence_length,
        prescribednodeidxs,
    )
end

"""
    allocate_panel_containers!(total_length, problem_dimensions::ProblemDimensions)

A helper function is assembling the prepost_container_cache.

# Arguments
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Returns
- `panels::NamedTuple` : A named tuple of named tuples containing the dimensions needed to reshape the cache with regards to the panel objects
"""
function allocate_panel_containers!(total_length, problem_dimensions::ProblemDimensions)
    body_vortex_panels = allocate_body_panel_container!(total_length, problem_dimensions)
    rotor_source_panels = allocate_rotor_panel_container!(total_length, problem_dimensions)
    wake_vortex_panels = allocate_wake_panel_container!(total_length, problem_dimensions)

    return (; body_vortex_panels, rotor_source_panels, wake_vortex_panels)
end

"""
    allocate_prepost_container_cache(paneling_constants::PanelingConstants)
    allocate_prepost_container_cache(problem_dimensions::ProblemDimensions)

Allocate the pre- and post-processing cache (used for intermediate calculations) based on paneling constants or problem dimensions.

# Arguments
- `paneling_constants::PanelingConstants` : a PanelingConstants object
OR
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Keyword Arguments
- `fd_chunk_size::Int=12` : chunk size to use for PreallocationTools caches.  Note that the automated chunk size for DuctAPE will always be the ForwardDiff threshold of 12 due to the size of the system, so it will be best to leave this at the default unless further development allows for chunk size selection for individual solvers.
- `levels::Int=1` : levels for nested duals.  Note that since ImplicitAD is being used for all solves, there should be no need for more than 1 level.


# Returns
- `prepost_container_caching::NamedTuple` : a Named Tuple containing:
  - `prepost_container_cache::PreallocationTools.DiffCache` : the cache
  - `prepost_container_cache_dims::NamedTuple` : a named tuple containing the dimensions used for reshaping the cache when needed.
"""
function allocate_prepost_container_cache(
    paneling_constants::PanelingConstants; fd_chunk_size=12, levels=1
)
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_prepost_container_cache(
        problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

function allocate_prepost_body_containers!(total_length, nbp, ncp, ndn, ncbn, nbodies)
    s = (2, nbp)
    l = lfs(s)
    Vtot_in = cache_dims!(total_length, l, s)

    Vtot_out = cache_dims!(total_length, l, s)

    Vtot_prejump = cache_dims!(total_length, l, s)

    vtot_body = cache_dims!(total_length, l, s)

    vtot_jump = cache_dims!(total_length, l, s)

    vtot_wake = cache_dims!(total_length, l, s)

    vtot_rotors = cache_dims!(total_length, l, s)

    s = (nbp,)
    l = lfs(s)
    Vtan_in = cache_dims!(total_length, l, s)

    Vtan_out = cache_dims!(total_length, l, s)

    cp_in = cache_dims!(total_length, l, s)
    cp_out = cache_dims!(total_length, l, s)

    s = (ncp,)
    l = lfs(s)
    vtan_casing_in = cache_dims!(total_length, l, s)

    vtan_casing_out = cache_dims!(total_length, l, s)

    casing_zpts = cache_dims!(total_length, l, s)

    cp_casing_in = cache_dims!(total_length, l, s)
    cp_casing_out = cache_dims!(total_length, l, s)

    s = (ndn - 1 - ncp,)
    l = lfs(s)
    vtan_nacelle_in = cache_dims!(total_length, l, s)

    vtan_nacelle_out = cache_dims!(total_length, l, s)

    nacelle_zpts = cache_dims!(total_length, l, s)

    cp_nacelle_in = cache_dims!(total_length, l, s)
    cp_nacelle_out = cache_dims!(total_length, l, s)

    s = (ncbn - 1,)
    l = lfs(s)
    vtan_center_body_in = cache_dims!(total_length, l, s)

    vtan_center_body_out = cache_dims!(total_length, l, s)

    center_body_zpts = cache_dims!(total_length, l, s)

    cp_center_body_in = cache_dims!(total_length, l, s)

    cp_center_body_out = cache_dims!(total_length, l, s)

    s = (ndn - 1,)
    l = lfs(s)
    duct_jump = cache_dims!(total_length, l, s)

    s = (ncbn - 1,)
    l = lfs(s)
    center_body_jump = cache_dims!(total_length, l, s)

    s = (nbp,)
    l = lfs(s)
    body_jump_term = cache_dims!(total_length, l, s)

    s = (nbodies,)
    l = lfs(s)
    body_inviscid_thrust = cache_dims!(total_length, l, s)

    body_viscous_drag = cache_dims!(total_length, l, s)

    body_thrust = cache_dims!(total_length, l, s)

    body_force_coefficient = cache_dims!(total_length, l, s)

    return (; casing_zpts, nacelle_zpts, center_body_zpts), #zpts
    (;
        Vtot_in,
        Vtot_out,
        Vtan_in,
        Vtan_out,
        Vtot_prejump,
        vtot_body,
        duct_jump,
        center_body_jump,
        body_jump_term,
        vtot_jump,
        vtot_wake,
        vtot_rotors,
        vtan_casing_in,
        vtan_casing_out,
        vtan_nacelle_in,
        vtan_nacelle_out,
        vtan_center_body_in,
        vtan_center_body_out,
    ), # vtan_tuple
    (;
        cp_in,
        cp_out,
        cp_casing_in,
        cp_casing_out,
        cp_nacelle_in,
        cp_nacelle_out,
        cp_center_body_in,
        cp_center_body_out,
    ), # cp_tuple
    body_inviscid_thrust,
    body_viscous_drag,
    body_thrust,
    body_force_coefficient
end

function allocate_prepost_container_cache(
    problem_dimensions::ProblemDimensions; fd_chunk_size=12, levels=1
)
    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        ncp,    # number of casing panels
        ndn,    # number of duct nodes
        ncbn,   # number of center_body nodes
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nwsn,   # number of nodes in each wake sheet
        nbodies,
    ) = problem_dimensions

    # - initialize - #
    total_length = [0]

    ### --- PRE-PROCESSING --- ###

    # - COORDINATES - #
    s = (2, ndn)
    l = lfs(s)
    rp_duct_coordinates = cache_dims!(total_length, l, s)

    s = (2, ncbn)
    l = lfs(s)
    rp_center_body_coordinates = cache_dims!(total_length, l, s)

    s = (2, nwsn, nws)
    l = lfs(s)
    wake_grid = cache_dims!(total_length, l, s)

    s = (nrotor,)
    l = lfs(s)
    rotor_indices_in_wake = cache_dims!(total_length, l, s)

    # - PANELS - #
    panels = allocate_panel_containers!(total_length, problem_dimensions)

    # - INDUCED VELS - #
    s = (nbp, nbn, 2)
    l = lfs(s)
    v_bb = cache_dims!(total_length, l, s)

    s = (nbp, nrotor * nws, 2)
    l = lfs(s)
    v_br = cache_dims!(total_length, l, s)

    s = (nbp, nwn, 2)
    l = lfs(s)
    v_bw = cache_dims!(total_length, l, s)

    # - LINSYS - #
    s = (nbp, nbn)
    l = lfs(s)
    AICn = cache_dims!(total_length, l, s)

    s = (2, nbn)
    l = lfs(s)
    AICpcp = cache_dims!(total_length, l, s)

    s = (nbp,)
    l = lfs(s)
    vdnb = cache_dims!(total_length, l, s)

    s = (2,)
    l = lfs(s)
    vdnpcp = cache_dims!(total_length, l, s)

    ##### ----- POST PROCESSING ----- #####

    ### --- ROTOR Post-Processing Cache --- ###

    # TODO: consider moving rotor stuff to its own function like body stuff
    s = (nrotor,)
    l = lfs(s)
    rotor_inviscid_thrust = cache_dims!(total_length, l, s)

    rotor_viscous_thrust = cache_dims!(total_length, l, s)

    rotor_thrust = cache_dims!(total_length, l, s)

    rotor_inviscid_torque = cache_dims!(total_length, l, s)

    rotor_viscous_torque = cache_dims!(total_length, l, s)

    rotor_torque = cache_dims!(total_length, l, s)

    rotor_inviscid_power = cache_dims!(total_length, l, s)

    rotor_viscous_power = cache_dims!(total_length, l, s)

    rotor_power = cache_dims!(total_length, l, s)

    rotor_CT = cache_dims!(total_length, l, s)

    rotor_CQ = cache_dims!(total_length, l, s)

    rotor_CP = cache_dims!(total_length, l, s)

    rotor_efficiency = cache_dims!(total_length, l, s)

    induced_efficiency = cache_dims!(total_length, l, s)

    s = (nbe, nrotor)
    l = lfs(s)
    rotor_inviscid_thrust_dist = cache_dims!(total_length, l, s)

    rotor_viscous_thrust_dist = cache_dims!(total_length, l, s)

    rotor_inviscid_torque_dist = cache_dims!(total_length, l, s)

    rotor_viscous_torque_dist = cache_dims!(total_length, l, s)

    rotor_inviscid_power_dist = cache_dims!(total_length, l, s)

    rotor_viscous_power_dist = cache_dims!(total_length, l, s)

    blade_normal_force_per_unit_span = cache_dims!(total_length, l, s)

    blade_tangential_force_per_unit_span = cache_dims!(total_length, l, s)

    cn = cache_dims!(total_length, l, s)

    ct = cache_dims!(total_length, l, s)

    cphi = cache_dims!(total_length, l, s)

    sphi = cache_dims!(total_length, l, s)

    ### --- BODY Post-Processing Cache --- ###
    zpts, vtan_tuple, cp_tuple, body_inviscid_thrust, body_viscous_drag, body_thrust, body_force_coefficient = allocate_prepost_body_containers!(
        total_length, nbp, ncp, ndn, ncbn, nbodies
    )

    ### --- TOTALS Post-Processing Cache --- ###
    s = (1,)
    l = lfs(s)
    total_thrust = cache_dims!(total_length, l, s)

    total_torque = cache_dims!(total_length, l, s)

    total_power = cache_dims!(total_length, l, s)

    total_CT = cache_dims!(total_length, l, s)

    total_CQ = cache_dims!(total_length, l, s)

    total_CP = cache_dims!(total_length, l, s)

    total_efficiency = cache_dims!(total_length, l, s)

    ideal_efficiency = cache_dims!(total_length, l, s)

    # return tuple of initialized cache and associated dimensions
    return (;
        prepost_container_cache=PreallocationTools.DiffCache(
            zeros(total_length[]), fd_chunk_size; levels=levels
        ),
        prepost_container_cache_dims=(;
            ### --- PRE --- ###
            rp_duct_coordinates,
            rp_center_body_coordinates,
            wake_grid,
            rotor_indices_in_wake,
            panels,
            ivb=(; v_bb, v_br, v_bw),
            AICn,
            AICpcp,
            vdnb,
            vdnpcp,
            ### --- Post --- ###
            # - ROTOR - #
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
            blade_loading_intermediate_containers=(; cn, ct, cphi, sphi),
            # - BODY - #
            zpts,
            vtan_tuple,
            cp_tuple,
            body_inviscid_thrust,
            body_viscous_drag,
            body_thrust,
            body_force_coefficient,
            # - TOTALS - #
            total_thrust,
            total_torque,
            total_power,
            total_CT,
            total_CQ,
            total_CP,
            total_efficiency,
            ideal_efficiency,
        ),
    )
end

"""
    allocate_solve_parameter_cache(
        solve_type::SolverOptionsType,
        paneling_constants::PanelingConstants;
        airfoils,
        fd_chunk_size=12,
        levels=1,
    )
    allocate_solve_parameter_cache(
        solve_type::SolverOptionsType,
        problem_dimensions::ProblemDimensions;
        airfoils,
        fd_chunk_size=12,
        levels=1
    )

Allocate the solve parameter cache for parameters passed into the solver(s).

# Arguments
- `solve_type::SolverOptionsType` : Solver options type used for dispatch
- `paneling_constants::PanelingConstants` : a PanlingConstants object used for sizing
- `airfoils::Airfoil` : an array of airfoil composite types (MUST be a composite type), such as one of the types available in DuctAPE.C4Blade.
OR
- `solve_type::SolverOptionsType` : Solver options type used for dispatch
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object used for sizing
- `airfoils::Airfoil` : an array of airfoil composite types (MUST be a composite type), such as one of the types available in DuctAPE.C4Blade.

# Keyword Arguments
- `fd_chunk_size::Int=12` : chunk size to use for PreallocationTools caches.  Note that the automated chunk size for DuctAPE will always be the ForwardDiff threshold of 12 due to the size of the system, so it will be best to leave this at the default unless further development allows for chunk size selection for individual solvers.
- `levels::Int=1` : levels for nested duals.  Note that since ImplicitAD is being used for all solves, there should be no need for more than 1 level.

# Returns
- `solve_parameter_caching::NamedTuple` : a Named Tuple containing:
  - `solve_parameter_cache::PreallocationTools.DiffCache` : the cache
  - `solve_parameter_cache_dims::NamedTuple` : a named tuple containing the dimensions used for reshaping the cache when needed.
"""
function allocate_solve_parameter_cache(
    solve_type::TS,
    paneling_constants::PanelingConstants,
    airfoils;
    fd_chunk_size=12,
    levels=1,
) where {TS<:InternalSolverOptions}

    # - Get problem dimensions - #
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_parameter_cache(
        solve_type, problem_dimensions, airfoils; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

"""
    allocate_solve_parameter_extras!(
        solver_options::SolverOptionsType, input_length, total_length
    )

Includes additional caching for various solvers.  Currently only does anything for SIAMFANLEOptions types.

# Arguments
- `input_length::Int` : the number of state variables in the solver
- `total_length::Vector{Int}` : a one-element vector used to store the total length in order to know how large of a cache to allocate.  Is updated in place.

# Returns
- `solve_parameter_extras::NamedTuple` : A named tuple containing dimensions related to extra caching parameters used in various solvers.
"""
function allocate_solve_parameter_extras!(
    solver_options::SIAMFANLEOptions, input_length, total_length
)
    s = (input_length,)
    l = lfs(s)
    resid_cache_vec = cache_dims!(total_length, l, s)

    s = (input_length,)
    l = lfs(s)
    jvp_cache_vec = cache_dims!(total_length, l, s)

    s = (input_length, max(2, solver_options.linear_iteration_limit + 1))
    l = lfs(s)
    krylov_cache_vec = cache_dims!(total_length, l, s)

    return (; resid_cache_vec, krylov_cache_vec, jvp_cache_vec)
end

function allocate_solve_parameter_extras!(
    solver_options::SolverOptionsType, input_length, total_length
)
    return (;)
end

"""
    allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

Allocate caches for airfoil struct fields.


# Arguments:
- `airfoils::Vector{Vector{Airfoil}}` : Airfoil objects, where each element of the outer vector is associated with a rotor and the elements of the inner vectors are the input airfoils to be interpolated across each rotor blade. Note that even if only one rotor is being used, the input must still be a vector of vectors.
- `total_length::Vector{Float}` : updated in-place, used for sizing the overall PreallocationTools cache.
- `number_of_blade_elements::Int` : number of blade elements for analysis (this is defined from the `num_wake_sheets` input in the PanelingConstants input).
- `number_of_rotors::Int` : the number of rotors (note that a stator is considered a rotor).

# Returns:
- `airfoil_cache::NamedTuple` : contains:
  - `airfoil_cache_dims::Matrix{NamedTuple}` : A matrix (num blade elem x num rotor) of named tuples containing the indices and shapes required for extracting the airfoils from the PreallocationTools cache.
  - `airfoil_constructors::Matrix{type.name.wrapper}` : Wrappers for re-constructing the airfoil objects from the cache.

Notes:
- Due to how DuctAPE interpolates the blade elements, we require all airfoil types and the sizes of the fields in those types to be identical across blades; but they may be different for each rotor. For example, you could use a CCBlade-like airfoil for the rotor and a DFDC-like airfoil for the stator, but all the rotor airfoils must be the same type with fields of the same size and likewise for the stator.
- Since the airfoil input structure is a vector of vectors, the vectors for each rotor may be different sizes as well. For example, you can define 6 airfoils to get interpolated for the rotor, and 2 for the stator (or however many you want).
"""
function allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

    airfoil_cache_dims = reshape(
        reduce(
            vcat,
            [
                [
                    ntt.namedtuple(
                        propertynames(airfoils[irotor][1]),
                        [
                            cache_dims!(
                                total_length,
                                length(getproperty(airfoils[irotor][1], p)),
                                size(getproperty(airfoils[irotor][1], p)),
                            ) for p in propertynames(airfoils[irotor][1])
                        ],
                    ) for _ in 1:nbe
                ] for irotor in 1:nrotor
            ],
        ),
        (nbe, nrotor),
    )

    airfoil_constructors = reshape(
        reduce(
            vcat,
            [
                [typeof(airfoils[irotor][1]).name.wrapper for _ in 1:nbe] for
                irotor in 1:nrotor
            ],
        ),
        (nbe, nrotor),
    )

    return (; airfoil_cache_dims, airfoil_constructors)
end

function allocate_solve_parameter_cache(
    solve_type::TS,
    problem_dimensions::ProblemDimensions,
    airfoils;
    fd_chunk_size=12,
    levels=1,
) where {TS<:InternalSolverOptions}
    (;
        nrotor,     # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbp,    # number of body panels
        nws,    # number of wake sheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
    ) = problem_dimensions

    # - initialize - #
    total_length = [0]

    # - Initial Guesses - #
    s = (nbe, nrotor)
    l = lfs(s)
    Gamr = cache_dims!(total_length, l, s)

    s = (nws, nrotor)
    l = lfs(s)
    sigr = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    gamw = cache_dims!(total_length, l, s)

    # save state dimensions
    state_dims = (; Gamr, sigr, gamw)

    # - Operating Point - #
    s = (1,)
    l = lfs(s)
    Vinf = cache_dims!(total_length, l, s)

    Minf = cache_dims!(total_length, l, s)

    rhoinf = cache_dims!(total_length, l, s)

    muinf = cache_dims!(total_length, l, s)

    asound = cache_dims!(total_length, l, s)

    Ptot = cache_dims!(total_length, l, s)

    Ttot = cache_dims!(total_length, l, s)

    s = (nrotor,)
    l = lfs(s)
    Omega = cache_dims!(total_length, l, s)

    # - Induced Velocities on Rotors - #
    s = (nrotor * nbe, nbn, 2)
    l = lfs(s)
    v_rb = cache_dims!(total_length, l, s)

    s = (nrotor * nbe, nrotor * nws, 2)
    l = lfs(s)
    v_rr = cache_dims!(total_length, l, s)

    s = (nrotor * nbe, nwn, 2)
    l = lfs(s)
    v_rw = cache_dims!(total_length, l, s)

    # - Induced Velocities on Wakes - #
    s = (nwp, nbn, 2)
    l = lfs(s)
    v_wb = cache_dims!(total_length, l, s)

    s = (nwp, nrotor * nws, 2)
    l = lfs(s)
    v_wr = cache_dims!(total_length, l, s)

    s = (nwp, nwn, 2)
    l = lfs(s)
    v_ww = cache_dims!(total_length, l, s)

    # - Linear System - #

    s = (nbn + 2, nbn + 2)
    l = lfs(s)
    A_bb = cache_dims!(total_length, l, s)

    s = (nbn + 2,)
    l = lfs(s)
    b_bf = cache_dims!(total_length, l, s)

    s = (nbp, nwn)
    l = lfs(s)
    A_bw = cache_dims!(total_length, l, s)

    s = (2, nwn)
    l = lfs(s)
    A_pw = cache_dims!(total_length, l, s)

    s = (nbp, nrotor * nws)
    l = lfs(s)
    A_br = cache_dims!(total_length, l, s)

    s = (2, nrotor * nws)
    l = lfs(s)
    A_pr = cache_dims!(total_length, l, s)

    # - Blade Elements - #
    s = (nrotor,)
    l = lfs(s)
    B = cache_dims!(total_length, l, s)

    Rtip = cache_dims!(total_length, l, s)

    Rhub = cache_dims!(total_length, l, s)

    is_stator = cache_dims!(total_length, l, s)

    s = (nbe, nrotor)
    l = lfs(s)
    chords = cache_dims!(total_length, l, s)

    twists = cache_dims!(total_length, l, s)

    stagger = cache_dims!(total_length, l, s)

    solidity = cache_dims!(total_length, l, s)

    rotor_panel_centers = cache_dims!(total_length, l, s)

    inner_fraction = cache_dims!(total_length, l, s)

    inner_airfoil = allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

    outer_airfoil = allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

    # - Wake constants - #
    s = (nwn,)
    l = lfs(s)
    wakeK = cache_dims!(total_length, l, s)

    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(
            zeros(total_length[]), fd_chunk_size; levels=levels
        ),
        solve_parameter_cache_dims=(;
            state_dims,
            Gamr,
            sigr,
            gamw,
            operating_point=(; Vinf, Minf, rhoinf, muinf, asound, Ptot, Ttot, Omega),
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(; v_wb, v_wr, v_ww),
            linsys=(; A_bb, b_bf, A_bw, A_pw, A_br, A_pr),
            blade_elements=(;
                B,
                is_stator,
                chords,
                twists,
                stagger,
                solidity,
                rotor_panel_centers,
                inner_fraction,
                Rtip,
                Rhub,
                inner_airfoil,
                outer_airfoil,
            ),
            wakeK,
        ),
    )
end

function allocate_solve_parameter_cache(
    solve_type::TS,
    paneling_constants::PanelingConstants,
    airfoils;
    fd_chunk_size=12,
    levels=1,
) where {TS<:ExternalSolverOptions}

    # - Get problem dimensions - #
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_parameter_cache(
        solve_type, problem_dimensions, airfoils; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

function allocate_solve_parameter_cache(
    solve_type::TS,
    problem_dimensions::ProblemDimensions,
    airfoils;
    fd_chunk_size=12,
    levels=1,
) where {TS<:ExternalSolverOptions}
    (;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbp,    # number of body paneallocate_solve_parameter_extrasheets (also rotor nodes)
        nbe,    # number of blade elements (also rotor panels)
        nws, #number of wake sheets
    ) = problem_dimensions

    # - initialize - #
    total_length = [0]

    # - Initial Guesses - #
    s = (nbe, nrotor)
    l = lfs(s)
    vz_rotor = cache_dims!(total_length, l, s)

    vtheta_rotor = cache_dims!(total_length, l, s)

    s = (nwp,)
    l = lfs(s)
    Cm_wake = cache_dims!(total_length, l, s)

    SIAMFANLE_cache_vecs = allocate_solve_parameter_extras!(
        solve_type, total_length, total_length
    )

    # save state dimensions
    state_dims = (; vz_rotor, vtheta_rotor, Cm_wake)

    # - Operating Point - #
    s = (1,)
    l = lfs(s)
    Vinf = cache_dims!(total_length, l, s)

    Minf = cache_dims!(total_length, l, s)

    rhoinf = cache_dims!(total_length, l, s)

    muinf = cache_dims!(total_length, l, s)

    asound = cache_dims!(total_length, l, s)

    Ptot = cache_dims!(total_length, l, s)

    Ttot = cache_dims!(total_length, l, s)

    s = (nrotor,)
    l = lfs(s)
    Omega = cache_dims!(total_length, l, s)

    # - Induced Velocities on Rotors - #
    s = (nrotor * nbe, nbn, 2)
    l = lfs(s)
    v_rb = cache_dims!(total_length, l, s)

    s = (nrotor * nbe, nrotor * nws, 2)
    l = lfs(s)
    v_rr = cache_dims!(total_length, l, s)

    s = (nrotor * nbe, nwn, 2)
    l = lfs(s)
    v_rw = cache_dims!(total_length, l, s)

    # - Induced Velocities on Wakes - #
    s = (nwp, nbn, 2)
    l = lfs(s)
    v_wb = cache_dims!(total_length, l, s)

    s = (nwp, nrotor * nws, 2)
    l = lfs(s)
    v_wr = cache_dims!(total_length, l, s)

    s = (nwp, nwn, 2)
    l = lfs(s)
    v_ww = cache_dims!(total_length, l, s)

    # - Linear System - #

    s = (nbn + 2, nbn + 2)
    l = lfs(s)
    A_bb = cache_dims!(total_length, l, s)

    s = (nbn + 2,)
    l = lfs(s)
    b_bf = cache_dims!(total_length, l, s)

    s = (nbp, nwn)
    l = lfs(s)
    A_bw = cache_dims!(total_length, l, s)

    s = (2, nwn)
    l = lfs(s)
    A_pw = cache_dims!(total_length, l, s)

    s = (nbp, nrotor * nws)
    l = lfs(s)
    A_br = cache_dims!(total_length, l, s)

    s = (2, nrotor * nws)
    l = lfs(s)
    A_pr = cache_dims!(total_length, l, s)

    # - Blade Elements - #
    s = (nrotor,)
    l = lfs(s)
    B = cache_dims!(total_length, l, s)

    Rtip = cache_dims!(total_length, l, s)

    Rhub = cache_dims!(total_length, l, s)

    is_stator = cache_dims!(total_length, l, s)

    s = (nbe, nrotor)
    l = lfs(s)
    chords = cache_dims!(total_length, l, s)

    twists = cache_dims!(total_length, l, s)

    stagger = cache_dims!(total_length, l, s)

    solidity = cache_dims!(total_length, l, s)

    rotor_panel_centers = cache_dims!(total_length, l, s)

    inner_fraction = cache_dims!(total_length, l, s)

    inner_airfoil = allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

    outer_airfoil = allocate_airfoil_cache(airfoils, total_length, nbe, nrotor)

    # - Wake constants - #
    s = (nwn,)
    l = lfs(s)
    wakeK = cache_dims!(total_length, l, s)

    return (;
        solve_parameter_cache=PreallocationTools.DiffCache(
            zeros(total_length[]), fd_chunk_size; levels=levels
        ),
        solve_parameter_cache_dims=(;
            state_dims,
            vz_rotor,
            vtheta_rotor,
            Cm_wake,
            operating_point=(; Vinf, Minf, rhoinf, muinf, asound, Ptot, Ttot, Omega),
            ivr=(; v_rb, v_rr, v_rw),
            ivw=(; v_wb, v_wr, v_ww),
            linsys=(; A_bb, b_bf, A_bw, A_pw, A_br, A_pr),
            blade_elements=(;
                B,
                is_stator,
                chords,
                twists,
                stagger,
                solidity,
                rotor_panel_centers,
                inner_fraction,
                Rtip,
                Rhub,
                inner_airfoil,
                outer_airfoil,
            ),
            wakeK,
            SIAMFANLE_cache_vecs...,
        ),
    )
end

"""
    allocate_solve_container_cache(
        solve_type::SolverOptionsType,
        paneling_constants::PanelingConstants;
        fd_chunk_size=12,
        levels=1,
    )
    allocate_solve_container_cache(
        solve_type::SolverOptionsType,
        problem_dimensions::ProblemDimensions;
        fd_chunk_size=12,
        levels=1,
    )

Allocate the solve cache (used for intermediate calculations) based on paneling constants or problem dimensions.

# Arguments
- `paneling_constants::PanelingConstants` : a PanelingConstants object
OR
- `problem_dimensions::ProblemDimensions` : a ProblemDimensions object

# Keyword Arguments
- `fd_chunk_size::Int=12` : chunk size to use for PreallocationTools caches.  Note that the automated chunk size for DuctAPE will always be the ForwardDiff threshold of 12 due to the size of the system, so it will be best to leave this at the default unless further development allows for chunk size selection for individual solvers.
- `levels::Int=1` : levels for nested duals.  Note that since ImplicitAD is being used for all solves, there should be no need for more than 1 level.

# Returns
- `solve_container_caching::NamedTuple` : a Named Tuple containing:
  - `solve_container_cache::PreallocationTools.DiffCache` : the cache
  - `solve_container_cache_dims::NamedTuple` : a named tuple containing the dimensions used for reshaping the cache when needed.
"""
function allocate_solve_container_cache(
    solve_type::InternalSolverOptions,
    paneling_constants::PanelingConstants;
    fd_chunk_size=12,
    levels=1,
)
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_container_cache(
        solve_type, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

function allocate_solve_container_cache(
    solve_type::InternalSolverOptions,
    problem_dimensions::ProblemDimensions;
    fd_chunk_size=12,
    levels=1,
)
    (;
        nrotor,  # number of rotors
        nbodies, # number of bodies
        nwn,     # number of wake nodes
        nwp,     # number of wake panels
        nbn,     # number of body nodes
        nbe,     # number of blade elements (also rotor panels)
        nws,     # number of wake sheets (also rotor panel edges)
    ) = problem_dimensions

    # - initialize - #
    total_length = [0]

    # Strengths
    s = (nbn + nbodies,)
    l = lfs(s)
    gamb = cache_dims!(total_length, l, s)
    rhs = cache_dims!(total_length, l, s)

    # Blade Element Values
    s = (nbe, nrotor)
    l = lfs(s)
    beta1 = cache_dims!(total_length, l, s)

    Cz_rotor = cache_dims!(total_length, l, s)

    Ctheta_rotor = cache_dims!(total_length, l, s)

    Cmag_rotor = cache_dims!(total_length, l, s)

    cl = cache_dims!(total_length, l, s)

    cd = cache_dims!(total_length, l, s)

    alpha = cache_dims!(total_length, l, s)

    reynolds = cache_dims!(total_length, l, s)

    mach = cache_dims!(total_length, l, s)

    vz_rotor = cache_dims!(total_length, l, s)

    vtheta_rotor = cache_dims!(total_length, l, s)

    # Circulation
    Gamma_tilde = cache_dims!(total_length, l, s)

    H_tilde = cache_dims!(total_length, l, s)

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltaGamma2 = cache_dims!(total_length, l, s)

    deltaH = cache_dims!(total_length, l, s)

    # Wake Velocities
    s = (nwp,)
    l = lfs(s)
    vz_wake = cache_dims!(total_length, l, s)

    vr_wake = cache_dims!(total_length, l, s)

    # State estimates
    s = (nbe, nrotor)
    l = lfs(s)
    Gamr_est = cache_dims!(total_length, l, s)

    s = (nws, nrotor)
    l = lfs(s)
    sigr_est = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    gamw_est = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    Cm_avg = cache_dims!(total_length, l, s)

    s = (nwp,)
    l = lfs(s)
    Cm_wake = cache_dims!(total_length, l, s)

    s = (nbe, nrotor)
    l = lfs(s)
    deltaG = cache_dims!(total_length, l, s)

    deltaG_prev = cache_dims!(total_length, l, s)

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltas = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    deltag = cache_dims!(total_length, l, s)

    deltag_prev = cache_dims!(total_length, l, s)

    # Convergence criteria
    s = (nrotor,)
    l = lfs(s)
    maxBGamr = cache_dims!(total_length, l, s)

    maxdeltaBGamr = cache_dims!(total_length, l, s)

    s = (1,)
    l = lfs(s)
    maxdeltagamw = cache_dims!(total_length, l, s)

    return (;
        solve_container_cache=PreallocationTools.DiffCache(
            zeros(total_length[]), fd_chunk_size; levels=levels
        ),
        solve_container_cache_dims=(;
            gamb,
            rhs,
            Gamr_est,
            sigr_est,
            gamw_est,
            beta1,
            alpha,
            reynolds,
            mach,
            vz_rotor,
            vtheta_rotor,
            Cz_rotor,
            Ctheta_rotor,
            Cmag_rotor,
            cl,
            cd,
            vz_wake,
            vr_wake,
            Cm_wake,
            Cm_avg,
            Gamma_tilde,
            H_tilde,
            deltaGamma2,
            deltaH,
            deltaG,
            deltaG_prev,
            deltag,
            deltag_prev,
            deltas,
            maxBGamr,
            maxdeltaBGamr,
            maxdeltagamw,
        ),
    )
end

function allocate_solve_container_cache(
    solve_type::TS, paneling_constants::PanelingConstants; fd_chunk_size=12, levels=1
) where {TS<:ExternalSolverOptions}
    problem_dimensions = get_problem_dimensions(paneling_constants)

    return allocate_solve_container_cache(
        solve_type, problem_dimensions; fd_chunk_size=fd_chunk_size, levels=levels
    )
end

function allocate_solve_container_cache(
    solve_type::TS, problem_dimensions::ProblemDimensions; fd_chunk_size=12, levels=1
) where {TS<:ExternalSolverOptions}
    (;
        nrotor, # number of rotors
        nwn,    # number of wake nodes
        nwp,    # number of wake panels
        nbn,    # number of body nodes
        nbe,    # number of blade elements (also rotor panels)
    ) = problem_dimensions

    # - initialize - #
    total_length = [0]

    # Strengths
    # TODO: is this always going to be 2?, may want to add nb for nbodies to problem dims and then do +nb
    s = (nbn + 2,)
    l = lfs(s)
    gamb = cache_dims!(total_length, l, s)
    rhs = cache_dims!(total_length, l, s)

    s = (nbe, nrotor)
    l = lfs(s)
    Gamr = cache_dims!(total_length, l, s)

    s = (nbe + 1, nrotor)
    l = lfs(s)
    sigr = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    gamw = cache_dims!(total_length, l, s)

    # Blade Element Values
    s = (nbe, nrotor)
    l = lfs(s)
    beta1 = cache_dims!(total_length, l, s)

    Cz_rotor = cache_dims!(total_length, l, s)

    Ctheta_rotor = cache_dims!(total_length, l, s)

    Cmag_rotor = cache_dims!(total_length, l, s)

    cl = cache_dims!(total_length, l, s)

    cd = cache_dims!(total_length, l, s)

    alpha = cache_dims!(total_length, l, s)

    reynolds = cache_dims!(total_length, l, s)

    mach = cache_dims!(total_length, l, s)

    # Circulation
    Gamma_tilde = cache_dims!(total_length, l, s)

    H_tilde = cache_dims!(total_length, l, s)

    s = (nbe + 1, nrotor)
    l = lfs(s)
    deltaGamma2 = cache_dims!(total_length, l, s)

    deltaH = cache_dims!(total_length, l, s)

    # Wake Velocities
    s = (nwp,)
    l = lfs(s)
    vz_wake = cache_dims!(total_length, l, s)

    vr_wake = cache_dims!(total_length, l, s)

    # State estimates
    s = (nbe, nrotor)
    l = lfs(s)
    vz_est = cache_dims!(total_length, l, s)

    vtheta_est = cache_dims!(total_length, l, s)

    s = (nwp,)
    l = lfs(s)
    Cm_est = cache_dims!(total_length, l, s)

    s = (nwn,)
    l = lfs(s)
    Cm_avg = cache_dims!(total_length, l, s)

    return (;
        solve_container_cache=PreallocationTools.DiffCache(
            zeros(total_length[]), fd_chunk_size; levels=levels
        ),
        solve_container_cache_dims=(;
            gamb,
            rhs,
            Gamr,
            sigr,
            gamw,
            beta1,
            alpha,
            reynolds,
            mach,
            Cz_rotor,
            Ctheta_rotor,
            Cmag_rotor,
            cl,
            cd,
            vz_wake,
            vr_wake,
            vz_est,
            vtheta_est,
            Cm_est,
            Cm_avg,
            Gamma_tilde,
            H_tilde,
            deltaGamma2,
            deltaH,
        ),
    )
end
