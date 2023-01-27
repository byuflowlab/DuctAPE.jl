#=

Functions regarding rotor aerodynamics

Authors: Judd Mehr,

=#

"""
"""
function calculate_induced_velocities()

    #TODO: vm comes from solving linear system

    #vtheta comes from gamma and such
    vtheta = 1.0 / (2.0 * pi * radial_positions) * (Gamma_tilde + 0.5 * B * Gamma)

    return nothing
end

"""
"""
function calculate_angle_of_attack(twist, Wm, Wtheta)

    # - Calculate Inflow Angle - #
    inflow = atan(Wm, -Wtheta)

    # - Calculate Angle of Attack - #
    alpha = twist - inflow

    return alpha
end

"""
"""
function calculate_inflow_velocities(Vinf, Omega, radial_position, vm, vtheta)
    Wm = Vinf + vm
    Wtheta = vtheta - Omega * radial_position
    Wmag = sqrt.(Wm .^ 2 .+ Wtheta .^ 2)
    return (Wm=Wm, Wtheta=Wtheta, Wmag=Wmag)
end

"""
"""
function search_polars(airfoil, alpha)
    return ccb.afeval(airfoil, alpha, 0.0, 0.0) #requires entries for Re and Ma, even though they aren't used.
end

"""
"""
function calculate_gamma_sigma(Vinf, blade_elements, vm=0.0, vtheta=0.0)
    TF = eltype(blade_elements.chords)

    # - Calculate Inflow Velocity - #
    if vm == 0.0 && vtheta == 0.0
        #first run through, no induced velocities

        inflow = calculate_inflow_velocities.(
            Ref(Vinf),
            Ref(blade_elements.omega),
            blade_elements.radial_positions,
            Ref(vm),
            Ref(vtheta),
        )[1]
    else
        inflow =
            calculate_inflow_velocities.(
                Ref(Vinf), Ref(blade_elements.omega), blade_elements.radial_positions
            )
    end

    # - Calculate Angle of Attack - #
    alphas = calculate_angle_of_attack.(blade_elements.twists, inflow.Wm, inflow.Wtheta)

    # - Look up lfit and drag data - #
    cl = zeros(TF, length(alphas))
    cd = zeros(TF, length(alphas))

    for i in 1:length(alphas)
        cl[i], cd[i] = search_polars(blade_elements.airfoils[i], alphas[i])
    end

    # - Calculate Gamma and Sigma - #
    Gamma = @.(0.5 * inflow.Wmag * blade_elements.chords * cl)
    Sigma = @.(0.5 * inflow.Wmag * blade_elements.chords * cd)

    return Gamma, Sigma
end

