#=

Functions regarding rotor aerodynamics

Authors: Judd Mehr,

=#

"""
DONE. clean up and check
"""
function calculate_induced_velocities(
    BGamma,
    Gamma_tilde,
    blade_elements,
    A_bodies_to_rotor,
    gamma_bodies,
    A_wake_to_rotor,
    gamma_wake,
    # A_rotor_to_rotor,
    # Sigma,
)

    # - Rename for Convenience - #
    nr = length(blade_elements)

    # - Initialize - #
    vm = similar(BGamma)
    vtheta = similar(BGamma)

    for i in 1:nr
        # - Add Body Induced Velocities - #
        vm[:, i] = A_bodies_to_rotor[i] * gamma_bodies

        # - Add Wake Induced Velocities - #
        for w in 1:length(gamma_wake[:, 1])
            vm[:, i] += A_wake_to_rotor[w, i] * gamma_wake[i, :]'
        end

        # # - Add Rotor Induced Velocities - #
        # for j in 1:nr
        #     vm[:, i] += A_rotor_to_rotor[i, j] * Sigma[i]
        # end

        #vtheta comes from gamma and such
        vtheta[:, i] =
            1.0 / (2.0 * pi * blade_elements[i].radial_positions) *
            (Gamma_tilde[:, i] + 0.5 * BGamma[:, i])
    end

    return (vm=vm, vtheta=vtheta)
end

"""
DONE. clean up and check
"""
function calculate_angle_of_attack(twist, Wm, Wtheta)

    # - Calculate Inflow Angle - #
    inflow = atan(Wm, -Wtheta)

    # - Calculate Angle of Attack - #
    alpha = twist - inflow

    return alpha
end

"""
DONE. clean up and check
"""
function calculate_inflow_velocities(blade_elements, Vinf, vm, vtheta)
    Wm = similar(vm)
    Wtheta = similar(vm)

    for i in 1:nr
        Wm[:, i] = Vinf .+ vm[:, i]
        Wtheta[:, i] =
            vtheta[:, i] .- blade_elements[i].omega .* blade_elements[i].radial_position
    end

    Wmag = sqrt.(Wm .^ 2 .+ Wtheta .^ 2)

    return (Wm=Wm, Wtheta=Wtheta, Wmag=Wmag)
end

"""
DONE. clean up and check
"""
function search_polars(airfoil, alpha)
    return ccb.afeval(airfoil, alpha, 0.0, 0.0) #requires entries for Re and Ma, even though they aren't used.
end

"""
DONE. clean up and check
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

