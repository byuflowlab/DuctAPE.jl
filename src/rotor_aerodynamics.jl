#=

Functions regarding rotor aerodynamics

Authors: Judd Mehr,

=#

"""
DONE. clean up and check
"""
function calculate_induced_velocities_on_rotors(
    BGamma, #matrix size num blade elem x num rotors
    Gamma_tilde, #matrix size num blade elem x num rotors
    blade_elements,
    vxd_bodies_to_rotor, #vector of matrices
    vrd_bodies_to_rotor, #vector of matrices
    gamma_bodies,
    vxd_wake_to_rotor, #matrix of matrices size num wakes x num rotors
    vrd_wake_to_rotor, #matrix of matrices size num wakes x num rotors
    gamma_wake, # matrix
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
        # TODO: need to double check the theory and implementation here
        vm[:, i] .= (vxd_bodies_to_rotor[i] .+ vrd_bodies_to_rotor[i]) * gamma_bodies

        # - Add Wake Induced Velocities - #
        for w in 1:length(gamma_wake[:, 1])
            vx = vxd_wake_to_rotor[w, i] * gamma_wake[w, :]
            vr = vrd_wake_to_rotor[w, i] * gamma_wake[w, :]
            vm[:, i] .+= sqrt(vx^2 + vr^2)
        end

        # # - Add Rotor Induced Velocities - #
        # for j in 1:nr
        #     vm[:, i] += A_rotor_to_rotor[i, j] * Sigma[i]
        # end

        #vtheta comes directly from circulation

        if i == 1
            vtheta[:, i] .=
                1.0 / (2.0 * pi * blade_elements[i].radial_positions) * (0.5 * BGamma[:, i])
        else
            vtheta[:, i] .=
                1.0 / (2.0 * pi * blade_elements[i].radial_positions) *
                (Gamma_tilde[:, i - 1] + 0.5 * BGamma[:, i])
        end
    end

    return (vm=vm, vtheta=vtheta)
end

"""
DONE.HAS TEST. clean up
twist is from plane of rotation
TODO: need to check upstream in the functions that the sign convention for the velocities matches expectations for this function.
"""
function calculate_angle_of_attack(twist, Wm, Wtheta)

    # - Calculate Inflow Angle - #
    inflow = atan.(Wm, -Wtheta)

    # - Calculate Angle of Attack - #
    alpha = twist - inflow

    # if typeof(alpha) == Float64
    #     println("section twist degrees: ", twist * 180 / pi)
    #     println("inflow angle degrees: ", inflow * 180 / pi)
    #     println("angle of attack degrees: ", alpha * 180 / pi)
    # end

    return alpha
end

"""
DONE.HAS TEST. clean up
vm and vtheta are matrices of same size as Gamma_tilde, etc. meaning there is a value at each of the blade elements on each of the rotors.
"""
function calculate_inflow_velocities(blade_elements, Vinf, vm, vtheta)
    nr = length(blade_elements)

    Wm = similar(vm)
    Wtheta = similar(vm)

    for i in 1:nr
        Wm[:, i] = Vinf .+ vm[:, i]
        Wtheta[:, i] =
            vtheta[:, i] .- blade_elements[i].omega .* blade_elements[i].radial_positions
    end

    Wmag = sqrt.(Wm .^ 2 .+ Wtheta .^ 2)

    return (Wm=Wm, Wtheta=Wtheta, Wmag=Wmag)
end

"""
DONE.HAS TEST. clean up
alpha is in radians
TODO: this is just a place holder. need to add more comprehensive function later
"""
function search_polars(airfoil, alpha)
    return ccb.afeval(airfoil, alpha, 0.0, 0.0) #requires entries for Re and Ma, even though they aren't used.
end

"""
DONE.HAS TEST. clean up
"""
function calculate_gamma_sigma(blade_elements, Vinf, vm=-1.0, vtheta=-1.0)

    # - Rename for Convenience - #
    TF = eltype(blade_elements[1].chords)

    # - Initialize Outputs - #
    Gamma = zeros(TF, blade_elements[1].num_radial_stations[1], length(blade_elements))
    Sigma = zeros(TF, blade_elements[1].num_radial_stations[1], length(blade_elements))

    # - Call function - #
    calculate_gamma_sigma!(Gamma, Sigma, blade_elements, Vinf, vm, vtheta)

    # - Return Initialized Values - #
    return Gamma, Sigma
end

"""
DONE.HAS TEST. clean up
"""
function calculate_gamma_sigma!(Gamma, Sigma, blade_elements, Vinf, vm=-1.0, vtheta=-1.0)

    # - Rename for Convenience - #
    TF = eltype(Gamma)

    # - Initialize reused vectors - #
    cl = zeros(TF, length(blade_elements[1].chords))
    cd = zeros(TF, length(blade_elements[1].chords))

    # - Calculate Inflow Velocity - #
    if vm == -1.0 && vtheta == -1.0
        #first run through, no induced velocities
        vm = zeros(TF, length(blade_elements[1].chords), length(blade_elements))
        vtheta = zeros(TF, length(blade_elements[1].chords), length(blade_elements))
    end

    inflow = calculate_inflow_velocities(blade_elements, Ref(Vinf), vm, vtheta)

    # - Loop through rotors - #
    for i in 1:length(blade_elements)

        # - Calculate Angle of Attack - #
        alphas =
            calculate_angle_of_attack.(
                blade_elements[i].twists, inflow.Wm[:, i], inflow.Wtheta[:, i]
            )

        # - Look up lfit and drag data - #
        for a in 1:length(alphas)
            cl[a], cd[a] = search_polars(blade_elements[i].airfoils[a], alphas[a])
        end

        # if eltype(cl) == Float64
        #     println("lift coeffs: ", cl)
        # end

        # - Calculate Gamma and Sigma - #
        @. Gamma[:, i] = 0.5 * inflow.Wmag[:, i] * blade_elements[i].chords * cl
        @. Sigma[:, i] =
            blade_elements[i].num_blades / (4 * pi * blade_elements[i].radial_positions) *
            inflow.Wmag[:, i] *
            blade_elements[i].chords *
            cd

        # if eltype(Gamma) == Float64
        #     println("Gamma: ", Gamma)
        # end

    end

    return nothing
end

