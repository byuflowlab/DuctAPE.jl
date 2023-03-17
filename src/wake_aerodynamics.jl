"""
    calculate_enthalpy_jumps(Γr, Ωr, num_blades)

Calculate enthalpy jump across each blade.
"""
function calculate_enthalpy_jumps(Γr, Ωr, num_blades)

    # number of blade elements (or radial positions) and rotors
    nr, nrotor = size(Γr)

    # compute cumulative sum of enthalpy jumps
    H_tilde = similar(Γr) .= 0
    for ir = 1:nr
        for irotor in 1:nrotor
            Ω = Ωr[irotor]
            Γ = Γr[ir, irotor]
            B = num_blades[irotor]
            H_tilde[ir, irotor] += Ω*B*Γ/(2.0*pi)
        end
    end

    return H_tilde
end

"""
    calculate_net_circulation(Γr, num_blades)

Calculate net circulation on each blade element
"""
function calculate_net_circulation(Γr, num_blades)

    # number of blade elements (or radial positions) and rotors
    ir, nrotor = size(Γr)

    # calculate net circulations
    Γ_tilde = similar(Γr) .= 0
    for ir in 1:nr
        for ir in 1:nrotor
            Γ = Γr[ir, irotor]
            B = num_blades[irotor]
            Γ_tilde[ir, irotor] += B*Γ
        end
    end

    return Γ_tilde
end

# """
# DONE. HAS TEST. clean up
# """
# function get_surface_velocity(body_vortex_strengths, duct_panels, wake_panels)
#     ## -- Set up Edge Velocities -- ##
#     x_edge = wake_panels[1].panel_center[:, 1]
#     _, duct_inner_idx = findmin(duct_panels.panel_center[:, 1])
#     v_surf = -reverse(body_vortex_strengths[1:duct_inner_idx])
#     return fm.akima(reverse(duct_panels.panel_center[1:duct_inner_idx, 1]), v_surf, x_edge)
# end

"""
    calculate_wake_velocities(Ax_bw, Ar_bw, Γb, Ax_ww, Ar_ww, Γw)

Calculate the magnitude of the meridional velocity.  This velocity is defined tangent to a
streamline in the x-r plane such that its magnitude is given by: `Vm = Vx + Vr`
"""
function calculate_wake_velocities(Ax_bw, Ar_bw, Γb, Ax_ww, Ar_ww, Γw, Ax_rw, Ar_rw, Σr)

    # number of streamwise and radial panels
    nx, nr = size(Γw)

    # initialize meridional velocity magnitude
    Vm = similar(Γw) .= 0

    # loop through each radial set of wake panels
    for ir in 1:nr

        # add body induced velocities
        Vm[:,ir] .+= Ax_bw[ir] * Γb[ir]
        Vm[:,ir] .+= Ar_bw[ir] * Γb[ir]

        # add wake induced velocities
        for iwake in 1:size(Γw, 1)
            Vm[:, ir] .+= Ax_ww[iwake, ir] * Γw[iwake, ir]
            Vm[:, ir] .+= Ar_ww[iwake, ir] * Γw[iwake, ir]
        end

        # add rotor induced velocities
        for irotor in 1:size(Σr, 1)
            Vm[:, ir] .+= Ax_rw[irotor, ir] * Σr[irotor, ir]
            Vm[:, ir] .+= Ar_rw[irotor, ir] * Σr[irotor, ir]
        end

    end

    return Vm
end

"""
    calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, Γ_tilde, H_tilde, rotor_indices)

Calculate wake vortex strengths
"""
function calculate_wake_vortex_strengths!(Γw, wake_panels, Vm, rotor_indices)

    # number of streamwise and radial panels
    nx, nr = size(Γw)

    # loop through each radial position
    for ir in 1:nr
        # loop through each streamwise position
        for ix in 1:nx

            # find upstream rotor
            irotor = findlast(idx -> idx <= x, rotor_indices)

            # compute average meridional velocities
            if ir == 1 || ir == nr
                # use panel velocity as average velocity
                Vm_avg = Vm[ir, ix]
            else
                # average velocities above and below the current panel
                Vm_avg = 1/2 * (Vm[ir-1, ix] + Vm[ir+1, ix])
            end

            # calculate the wake vortex strength
            if Vm_avg <= 0.0
                Γw[ir,ix] = 0.0
            else
                r = wake_panels[ir].panel_center[ix, 2]
                K = -1/(8*pi^2*r^2)
                ΔΓ2 = Γ_tilde[ir+1, irotor]^2 - Γ_tilde[ir, irotor]^2
                ΔH = H_tilde[ir+1, irotor] - H_tilde[ir, irotor]
                Γw[ir,ix] = (K*ΔΓ2 + ΔH)/Vm_avg
            end
        end
    end

    return Γw
end

"""


"""
function initialize_wake_vortex_strengths(rotor_circulation_strengths, params)

    # TODO: Update this function

    # - Calculate Enthalpy Jumps - #
    H_tilde = calculate_enthalpy_jumps(rotor_circulation_strengths, params.blade_elements)

    # - Calculate Net Circulation - #
    BGamma, Γ_tilde = calculate_net_circulation(
        rotor_circulation_strengths, params.blade_elements
    )

    # - Get Vm values from thrust method in dfdc - #
    Vm = vm_from_thrust(
        params.freestream,
        params.blade_elements,
        BGamma,
        params.wake_panels,
        params.rotor_indices,
    )

    # - Initialize matrix of wake vortex strengths - #
    wake_vortex_strengths = similar(Vm)

    # - Use the normal function for calculating the vortex panel strengths - #
    calculate_wake_vortex_strengths!(
        wake_vortex_strengths,
        params.wake_panels,
        Vm,
        params.freestream.Vinf,
        Γ_tilde,
        H_tilde,
        params.rotor_indices,
    )

    # println(wake_vortex_strengths)
    return wake_vortex_strengths
end

"""
"""
function vm_from_thrust(freestream, blade_elements, BGamma, wake_panels, rotor_indices)

    # - Initialize - #

    TF = eltype(BGamma)
    nx = length(wake_panels[1].panel_center[:, 1])
    nw = length(wake_panels)

    # average meridional velocity approxmiation
    Vm = freestream.Vinf * ones(TF, blade_elements[1].num_radial_stations[1] - 1, nx)

    # accumulated axial induced velocity
    v_in = freestream.Vinf

    for i in 1:length(blade_elements)
        # r = findlast(idx -> idx <= x, rotor_indices)
        r = rotor_indices[i]

        ## -- Calculate Total Thrust of the Rotor -- ##

        annular_area = 0.0
        bg_avg = 0.0

        for j in 1:nw

            # - Blade section annular area - #
            dA =
                pi * (
                    blade_elements[i].radial_positions[j + 1]^2 -
                    blade_elements[i].radial_positions[j]^2
                )

            annular_area += dA

            # - Area Weighted Circulation - #
            # this is the numerator, we divide by the total annular area a few lines down
            bg_avg += dA * (BGamma[j + 1, r] + BGamma[j, r]) / 2.0
        end

        # disk area used in thrust calculation
        disk_area = pi * blade_elements[i].radial_positions[end]^2

        # - Total Thrust - #
        T =
            freestream.rho *
            bg_avg *
            (disk_area / annular_area) *
            blade_elements[i].omega[1] / (2.0 * pi)

        ## -- Calculate axial induced velocity of rotor -- ##

        #axial induced velocity
        #includes contributions from prior rotors
        if v_in^2 + 2.0 * T / (freestream.rho * disk_area) < 0.0
            vx = 0.1 * freestream.Vinf
        else
            vx = sqrt(v_in^2 + 2.0 * T / (freestream.rho * disk_area)) - v_in / 2.0
        end

        # increase v_in for future rotors by the inducement of this rotor
        v_in = vx

        #TODO: not sure what to do here, this doesn't converge.

        # # average meridional velocity behind rotor
        # if i < length(blade_elements)
        #     Vm[:, r:rotor_indices[i + 1]] .= freestream.Vinf + v_in
        # else
        #     Vm[:, r:end] .= freestream.Vinf + v_in
        # end

        # simply set the wake just behind the rotor and leave the rest as vinf? this does converge.
        Vm[:, r] .= freestream.Vinf + v_in

        # what about a decrease values behind rotors? this does converge, but converges to different solution than just leaving most things vinf except at the rotors
        # does this give the solver more differences to work with?
        if i < length(blade_elements)
            for x in r:(rotor_indices[i + 1] - 1)
                Vm[:, x + 1] .= Vm[:, x] .* 0.99
            end
        else
            for x in r:(nx - 1)
                Vm[:, x + 1] .= Vm[:, x] .* 0.99
            end
        end
    end

    return Vm
end
