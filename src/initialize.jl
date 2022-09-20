#=
Functions for initialization of the ducted rotor system before the newton solve

Authors: Judd Mehr,
=#

"""

Generates Wake Grid and System Paneling
"""
function initialize_geometry(ductgeometry, rotors, gridoptions)

    # Create Wake Grid
    wakegrid, updatedrotors = generate_wake_grid(ductgeometry, rotors, gridoptions)

    # Create Panels
    systempanels = generate_panel_system(ductgeometry, updatedrotors, wakegrid)

    return wakegrid, systempanels, updatedrotors
end

"""

Sets initial flow conditions for the system.

Runs a panel method on mirrored duct geometry to get initial guess for flow inside duct.
Runs CCBlade using panel method inflow to get first guess for rotor performance.
Then initializes guesses for the blade circulation, wake vorticity, and rotor source strengths.
"""
function initialize_flow_data(systempanels)

    # Run Panel Method to get first guess on flow field.

    # Run CCBlade to get W and cl values

    # Set Blade Circulation

    # Set Wake Enthalpy

    # Set Wake Vorticity

    # Set Rotor Source Strengths

    # Generate Grid Flow Data Object

    # Generate Panel Strengths Object

    return gridflowdata, panelstrengths
end

"""

is this needed in rework?
"""
function initialize_system()

    # Initialize Rotor and Wake Aerodynamics

    system_aero, rotor_velocities, initial_vaxial_average = initialize_system_aerodynamics(
        rotors, blades, wakegrid, freestream; niter=10, rlx=0.5
    )

    # Initialize average V_m
    vm_average = initialize_vm_average(initial_vaxial_average)

    #see calculate gamma_theta_i values on wake vortex sheet panesl
    gamma_theta_wakes = calculate_gamma_theta(system_aero, vm_average)

    #NOTE: gamth = gth after this point in dfdc
    #
    #see GAMSOLV
    #
    # Create some sort of system object, make sure convergence flag is defaulted to false.
    #TODO: figure out what to return
    return system
end
