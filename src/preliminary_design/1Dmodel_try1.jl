#=
Attempt at a "1D" model in order to provide a reasonable starting point for design/optimization

Goal for option 1: input open rotor design, output duct inlet and outlet areas
Goal for option 2: input some performance metric (say, thrust), and some sizing guide (say rotor diameter), and some operating conditions (say, Vinf and Omega), and output reasonable duct inlet and outlet as well as some rotor geometry (chord, twist, camber?)

=#

#=
#-----------------------------------#
# Potentially useful relationships  #
#-----------------------------------#

Thrust = mdot * (exit_velocity - freestream_velocity) (conservation of momentum)

=#


######################################################################
#                                                                    #
#                         1D Model Option 1                          #
#                                                                    #
######################################################################

#=

Known (or can be approximated): Rtip, Rhub, chord, twist, airfoil parameters, reynolds number, rotation rate, number of blades, open rotor performance data (CT, CP, CQ, eta), freestream conditions (Vinf, rhoinf, muinf, asound).

Want to Find: inlet area, outlet area for reasonable ducted performance (at least as much thrust and efficiency as without the duct if possible)

=#

"""
losses due to the inlet?  what losses? what counts as the inlet?
"""
function inlet_loss() end

"""
losses between stations inside duct, e.g. between the inlet and fan face.
"""
function duct_loss() end

function fanface_velocity()
    return Vinf - inlet_loss() - duct_loss()
    # ?Q? what about velocity increases due to duct? are losses negative?
end

"""
not sure what's happening here, but seems like we might need something along these lines
"""
function pressure_from_open_ct(CT, disk_area)
    T = CT * rho * n^2 * disk_area^4 # is this valid
    dp = T / disk_area

    return dp
end


"""
probably need a multi-dispatch here and then can handle both options maybe.
rotor effects on the flow, takes in velocity going in, and gives velocity coming out of the fan?
"""
function rotor_effects() end

"""
stators only add losses? what losses? aren't there some gains by de-rotating the flow?
"""
function stator_effects() end

"""
probably going to need this in some form
"""
function converve_mass() end

"""
probably going to need this in some form
"""
function conserve_momentum() end

"""
function to call the other functions and find what reasonable inlet and outlet areas might be
"""
function find_io_areas(
    Rtip, # rotor tip radius
    Rhub, # rotor hub radius
    radial_stations, # radial stations between Rhub and Rtip
    chord, # chord distribution
    twist, # twist distribution (better for this to be stagger?)
    ambient_static_pressure=101325.0, # Pa standard at sea level
)

    # solver variable initial guesses, just set to disk annulus area for now
    inlet_area = pi * (Rtip^2 - Rhub^2)
    outlet_area = inlet_area

    #TODO: choose a reasonable solver approach (NLSolve?)

    return res.zero
end

######################################################################
#                                                                    #
#                         1D Model Option 2                          #
#                                                                    #
######################################################################
