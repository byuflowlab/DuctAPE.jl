#=

Solver Functions

Authors: Judd Mehr,

Procedure:
Iteration setup:
1. Define paneling on grid streamsurfaces, define source drag panels
2. Evaluate at body panels, rotor blade stations, source drag panels: axij , arij , aij, bxij , brij ,bij
3. Set initial guess for Γk
4. Set corresponding Γ ̃ and H ̃ fields
5. Set initial guess for γi using (41) and (42)
6. Set initial σi = 0

One Newton iteration:
1. Using current Γk, σi, evaluate γi, vxi , vri , vθi , V⃗i and derivatives w.r.t. Γk, σi 2. Evaluate residuals of equations (75), (73), (74), and derivatives w.r.t. Γk, σi
3. Solve Newton system for δΓk, δσi
4. Update Γk, σi

=#

#############################
##### ----- TYPES ----- #####
#############################

"""
    OperatingConditions{TVI,TVR,TF}

**Fields:**
 - `vinf::Array{Float}` : Array of freestream velocities
 - `vref::Array{Float}` : Array of reference velocities
 - `rho::Float` : air density value
 - `vso::Float` : speed of sound value
 - `mu::Float` : air viscosity value
"""
struct OperatingConditions{TVI,TVR,TF}
    vinf::TVI
    vref::TVR
    rho::TF
    vso::TF
    mu::TF
    # boundarylayer::TB
end

"""
"""
struct Outputs{TTt,TTd,TTr,TPt,TPr,TQt,TQr}
    totalthrust::TTt
    ductthrust::TTd
    rotorthrusts::TTr #array
    totalpower::TPt
    rotorpowers::TPr #array
    totaltorque::TQt
    rotortorques::TQr #array
end

#######################################
##### ----- ITERATION SETUP ----- #####
#######################################

"""
need to get all the a's and b's (coefficient matrices) for the walls and rotors.
Not quite sure where this is done in dfdc, but the functions talking about pointers might be a good place to start (e.g. dfdcsubs.f line 1770)
"""
function initialize_linear_system() end

####################################
##### ----- NEWTON SOLVE ----- #####
####################################

"""
probably don't need to do things the way they are done in dfdc. There are better ways in julia.  Probably use the LinearSolve.jl package and whatever other packages convenient to get the derivatives as needed.
"""
function solve_linear_system() end
