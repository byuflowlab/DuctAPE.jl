using FiniteDiff
using ForwardDiff
# using ReverseDiff
using FLOWFoil
const ff = FLOWFoil
using CCBlade
const ccb = CCBlade

@testset "Derivatives" begin

    # set up inputs
    #---------------------------------#
    #         ROTOR Geometry          #
    #---------------------------------#

    # Blade Tip Radius, in meters
    Rtip = 10 / 2.0 * 0.0254  # inches to meters

    tip_gap = 0.0

    # Blade Hub radius, in meters
    Rhub = 0.10 * Rtip

    # number of blades
    B = 2

    # x position of rotor
    xrotor = 0.25

    # Blade section non-dimensional radial positions, chords lengths, and local twists angles in degrees
    propgeom = [
        0.15 0.130 32.76
        0.20 0.149 37.19
        0.25 0.173 33.54
        0.30 0.189 29.25
        0.35 0.197 25.64
        0.40 0.201 22.54
        0.45 0.200 20.27
        0.50 0.194 18.46
        0.55 0.186 17.05
        0.60 0.174 15.97
        0.65 0.160 14.87
        0.70 0.145 14.09
        0.75 0.128 13.39
        0.80 0.112 12.84
        0.85 0.096 12.25
        0.90 0.081 11.37
        0.95 0.061 10.19
        # 1.00 0.041 8.99
    ]

    # extract non-dimensional radial positions
    rnondim = propgeom[:, 1]
    # Dimensionalize chords
    chords = propgeom[:, 2] * Rtip
    # convert twists to radians
    twists = propgeom[:, 3] * pi / 180

    # use a NACA 4412 airfoils
    airfoils = fill(ccb.AlphaAF("test/data/naca4412.dat"), length(rnondim))

    #---------------------------------#
    #         Paneling Options        #
    #---------------------------------#
    nwake_sheets = 15

    # non-dimensional wake length
    wake_length = 1.0

    # number of panels between discrete points
    # npanels = [10; 5; 25]
    npanels = [10; 25]

    nhub_inlet = 20
    nduct_inlet = 20

    #---------------------------------#
    #       Operating Conditions      #
    #---------------------------------#

    #Vinf
    Vinf = 5.0

    # rotor rotation rate in rad/s
    Omega = 5400 * pi / 30  # convert from RPM to rad/s

    # freestream conditions
    rho = 1.225 #kg/m^3
    mu = 1.81e-5 # kg/(m⋅s)
    asound = 341.0 #m/s

    #---------------------------------#
    #      Define BODY Coordinates    #
    #---------------------------------#

    # - Duct Coordinates - #
    # use duct coordinates from FLOWFoil validation cases
    include("data/naca_662-015.jl")
    duct_coordinates = [x_duct r_duct] ./ 2.0

    # - Hub Coordinates - #
    # use hub coordinates from FLOWFoil validation cases
    # include("../test/data/bodyofrevolutioncoords.jl")
    # hub_coordinates = [x_hub[1:(end - 1)] ./ 2.0 r_hub[1:(end - 1)] * Rhub / maximum(r_hub)]
    hub_coordinates = nothing

    #---------------------------------#
    #          Define Inputs          #
    #---------------------------------#

    # Rotor Parameters
    rotor_parameters = [(;
        xrotor,
        nwake_sheets,
        r=rnondim,
        chords,
        twists,
        airfoils,
        Rtip,
        Rhub,
        tip_gap,
        B,
        Omega,
    )]

    # Paneling Parameters
    paneling_constants = (; npanels, nhub_inlet, nduct_inlet, wake_length, nwake_sheets)

    # Freestream Parameters
    freestream = (; rho, mu, asound, Vinf)

    # initialize various inputs used in analysis
    inputs = dt.precomputed_inputs(
        duct_coordinates,
        # hub_coordinates,
        nothing,
        paneling_constants,
        rotor_parameters,
        freestream;
        finterp=FLOWMath.linear,
    )

    # initialize states
    states = dt.initialize_states(inputs)


    p = (converged=[false],)

    # ForwardDiff Jacobian
    rfor = copy(states)
    rwrapfor(r, states) = dt.residual!(r, states, inputs, p)
    fordiff_j = ForwardDiff.jacobian(rwrapfor, rfor, states)

    # FiniteDiff Jacobian
    rwrapfin(r, states) = dt.residual!(r .= NaN, states, inputs, p)
    findiff_j = similar(fordiff_j) .= 0.0
    FiniteDiff.finite_difference_jacobian!(findiff_j, rwrapfin, states)

    @test isapprox(findiff_j, fordiff_j; atol=1e-3)
end
