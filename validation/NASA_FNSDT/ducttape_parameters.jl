#=
This is a file containing all the parameters that get input into DuctAPE.
It doesn't need to be its own file, but sometimes the parameters can get quite lengthy, especially if you include the duct and hub coordinates.
=#

#---------------------------------#
#   Duct + Centerbody Geometry    #
#---------------------------------#
#= I'll first include the duct and hub coordinates, which in this case I have stashed down a few levels where they were saved from parsing the NASA geometry.
TODO: need to figure out the hub trailing edge, may need to modify (the readme indicates that this shouldn't be a problem) and keep a different file closer to this one.
Note that the geometry for DuctAPE is required to be clockwise, such that the duct geometry starts at the trailing edge, proceeds along the inner surface, then around the leading edge, along the outer surface, and ends back at the trailing edge.
Similarly, the centerbody geometry starts at the leading edge and proceeds to the trailing edge
=#
include(
    "../geometry_parsing/parsed_geometry_files/extracted_duct_and_centerbody_geometry_in_julia_arrays.jl",
)
hub_coordinates = hub_coordinates_closed_TE

# - Get rotor tip and hub information from the body data - #
# need to spline things and get the center body and inner duct wall r locations at the rotor x location
_, leid = findmin(duct_coordinates[:, 1])
inner_duct = reverse(duct_coordinates[1:leid, :]; dims=1)

xrotor = 0.0 # axial position of rotor stacking axis (dimensional)
Rtip_rotor = FLOWMath.akima(inner_duct[:, 1], inner_duct[:, 2], xrotor)
Rhub_rotor = FLOWMath.akima(hub_coordinates[:, 1], hub_coordinates[:, 2], xrotor)

xstator = 7.26 * 0.0254 #location of stator stacking axis
Rtip_stator = FLOWMath.akima(inner_duct[:, 1], inner_duct[:, 2], xstator)
Rhub_stator = FLOWMath.akima(hub_coordinates[:, 1], hub_coordinates[:, 2], xstator)

#---------------------------------#
#    Rotor + Stator Parameters    #
#---------------------------------#
#=
The rotor parameters are contained in a vector of named tuples, such that each named tuple (the name matters! DuctAPE expects the naming conventions presented here) has the information for a single rotor.  These should be included in the order they appear axially, e.g. rotor then stator
Note that the rotor parameters take in twist, not stagger, stagger is calculated from twist internally.  In our case, stagger is the angle from the axial direction, and twist is the angle from the plane of rotation, so they are just off by (90-theta) degrees
=#

##### ----- ROTOR ----- #####
B_rotor = 22 # number of rotor blades
# RPM_rotor = 7808.0 #mechanical RPM in hotwire data readme
RPM_rotor = 1000 #mechanical RPM in hotwire data readme
RPMc_rotor = 12657.0 * 0.61 #corrected RPM from aero paper
Omega_rotor = RPM_rotor * pi / 30.0 #TODO: which RPM to use?
# Omega_rotor = RPMc_rotor * pi / 30.0

# sections to include
rotorids2include = 1:25

rotor_airfoils = []

sidxoverwrite = 1
for sidx in rotorids2include
    include(
        "xrotor_airfoil_parameter_approximations/rotor_section$(isnothing(sidxoverwrite) ? sidx : sidxoverwrite)_dfdc_params.jl",
    )
    push!(
        rotor_airfoils,
        dt.DFDCairfoil(;
            alpha0=alpha0,
            clmax=clmax,
            clmin=clmin,
            dclda=dclda,#
            dclda_stall=0.1,# dclda_stall,
            dcl_stall=dcl_stall,
            cdmin=cdmin,
            clcdmin=clcdmin,
            dcdcl2=dcdcl2,
            cmcon=0.0, #cmom
            Re_ref=re_ref,
            Re_exp=re_exp,
            mcrit=mcrit,
        ),
    )
end

# read in rotor geometry
include("rotor_geometry.jl")

rotor_parameters = (;
    xrotor=xrotor,
    r=rotor_rct[rotorids2include, 1] ./ Rtip_rotor, #non-dimensionalize
    chords=rotor_rct[rotorids2include, 2],
    twists=rotor_rct[rotorids2include, 3],
    airfoils=rotor_airfoils,
    Rtip=Rtip_rotor,
    Rhub=Rhub_rotor,
    tip_gap=0.0,
    B=B_rotor,
    Omega=Omega_rotor,
    fliplift=false,
)

# ##### ----- STATOR ----- #####
B_stator = 54 # number of stator blades

# sections to include
# statorids2include = 1:25
statorids2include = 15

# statorids2include = [
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# 10
# # 11 # can't get a good fit for this one.
# 12
# 13
# 14
# 15
# 16
# 17
# 18
# 19
# 20
# 21
# 22
# 23
# 24
# 25
# ]

stator_airfoils = []

for sidx in statorids2include
    include("xrotor_airfoil_parameter_approximations/stator_section$(sidx)_dfdc_params.jl")
    push!(
        stator_airfoils,
        dt.DFDCairfoil(;
            alpha0=alpha0,
            clmax=clmax,
            clmin=clmin,
            dclda=dclda,#
            dclda_stall=0.1,# dclda_stall,
            # dclda_stall= dclda_stall,
            dcl_stall=dcl_stall,
            cdmin=cdmin,
            clcdmin=clcdmin,
            dcdcl2=dcdcl2,
            cmcon=0.0, #cmom
            Re_ref=re_ref,
            Re_exp=re_exp,
            mcrit=mcrit,
        ),
    )
end

# read in stator geometry files
include("stator_geometry.jl")

if length(stator_airfoils) == 1
    stator_airfoils = fill(stator_airfoils[1], length(stator_rct[:, 1]))
    statorids2include = 1:length(stator_airfoils)
end

stator_parameters = (;
    xrotor=xstator,
    r=stator_rct[statorids2include, 1] ./ Rtip_stator, #only used for interpolation
    chords=stator_rct[statorids2include, 2],
    twists=stator_rct[statorids2include, 3],
    airfoils=stator_airfoils,
    Rtip=Rtip_stator, # unused right now, duct position set by leading rotor
    Rhub=Rhub_stator, # unused right now, duct position set by leading rotor
    tip_gap=0.0, #unused right now
    B=B_stator,
    Omega=0.0, # USED, make sure stator is set to zero rotation rate
    fliplift=true,
)

# - put rotor and stator parameters together into a vector - #
rotorstator_parameters = [rotor_parameters, stator_parameters]
# rotorstator_parameters = [rotor_parameters]

#---------------------------------#
#         Paneling Constants      #
#---------------------------------#
#=
DuctAPE repanels the duct and center body geometries in such a way to allow nice alignment of the body and wake panels.  In addition, a set number of panels between discrete locations needs to be maintained for gradient-based optimization purposes.
For that reason, we need to tell DuctAPE how many panels to put between inlet and rotor, rotor and stator, stator and trailing edge, and trailing edge and wake.  If the duct and center body trailing edges don't line up, then we also need to prescribe how many panels go between trailing edges.
We also need to tell DuctAPE how long to make the wake aft of the trailing edge.  Typically, 1 chord length is used, but perhaps 1 rotor radius might be better if that length is longer than the duct.
In addition, the same x-spacing used on the inner duct wall panels is used for the outer duct wall panels.
=#

npanref = 1.0 # use this to manually refine things uniformly

paneling_constants = (;
    # wake_length=1.0, # OPEN
    wake_length=0.025, # CLOSED
    nwake_sheets=ceil(Int, npanref * 11), #number of wake sheets
    nduct_inlet=ceil(Int, npanref * 13), #number of panels to use between the duct leading edge and first rotor
    nhub_inlet=ceil(Int, npanref * 10), #number of panels between centerbody leading edge and first rotor
    npanels=[
        ceil(Int, npanref * 5), #panels between rotor and stator
        ceil(Int, npanref * 10), #panels between stator and duct trailing edge
        ceil(Int, npanref * 35), #CLOSED panels between duct and center body trailing edge
        # ceil(Int, npanref * 10), #OPEN panels between duct and center body trailing edge
        ceil(Int, npanref * 1), #CLOSED panels between trailing edge and end of wake
        # ceil(Int, npanref * 25), #OPEN panels between trailing edge and end of wake
    ],
)

#---------------------------------#
#            Freestream           #
#---------------------------------#
Ma = 0.05 # NASA aero paper claimed all runs were done at Ma=0.05
asound = 341.0
Vinf = Ma * asound
rhoinf = 1.225 # choose values that match asound
muinf = 1.81e-5 # choose values that match asound
freestream = (; Vinf, muinf, rhoinf, asound)

#---------------------------------#
#      Reference Parameters       #
#---------------------------------#
#=
Reference parameters are used in the post processing functions
=#
reference_parameters = (; Vref=Vinf, Rref=Rtip_rotor)

#######################################################################
##                                                                    #
##                         Write a DFDC File                          #
##                                                                    #
#######################################################################
#filename = "nasa.case"

#op_data = (; rhoinf, muinf, Vso=asound, Vinf=Vinf, Vref=Vinf, Alt=0.0, RPM=RPM_rotor)

#wake_data = (; nwake=5, xwake=0.1, rlx_wake="F\n", nwake_sheets=15)

#rct = rotor_rct

#airfoil_data = []
#for sidx in rotorids2include
#    if sidx % 2 == 0
#        include(
#            "xrotor_airfoil_parameter_approximations/rotor_section$(sidx)_dfdc_params.jl"
#        )
#        push!(
#            airfoil_data,
#            (;
#                xisection=rct[sidx, 1] ./ Rtip_rotor,
#                alpha0,
#                clmax,
#                clmin,
#                dclda,
#                dclda_stall,
#                dcl_stall,
#                cdmin,
#                clcdmin,
#                dcdcl2,
#                cmcon=0.0, #cmom
#                Re_ref=re_ref,
#                Re_exp=re_exp,
#                mcrit,
#            ),
#        )
#    end
#end
#airfoil_data = [airfoil_data]

#include("xrotor_airfoil_parameter_approximations/rotor_section15_dfdc_params.jl")
#airfoil_data = [(;
#    xisection=0.0,
#    alpha0,
#    dclda,
#    clmax,
#    clmin,
#    dclda_stall,
#    dcl_stall,
#    cmcon=0.0, #cmom
#    mcrit,
#    cdmin,
#    clcdmin,
#    dcdcl2,
#    Re_ref=re_ref,
#    Re_exp=re_exp,
#)]

#rotor_data = [(;
#    # naf=length(airfoil_data[1]),
#    naf=1,
#    xrotor,
#    B=B_rotor,
#    r=rct[:, 1],
#    chord=rct[:, 2],
#    twist=rct[:, 3] * 180 / pi,
#)]
