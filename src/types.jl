
##################################
####    Inputs and Options    ####
##################################

"""
"""
struct System{FS,GO,DG,WGG,RG,BD,PG,BA,G,S,PS,RV,GT}
    #Parameters
    freestream::FS
    gridoptions::GO
    #Geometry
    ductgeometry::DG
    wakegridgeometry::WGG
    rotorgeometries::RG
    bladedimensions::BD
    panelgeometries::PG
    #Aerodynamics
    bladeaerodynamics::BA
    Gammas::G
    Sigmas::S
    panelstrengths::PS
    rotorvelocities::RV
    gridthermodynamics::GT
end

"""
    Freestream{TVI,TVR,TF}

**Fields:**
 - `vinf::Float` : Freestream velocities
 - `vref::Float` : Reference velocities
 - `rho::Float` : Air density value
 - `vso::Float` : Speed of sound value
 - `mu::Float` : Air dynamic viscosity value
"""
struct Freestream{TVI,TVR,TF}
    vinf::TVI
    vref::TVR
    rho::TF
    vso::TF
    mu::TF
end

"""
    GridOptions{TF,TI}

**Fields:**
 - `num_radial_stations::Integer` : Number of radial stations (equal to number of rotor blade elements used in analysis)
 - `inlet_length::Float` : inlet length (unused)
 - `wake_length::Float` : length of wake behind duct relative to chord length
 - `wake_expansion_factor::Float` : expansion factor to apply to wake grid generation
"""
struct GridOptions{TF,TI}
    num_radial_stations::TI
    inlet_length::TF
    wake_length::TF
    wake_expansion_factor::TF
    # wakerelax::TB
end

###############################
####    System Geometry    ####
###############################

"""
    DuctGeometry{TA,TF}

**Fields:**
 - `wallinnerxcoordinates::Array{Float}` : x coordinates of inner (lower) wall geometry
 - `wallinnerrcoordinates::Array{Float}` : r coordinates of inner (lower) wall geometry
 - `wallinnerspline::FLOWMath.Akima` : Spline of inner coordinates of duct wall
 - `wallouterxcoordinates::Array{Float}` : x coordinates of outer (upper) wall geometry
 - `wallouterrcoordinates::Array{Float}` : r coordinates of outer (upper) wall geometry
 - `wallouterspline::FLOWMath.Akima` : Spline of outer coordinates of duct wall
 - `hubxcoordinates::Array{Float}` : x coordinates of hub geometry
 - `hubrcoordinates::Array{Float}` : r coordinates of hub geometry
 - `hubspline::FLOWMath.Akima` : Spline of hub coordinates
 - `LEx::Float` : x-position of leading edge
 - `TEx::Float` : x-position of trailing edge
 - `chord::Float` : chord length
 - `wallbluntTE::Bool` : flag for blunt trailing edge on wall
 - `hubbluntTE::Bool` : flag for blunt trailing edge on hub, unused
"""
struct DuctGeometry{TA,TSDi,TSDo,TSH,TF,TB}
    wallinnerxcoordinates::TA
    wallinnerrcoordinates::TA
    wallinnerspline::TSDi
    wallouterxcoordinates::TA
    wallouterrcoordinates::TA
    wallouterspline::TSDo
    hubxcoordinates::TA
    hubrcoordinates::TA
    hubspline::TSH
    LEx::TF
    TEx::TF
    chord::TF
    wallbluntTE::TB
    hubbluntTE::TB
end

"""
    RotorGeometry{TF, TI, TA, TC, TT, TAF, TR, TM}

**Fields:**
 - `xlocation::Float` : x location of rotor plane, non-dimensional based on duct chord (max TE location - min LE location of hub/wall)
 - `numblades::Int` : number of rotor blades
 - `radialstations::Array{Float}` : array of radial stations defining rotor blade, non-dimensional with hub=0, tip=1
 - `tipgap::Float` : gap between blade tip and duct wall (not implemented yet)
 - `chords::Array{Float}` : array of chord lengths at radial stations defining rotor blade, non-dimensional based on blade tip radius
 - `twists::Array{Float}` : array of twist values (in degrees) at radial stations defining rotor blade
 - `skews::Array{Float}` : array of skew values (similar to sweep) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only, the solver can't use this information)
 - `rakes::Array{Float}` : array of rake values (similar to dihedral) at radial stations defining rotor blade, non-dimensional based on rotor tip radius. (note: this is for reference only right now. it may be implemented into the grid initialization functions later.)
 - `airfoils::Array{Airfoil}` : array of airfoil data objects at radial stations defining rotor blade
 - `solidity:Array{Float}` : array of rotor solidity at radial stations defining rotor blade, chord/distance between blade sections
 - `rpm::Float` : RPM of rotor
"""
struct RotorGeometry{TF,TI,TR,TG,TC,TT,TSk,TRa,TAF,TSo,TRpm}
    xlocation::TF
    numblades::TI
    radialstations::TR
    tipgap::TG
    chords::TC
    twists::TT
    skews::TSk
    rakes::TRa
    airfoils::TAF
    solidities::TSo
    rpm::TRpm
end

"""
    BladeDimensions{TF, TA}

**Fields:**
 - `rhub::Float` : hub radius (dimensional)
 - `rtip::Float` : tip radius (dimensional)
 - `rdim::Array{Float}` : array of dimensional radial stations
 - `cdim::Array{Float}` : array of dimensional chords
 - `tdim::Array{Float}` : array of twists (already dimensional in rotorgeometry)
 - `sweptannulus::Float` : area of blade swept annulus
 - `sweptarea::Float` : area of blade tip swept disk
"""
struct BladeDimensions{TF,TR,TC,TT}
    hubr::TF
    tipr::TF
    rdim::TR
    cdim::TC
    tdim::TT
    sweptannulus::TF
    sweptarea::TF
end

"""
    WakeGridGeometry{TF,TI,TA,TW,TH}

Wake grid geometry object

**Fields:**
 - `x_grid_points::Matrix{Float}` : 2D Array of x grid points
 - `r_grid_points::Matrix{Float}` : 2D Array of radial grid points
 - `nx::Int` : number of x stations
 - `nr::Int` : number of radial stations
 - `wallTEidx::Int` : index of duct wall trailing edge x location
 - `hubTEidx::Int` : index of hub wall trailing edge x location
 - `rotoridxs::Array{Int}` : array of indices of rotor x locations
"""
struct WakeGridGeometry{TF,TI,TA}#,TW,TH}
    x_grid_points::TF
    r_grid_points::TF
    nx::TI
    nr::TI
    wallTEidx::TI
    hubTEidx::TI
    rotoridxs::TA
    # wall_xstations::TW
    # hub_xstations::TH
end

"""
    Panels{TPEx,TPEr,TPC,TPN,TPT}

**Fields:**
 - `panel_edges_x::Array{Tuple{Float, Float}}` : Array of sets of x locations for panel edges
 - `panel_edges_r::Array{Tuple{Float, Float}}` : Array of sets of r locations for panel edges
 - `panel_centers::Array{Tuple{Float, Float}}` : Array of sets of x,r locations for panel centers
 - `panel_normals::Array{Tuple{Float, Float}}` : Array of panel unit normal vectors
 - `panel_tyes::Array{String}` : Array of panel types (for use in assembling linear system)
"""
struct Panels{TPEx,TPEr,TPC,TPN,TPT}
    panel_edges_x::TPEx
    panel_edges_r::TPEr
    panel_centers::TPC
    panel_normals::TPN
    panel_types::TPT
end

"""
    PanelGeometries{TD,TH,TW,TR}

**Fields:**
 - `wall_panels::DuctTAPE.Panels` : panels defining duct wall airfoil
 - `hub_panels::DuctTAPE.Panels` : panels defining hub
 - `wake_panels::DuctTAPE.Panels` : panels defining rotor wake vortex sheets
 - `rotor_source_panels::DuctTAPE.Panels` : panels defining rotor drag source panels
"""
struct PanelGeometries{TD,TH,TW,TR}
    wall_panels::TD
    hub_panels::TH
    wake_panels::TW
    rotor_source_panels::TR
    #drag_panels::TD
end

###################################
####    System Aerodyanmics    ####
###################################

"""
    SystemAero{TG,TH,TS,TC,TR}

**Fields:**
 - `b_gamma_grid::Matrix{Float}` : B*Γ values on wake grid
 - `delta_enthalpy_grid::Matrix{Float}` : ΔH (enthalpy) values on wake grid
 - `delta_entropy_grid::Matrix{Float}` : ΔS (entropy) values on wake grid
 - `b_circ_rotors::Array{Array{Float}}` : B*γ values at rotor blades (one array per rotor in increasing order of x location)
 - `rotor_source_strengths::Array{Array{Float}}` : σ (source strength) values at rotor blades (one array per rotor in increasing order of x location)
 - `control_point_velocities::Matrix{Float}` : velocities at control points
"""
struct SystemAero{TG,TH,TS,TC,TR,TV}
    b_gamma_grid::TG
    delta_enthalpy_grid::TH
    delta_entropy_grid::TS
    b_circ_rotors::TC
    rotor_source_strengths::TR
    control_point_velocities::TV
end

"""
    BladeAero{TRe,TMa,TCl,TCd,TM,TG,TW,TVa,TVt}

Blade aerodynamic values.
TODO: unsure where to put these, or if they even need to be put anywhere.

**Fields:**
 - `reynolds::Array{Float}` : local section reynolds numbers
 - `mach::Array{Float}` : local section mach numbers
 - `cl::Array{Float}` : local section coefficients of lift
 - `cd::Array{Float}` : local section coefficients of drag
 - `cm::Array{Float}` : local section coefficients of moment
 - `Gamma::Array{Float}` : local section circulations
 - `W::Array{Float}` : local inflow velocity
 - `vax::Array{Float}` : local axial velocity
 - `vtan::Array{Float}` : local tangential velocity
"""
struct BladeAero{TRe,TMa,TCl,TCd,TM,TG,TW,TVa,TVt}
    reynolds::TRe
    mach::TMa
    cl::TCl
    cd::TCd
    cm::TM
    Gamma::TG
    W::TW
    vax::TVa
    vtan::TVt
end

"""
    RotorVelocities{TA}

**Fields:**
 - `induced_axial_velocities::Array{Float}` : local section induced axial velocities
 - `induced_radial_velocities::Array{Float}` : local section induced radial velocities
 - `induced_tangential_velocities::Array{Float}` : local section induced tangential (circumferential) velocities
 - `absolute_axial_velocities::Array{Float}` : local section absolute axial velocities
 - `absolute_radial_velocities::Array{Float}` : local section absolute radial velocities
 - `absolute_tangential_velocities::Array{Float}` : local section absolute tangential velocities
 - `relative_axial_velocities::Array{Float}` : local section relative axial velocities
 - `relative_radial_velocities::Array{Float}` : local section relative radial velocities
 - `relative_tangential_velocities::Array{Float}` : local section relateive tangential velocities
"""
struct RotorVelocities{TA}
    induced_axial_velocities::TA
    induced_radial_velocities::TA
    induced_tangential_velocities::TA
    absolute_axial_velocities::TA
    absolute_radial_velocities::TA
    absolute_tangential_velocities::TA
    relative_axial_velocities::TA
    relative_radial_velocities::TA
    relative_tangential_velocities::TA
end

#######################
####    Outputs    ####
#######################
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

"""
"""
struct Forces{TT,TQ}
    thrust::TT
    torque::TQ
end
