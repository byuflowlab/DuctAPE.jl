#=
Set up type drafts until better locations can be found
=#

#TODO: need to decide on units vs non-dimensional approach
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

"""
    DragObject{TX,TR,TA}

**Fields:**
 - `xs::Array{Float}` : array of x locations of drag object
 - `rs::Array{Float}` : array of r locations of drag object
 - `cdas::Array{Float}` : array of drag areas for drag object
end
"""
struct DragObject{TX,TR,TA}
    xs::TX
    rs::TR
    cdas::TA
end
