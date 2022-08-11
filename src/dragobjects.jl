#=
Types and Functions related to stationary drag objects

Authors: Judd Mehr,
=#

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
