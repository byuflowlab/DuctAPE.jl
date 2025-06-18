#TODO: only used somewhere in docs, remove when you can
function cosine_spacing(N=80)
    return scaled_cosine_spacing(N, 1.0, 0.0; mypi=pi)
end

# TODO: probably doesn't need to be so complicated, see if nominal cosine spacing can be used and move it to the file it's used in.
function scaled_cosine_spacing(N, scale, transform; mypi=pi)
    return transform .+ scale * [0.5 * (1 - cos(mypi * (i - 1) / (N - 1))) for i in 1:N]
end

# TODO: if this is only used in post-processing, probably move it there.
"""
    split_bodies!(
        casing_vec, nacelle_vec, cb_vec, casing_zpts, nacelle_zpts, cb_zpts, vec, controlpoint
    )

Split full duct and center body surface distributions into casing, nacelle, and center body.

# Arguments
- `casing_vec::Vector{Float}` : vector for casing values to be updated in place.
- `nacelle_vec::Vector{Float}` : vector for nacelle values to be updated in place.
- `cb_vec::Vector{Float}` : vector for center body values to be updated in place.
- `casing_zpts::Vector{Float}` : axial points for casing to be updated in place.
- `nacelle_zpts::Vector{Float}` : axial points for nacelle to be updated in place.
- `cb_zpts::Vector{Float}` : axial points for center body to be updated in place.
- `vec::Vector{Float}` : vector of surface distrubtion values (e.g. pressues, or velocities, etc.)
- `controlpoint::Matrix{Float}` : control points from which the zpts are taken
"""
function split_bodies!(
    casing_vec, nacelle_vec, cb_vec, casing_zpts, nacelle_zpts, cb_zpts, vec, controlpoint
)

    # get dimensions
    casing_ids = 1:size(casing_vec, 1)
    nacelle_ids = (casing_ids[end] + 1):(casing_ids[end] + size(nacelle_vec, 1))
    cb_ids = (nacelle_ids[end] + 1):(nacelle_ids[end] + size(cb_vec, 1))

    casing_vec .= @view(vec[casing_ids, :])
    nacelle_vec .= @view(vec[nacelle_ids, :])
    cb_vec .= @view(vec[cb_ids, :])
    casing_zpts .= @view(controlpoint[1, casing_ids])
    nacelle_zpts .= @view(controlpoint[1, nacelle_ids])
    cb_zpts .= @view(controlpoint[1, cb_ids])

    return casing_vec, nacelle_vec, cb_vec, casing_zpts, nacelle_zpts, cb_zpts
end

# only used in a test: TODO: update test and remove this function
function split_bodies(vec, controlpoint, endpanelidxs; duct=true, center_body=true)
    # get type of vector for consistent outputs
    TF = eltype(vec)

    #check if duct is used
    if !duct
        #center_body only
        return TF[], TF[], vec, TF[], TF[], controlpoint[:, 1]
    else
        # split duct into inner and outer
        ndpan = Int(endpanelidxs[2, 1])
        # get duct leading edge index. assumes duct comes first in vector
        _, leidx = findmin(controlpoint[1, 1:ndpan])
        if !center_body
            #duct only
            return vec[1:leidx, :],
            vec[(leidx + 1):ndpan, :],
            TF[],
            controlpoint[1, 1:leidx],
            controlpoint[1, (leidx + 1):ndpan],
            TF[]
        else
            #duct and center_body
            return vec[1:leidx, :],
            vec[(leidx + 1):ndpan, :],
            vec[(ndpan + 1):end, :],
            controlpoint[1, 1:leidx],
            controlpoint[1, (leidx + 1):ndpan],
            controlpoint[1, (ndpan + 1):end]
        end
    end

    # shouldn't get to this point...
    return nothing
end
