using DelimitedFiles
using LinearAlgebra
using FLOWMath
const fm = FLOWMath
include("plots_default_new.jl")
pyplot()

######################################################################
#                                                                    #
#                           FUNCTIONS USED                           #
#                                                                    #
######################################################################

#---------------------------------#
#      2D (Duct, Centerbody)      #
#---------------------------------#
"""
    parse2d(file; TF=Float64)

Parse .xr files from NASA data files
"""
function parse2d(file; TF=Float64)
    # need this since they are initialized inside the loop
    local x
    local r

    open(file) do f
        # files tell us how many coordinates there are
        N = parse(Int, readline(f))
        # initialize output arrays
        x = zeros(N)
        r = zeros(N)

        # loop through and parse lines to get x and r coordinates
        for (i, line) in enumerate(eachline(f))
            parts = split(line)
            x[i] = parse(TF, parts[1])
            r[i] = parse(TF, parts[2])
        end
    end

    # return x,r coordinates separately
    return x, r
end

"""
    write2d(filename, varname, coordinates, writetype)

write data parsed from NASA .xr file into julia arrays for easy loading
"""
function write2d(filename, varname, coordinates, writetype)
    f = open(filename, writetype)

    write(f, varname * " = [\n")
    for (i, c) in enumerate(eachrow(coordinates))
        write(f, "    $(c[1]) $(c[2])\n")
    end
    write(f, "]\n\n")
    close(f)

    return nothing
end

#---------------------------------#
#        3D (rotor, stator)       #
#---------------------------------#
"""
    parse3d(file; TF=Float64)

parse .xyz files into a 3d array of #sections, #coordinates, and xyz
"""
function parse3d(file; TF=Float64)
    lines = DelimitedFiles.readdlm(file)
    ncoord = lines[1, 1]
    nsec = lines[1, 2]
    sections = zeros(nsec, ncoord, 3)
    for i in 1:nsec
        sections[i, :, :] .= lines[((i - 1) * ncoord + 2):(ncoord * i + 1), :]
    end

    return sections
end

"""
    rotor_edges(sections; edge="LE")

Find the leading or trailing edges of the raw rotor sections, meaning the furthest foreward and aft axial coordinate locations.

returns the associated x and r (not z) coordinates along with the associated coordinate index
"""
function rotor_edges(sections; edge="LE")
    ns, nc, _ = size(sections)

    rotor_x = zeros(ns)
    rotor_r = zeros(ns)
    rotor_id = zeros(Int, ns)

    for is in 1:ns
        # get the section coordinates
        coordinates = sections[is, :, :]

        # find foremost point in the section
        if edge == "LE"
            rotor_x[is], rotor_id[is] = findmin(coordinates[:, 1])
        else
            rotor_x[is], rotor_id[is] = findmax(coordinates[:, 1])
        end
        # get associated radial position
        rotor_r[is] = sqrt(coordinates[rotor_id[is], 2]^2 + coordinates[rotor_id[is], 3]^2)
    end

    return rotor_x, rotor_r, rotor_id
end

"""
    interpolate_sections(sections, zte; nsections=25, dimrange=nothing, cutoff=0.05)

interpolates the give sections at the chosen number of sections, not including those cutoff on the ends due to the mismatch with the body walls.

nominally takes in the sections, the section trailing edge z locations, outputs interpolated x,y,z data

the user may also choose to input manually the dimensional range at which to sample.  In this case, one could use the r coordinates, as long as zte is also input as the r coordinates.
"""
function interpolate_sections(sections, zte; nsections=25, dimrange=nothing, cutoff=0.05)
    ns, nc, nd = size(sections)

    # - Define the dimensional range at which to sample, and initizliae the outputs - #
    if dimrange == nothing
        # linearly space if nothing is input
        interprange = range(cutoff, 1.0 - cutoff, nsections)
        dimrange =
            interprange .* (maximum(sections[:, 1, 3]) - minimum(sections[:, 1, 3])) .+
            sections[1, 1, 3]
        interp_sections = zeros(nsections, nc, nd)
    else
        # just use user inputs otherwise
        interp_sections = zeros(length(dimrange), nc, nd)
    end

    for is in 1:nsections
        # for each interpolated section, find the raw sections right after and right before
        io = min(ns, searchsortedfirst(zte, dimrange[is]))
        ii = max(1, io - 1)
        inner_fraction = (dimrange[is] - zte[ii]) / (zte[io] - zte[ii])

        # for each coordinate, interpolate linearly between the two nearest data stations
        for ic in 1:nc
            interp_sections[is, ic, 3] = fm.linear(
                [0.0; 1.0], [sections[ii, ic, 3]; sections[io, ic, 3]], inner_fraction
            )
            interp_sections[is, ic, 2] = fm.linear(
                [0.0; 1.0], [sections[ii, ic, 2]; sections[io, ic, 2]], inner_fraction
            )
            interp_sections[is, ic, 1] = fm.linear(
                [0.0; 1.0], [sections[ii, ic, 1]; sections[io, ic, 1]], inner_fraction
            )
        end
    end

    return interp_sections
end

"""
    approximate_streamlines(duct_inner, hub, xstream, rotor_x, rotor_r, rotor_id)

approximate streamlines based on a simple conservation of mass relative to the duct and center body geometry and relative radial position of anchor point (blade element leading or trailing edge location).
"""
function approximate_streamlines(duct_inner, hub, xstream, rotor_x, rotor_r, rotor_id)

    # - dimensions - #
    ns = length(rotor_x)

    # spline body geometry
    dsp = fm.Akima(duct_inner[:, 1], duct_inner[:, 2])
    hsp = fm.Akima(hub[:, 1], hub[:, 2])
    # initialize axial and radial coordinates of streamlins
    rstream = zeros(ns, N)

    for is in 1:ns
        sx = rotor_x[is]
        sr = rotor_r[is]

        # get reference area at this section
        ref_area = pi * (dsp(sx)^2 - hsp(sx)^2)

        # get the section area associated with this section
        section_area = pi * (sr^2 - hsp(sx)^2)
        # section_area = pi * (rgrid[1, 2:end] .^ 2 .- rgrid[1, 1:(end - 1)] .^ 2)

        # - loop through x locations from hub le to rotor section le
        for (ix, x) in enumerate(xstream[is, :])
            # get duct area of axial point of interest
            local_area = pi * (dsp(x)^2 - hsp(x)^2)

            # calclulate relative local section area
            local_section_area = local_area * section_area / ref_area

            # get radial position
            rstream[is, ix] = sqrt(local_section_area / pi + hsp(x)^2)
        end
    end

    return xstream, rstream
end

"""
    convert_to_mrt(sections; cutoff=4, conversion_factor=0.0254) #inch to meter

convert x,y,z coordinates to m, r*theta coordinates
also normalize, de-rotate, and return chord and stagger.

takes in sections, optionally cuts of closed blunt trailing edge
de-rotates, noramlizes, converts to m, r*theta
returns normalzied m coordinates, normlaized r*theta coordinates, stagger, and chord
"""
function convert_to_mrt(sections; cutoff=4, conversion_factor=1.0) #inch to meter
    # - get dimensions - #
    ns = length(sections[:, 1, 1])
    nc = length(sections[1, :, 1]) - cutoff * 2
    half = ceil(Int, nc / 2)
    mnorm = zeros(ns, nc)
    rthetanorm = zeros(ns, nc)
    stagger = zeros(ns)
    chord = zeros(ns)

    # - Loop through sections - #
    for is in 1:ns
        # get coordinates for this section and convert lengths if desired
        coordinates = sections[is, (cutoff + 1):(end - cutoff), :] .* conversion_factor

        ## -- r*theta coordinate conversion -- ##
        # get r values at each point using pythagorean theorem
        r = sqrt.(coordinates[:, 2] .^ 2 .+ coordinates[:, 3] .^ 2)

        # - use circle chord formula to get theta coordinate at each point. - #
        #p0 is at z=r and theta=zero for each axial location
        p0 = [coordinates[:, 1] zeros(nc) r .* ones(nc)]
        #a is the chord length from p0 to the coordinate
        a = [sign(coordinates[ic, 2]) * norm(coordinates[ic, :] - p0[ic, :]) for ic in 1:nc]
        # theta comes from the chord, angle, radius relation for circles
        theta = 2.0 .* asin.(a ./ (2.0 .* r))

        # r*theta coordinates
        rtheta = r .* theta

        # move so that leading edge is at zero
        rtheta .-= rtheta[half]

        ## -- m coordinate conversion -- ##
        # - split into halves - #
        # split into halves based on furthers forward point first
        # get minimum edge as index of minimum x location.
        minx, minidx = findmin(coordinates[:, 1])

        # split "halves"
        h1 = coordinates[1:minidx, :]
        h2 = coordinates[minidx:end, :]

        # reorder to work for splines
        rev1 = false
        rev2 = false
        if h1[1, 1] > h1[end, 1]
            reverse!(h1; dims=1)
            rev1 = true
        end
        if h2[1, 1] > h2[end, 1]
            reverse!(h2; dims=1)
            rev2 = true
        end

        # get lengths between coordinates, starting at zero for leading edge coordinate
        len(p1, p2) = norm(p2 - p1)
        h1len = zeros(length(h1[:, 1]))
        h2len = zeros(length(h2[:, 1]))
        h1m = [h1[:, 1] h1[:, 3]]
        h2m = [h2[:, 1] h2[:, 3]]
        h1len[2:end] = [len(h1m[ic + 1, :], h1m[ic, :]) for ic in 1:(length(h1[:, 1]) - 1)]
        h2len[2:end] = [len(h2m[ic + 1, :], h2m[ic, :]) for ic in 1:(length(h2[:, 1]) - 1)]

        # do cumulative sum starting with first coordinate at zero
        h1m = cumsum(h1len)
        h2m = cumsum(h2len)

        # re-reverse whichever half got reversed before.
        if rev1
            reverse!(h1m)
        end
        if rev2
            reverse!(h2m)
        end

        # put m coordinates back together, remembering the repeated leading edge point when you split things before.
        m = [h1m[1:(end - 1)]; h2m]

        # shift so that halfway point is at zero after getting m's.
        m .-= m[half]

        ## -- Stagger -- ##

        # Assemble for rotation
        mrt = [m rtheta]

        TE = (mrt[1, :] .+ mrt[end, :]) / 2.0
        chordvec = TE .- mrt[half, :]
        axvec = [1.0; 0.0]
        tempchord = norm(chordvec)
        stagger[is] = acos(dot(chordvec, axvec) / (tempchord))
        if chordvec[2] > 0.0
            stagger[is] *= -1.0
        end

        # - De-rotate Coordinates - #
        # get rotation matrix
        R = Rmatrix(stagger[is])

        #apply rotation
        rot = R * mrt'
        rotmrt = rot'

        ## -- Chord -- ##
        # set TE to be midpoint of 1st and last point
        # get chord lengths from LE and TE

        rotTE = (rotmrt[1, :] .+ rotmrt[end, :]) / 2.0
        chord[is] = rotTE[1] - rotmrt[half, 1]

        mnorm[is, :] = rotmrt[:, 1] / chord[is]
        rthetanorm[is, :] = rotmrt[:, 2] / chord[is]
    end

    return mnorm, rthetanorm, stagger, chord
end

function Rmatrix(angle)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end

"""
    write_mrt(m, rt, stagger, chord, fileprefix; filetype=".csv", savepath="")

writes m, r*theta coordinates to csv files along with a file conatining chord and stagger information
"""
function write_mrt(m, rt, stagger, chord, fileprefix; filetype=".csv", savepath="")

    # get sizes
    nc = size(m, 2)
    ns = length(chord)

    #write chord length and stagger angle
    chordstaggerfile = savepath * fileprefix * "_chord_stagger" * filetype
    g = open(chordstaggerfile, "w")
    if filetype == ".csv"
        write(g, "chord (m), stagger (rad)\n")
    else
        write(g, "#chord (m), stagger (rad)\n")
        write(g, "chord_stagger = [\n")
    end

    # loop through sections
    for is in 1:ns
        # - write chord and stagger - #
        if filetype == ".csv"
            write(g, "$(chord[is]), $(stagger[is])\n")
        else
            write(g, "$(chord[is]) $(stagger[is])\n")
        end

        # - write blade section coordinates - #
        #get file name for section
        coordinatefile = savepath * fileprefix * "_section_$(is)_m_rtheta" * filetype

        #open files
        f = open(coordinatefile, "w")
        # header
        if filetype == ".csv"
            write(f, "m, r-theta\n")
        else
            write(f, "#m, r-theta\n")
            write(f, "mr = [\n")
        end

        #loop through m and rtheta coordinates
        for ic in 1:nc
            if filetype == ".csv"
                write(f, "$(m[is,ic]), $(rt[is,ic])\n")
            else
                write(f, "$(m[is,ic]) $(rt[is,ic])\n")
            end
        end

        # close coordinate file
        if filetype != ".csv"
            write(f, "]")
        end
        close(f)
    end

    # close seciton file
    if filetype != ".csv"
        write(g, "]")
    end
    close(g)

    return nothing
end

function write_dtrotor(rs, chords, staggers, varname, filename)
    ns, nc = size(rs)
    f = open(filename, "w")
    write(f, "# $varname: radial position, chord, twist\n\n")
    write(f, varname * " = [\n")

    for is in 1:ns
        r = rs[is, :]
        rmid = r[ceil(Int, length(r) / 4)]
        twist = pi / 2.0 - staggers[is]
        write(f, "$rmid $(chords[is]) $twist\n")
    end
    write(f, "]")
    close(f)

    return nothing
end

######################################################################
#                                                                    #
#                         FUNCTIONS Not USED                         #
#                                                                    #
######################################################################
#= these were all written as first attempts. may be useful later, but not used for now

function transform3d(sections; cut=3)
    copysections = copy(sections)
    transformed_sections = sections[:, cut:(end - cut), :]

    ncoord = length(sections[1, :, 1])
    half = ceil(Int, ncoord / 2)

    ns = length(sections[:, 1, 1])
    for is in 1:ns
        trans = sections[is, half, :]
        transformed_sections[is, :, 1] .-= trans[1]
        transformed_sections[is, :, 2] .-= trans[2]
    end

    return transformed_sections
end

function rotated_surface_check(coordinates)

    # get leading edge as index of minimum x location.
    minx, minidx = findmin(coordinates[:, 1])

    # split "halves"
    h1 = coordinates[1:minidx, :]
    h2 = coordinates[minidx:end, :]

    # reorder to work for splines
    rev1 = false
    rev2 = false
    if h1[1, 1] > h1[end, 1]
        reverse!(h1; dims=1)
        rev1 = true
    end
    if h2[1, 1] > h2[end, 1]
        reverse!(h2; dims=1)
        rev2 = true
    end

    # spline "halves" based on x-location
    spy1 = fm.Akima(h1[:, 1], h1[:, 2])
    spz1 = fm.Akima(h1[:, 1], h1[:, 3])
    spy2 = fm.Akima(h2[:, 1], h2[:, 2])
    spz2 = fm.Akima(h2[:, 1], h2[:, 3])

    #sanity check
    x = (maximum(coordinates[:, 1]) + minimum(coordinates[:, 1])) / 2.0
    y1check = spy1(x)
    z1check = spz1(x)
    r1check = sqrt(y1check^2 + z1check^2)
    y2check = spy2(x)
    z2check = spz2(x)
    r2check = sqrt(y2check^2 + z2check^2)
    println("section $is cylinder check:", isapprox(r1check, r2check; atol=1e-5))
    if !isapprox(r1check, r2check; atol=1e-5)
        println("\tr1check = ", r1check)
        println("\tr2check = ", r2check)
    end

    return nothing
end

function extract3d_axial(sections, sample_r; finterp=fm.akima)
    nrext = length(sample_r)
    npts = size(sections, 2)
    extracts = zeros(nrext, npts, 2)

    for j in 1:npts
        extracts[:, j, 1] = finterp(sections[:, j, 3], sections[:, j, 1], sample_r)
        extracts[:, j, 2] = finterp(sections[:, j, 3], sections[:, j, 2], sample_r)
    end

    return extracts
end

function parameterize(vals, maxr, minr)
    ranger = maxr - minr
    return (vals .- minr) / ranger
end

function write3d(filename, coordinates, writetype)
    f = open(filename, writetype)

    # write(f, varname * " = [\n")
    for (i, c) in enumerate(eachrow(coordinates))
        write(f, "    $(c[1]) $(c[2]) $(c[3])\n")
    end
    # write(f, "]\n\n")
    close(f)

    return nothing
end

function write_normed_3d(extracts, sample_r, filename; cutoff=4, flip=false)
    half = ceil(Int, size(extracts, 2) / 2)
    ranger = maximum(sample_r) - minimum(sample_r)
    nondimr = round.(collect((sample_r .- minimum(sample_r)) / ranger), digits=2)
    chord = zeros(size(nondimr))
    stagger = zeros(size(nondimr))

    f = open(filename, "w")
    write(f, "Rtip = $(maximum(sample_r))\n")
    write(f, "Rhub = $(minimum(sample_r))\n\n")

    #get chord and stagger angle
    #and write radial position, chord, and stagger

    write(
        f,
        "#non-dimensional radius between hub and tip, chord normalized by tip radius, stagger in degrees\n",
    )
    write(f, "rcs = [\n")
    for i in 1:size(extracts, 1)
        chord[i], stagger[i] = calc_chord_stagger(extracts[i, :, :])
        write(f, "$(nondimr[i]) $(chord[i]) $(stagger[i])\n")
    end
    write(f, "]\n\n")
    close(f)

    #write
    for i in 1:size(extracts, 1)
        varname = "r$(round(Int,nondimr[i]*100))"
        coordinates = extracts[i, :, :]
        R = Rmatrix(stagger[i])
        S = chord[i]
        T = coordinates[half, :]
        coordinates[:, 1] .-= T[1]
        coordinates[:, 2] .-= T[2]
        coordinates ./= S
        coordinates = R * coordinates'
        coordinates = coordinates[:, cutoff:(end - cutoff)]'
        if flip
            coordinates[:, 2] .*= -1.0
        end
        write2d(filename, varname, coordinates, "a")
    end

    return nothing
end

function calc_chord_stagger(coordinates)
    nc = size(coordinates, 1)
    half = ceil(Int, nc / 2)
    chordvec = coordinates[1, :] .- coordinates[half, :]
    chord = norm(chordvec)
    axvec = [1.0; 0.0]
    stagger = acosd(dot(chordvec, axvec) / (chord))
    if chordvec[2] > 0.0
        stagger *= -1.0
    end

    return chord, stagger
end
=#
