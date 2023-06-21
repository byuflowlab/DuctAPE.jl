"""
"""
function write_openvsp_bem(; savepath, filename, header, data)
    write_openvsp_bem(;
        # File Info
        savepath,
        filename,
        # Header
        rotor_name=header.rotor_name,
        num_sections=header.num_sections,
        num_blades=header.num_blades,
        diameter=header.diameter,
        beta=header.beta,
        feather=header.feather,
        precone=header.precone,
        center=header.center,
        normal=header.normal,
        # Blade Data
        radius=data.radius,
        chord=data.chord,
        twist=data.twist,
        rake=data.rake,
        skew=data.skew,
        sweep=data.sweep,
        thickness=data.thickness,
        cli=data.cli,
        axial=data.axial,
        tangential=data.tangential,
    )
    return nothing
end

"""
    write_openvsp_bem(; kwargs)

Write an OpenVSP .bem file for visualizing rotor geometry.

**Keyword Arguments:**
## File Info
- `savepath::String` : save path
- `filename::String` : file name, will allways be appended by ".bem"
## Header Info
- `header::NamedTuple` : named tuple containing the header info (items below).
--OR--
- `num_sections::Int` :
- `num_blades::Int` :
- `diameter::Float64` :
- `beta::Float64` : (degrees)
- `feather::Float64` : (degrees)
- `precone::Float64` : (degrees)
- `center::Vector{Float64}` :
- `normal::Vector{Float64}` :
## Blade Data
- `data::NamedTuple` : named tuple containing the blade data (items below).
--OR--
- `radius::Vector{Float64}` :
- `chord::Vector{Float64}` :
- `twist::Vector{Float64}` :
- `rake::Vector{Float64}` :
- `skew::Vector{Float64}` :
- `sweep::Vector{Float64}` :
- `thickness::Vector{Float64}` :
- `cli::Vector{Float64}` :
- `axial::Vector{Float64}` :
- `tangential::Vector{Float64}` :
"""
function write_openvsp_bem(;
    # File Info
    savepath="",
    filename="test_prop",
    # Header
    rotor_name="Test Prop",
    num_sections=41,
    num_blades=3,
    diameter=5.0,
    beta=20.0,
    feather=0.0,
    precone=0.0,
    center=[0.0; 0.0; 0.0],
    normal=[-1.0; 0.0; 0.0],
    # Blade Data
    radius=range(0.0, 1.0, num_sections),
    chord=ones(num_sections),
    twist=zeros(num_sections),
    rake=zeros(num_sections),
    skew=zeros(num_sections),
    sweep=zeros(num_sections),
    thickness=0.05 * ones(num_sections),
    cli=zeros(num_sections),
    axial=zeros(num_sections),
    tangential=zeros(num_sections),
)

    # check that the inputs are the right lengths
    @assert length(center) == length(normal) == 3
    @assert allequal(
        [
            num_sections
            length(radius)
            length(chord)
            length(twist)
            length(rake)
            length(skew)
            length(sweep)
            length(thickness)
            length(cli)
            length(axial)
            length(tangential)
        ],
    )

    # open file
    f = open(savepath * filename * ".bem", "w")

    # Write Header
    write(f, rotor_name * "\n")
    write(f, "Num_Sections: $(num_sections)\n")
    write(f, "Num_Blade: $(num_blades)\n")
    write(f, "Diameter: $(diameter)\n")
    write(f, "Beta 3/4 (deg): $(beta)\n")
    write(f, "Feather (deg): $(feather)\n")
    write(f, "Pre_Cone (deg): $(precone)\n")
    write(f, "Center: $(center[1]), $(center[2]), $(center[3])\n")
    write(f, "Normal: $(normal[1]), $(normal[2]), $(normal[3])\n")
    write(f, "\n")

    # write column headings for data
    write(
        f,
        "Radius/R, Chord/R, Twist (deg), Rake/R, Skew/R, Sweep, t/c, CLi, Axial, Tangential\n",
    )

    # Loop through data
    for i in 1:num_sections
        write(
            f,
            "$(radius[i]), $(chord[i]), $(twist[i]), $(rake[i]), $(skew[i]), $(sweep[i]), $(thickness[i]), $(cli[i]), $(axial[i]), $(tangential[i]),\n",
        )
    end

    # close file
    close(f)

    return nothing
end
