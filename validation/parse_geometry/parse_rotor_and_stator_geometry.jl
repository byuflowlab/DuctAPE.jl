# Get streamlines first, which will also get the interpolated geometries
include("approximate_streamlines_through_rotor_and_stator.jl")

#---------------------------------#
#              Rotor              #
#---------------------------------#
m_rotor, rt_rotor, stagger_rotor, chord_rotor = convert_to_mrt(
    interp_rotor_sections; cutoff=4
)
write_mrt(m_rotor, -rt_rotor, stagger_rotor, chord_rotor, savepath * "rotor")
write_mrt(
    m_rotor,
    -rt_rotor,
    stagger_rotor,
    chord_rotor,
    savepath * "../jls/rotor";
    filetype=".jl",
)
write_dtrotor(
    sqrt.(interp_rotor_sections[:, :, 2] .^ 2 .+ interp_rotor_sections[:, :, 3] .^ 2),
    chord_rotor,
    stagger_rotor,
    "rotor_rct",
    "../validation_scripts/rotor_geometry.jl",
)

plot(; legend=:outertopright, aspectratio=1)
for i in 1:24:length(m_rotor[:, 1])
    plot!(
        m_rotor[i, :] * chord_rotor[i] ./ i2m,
        -rt_rotor[i, :] * chord_rotor[i] ./ i2m;
        label="$i",
    )
end
savefig(savepath * "../figures/rotorsecitonsanitycheck.pdf")

#---------------------------------#
#              Stator             #
#---------------------------------#
m_stator, rt_stator, stagger_stator, chord_stator = convert_to_mrt(
    interp_stator_sections; cutoff=5
)
write_mrt(m_stator, rt_stator, -stagger_stator, chord_stator, savepath * "stator")
write_mrt(
    m_stator,
    rt_stator,
    -stagger_stator,
    chord_stator,
    savepath * "../jls/stator";
    filetype=".jl",
)
write_dtrotor(
    sqrt.(interp_stator_sections[:, :, 2] .^ 2 .+ interp_stator_sections[:, :, 3] .^ 2),
    chord_stator,
    -stagger_stator,
    "stator_rct",
    "../validation_scripts/stator_geometry.jl",
)

plot(; legend=:outertopright, aspectratio=1)
for i in 1:24:length(m_stator[:, 1])
    plot!(
        m_stator[i, :] * chord_stator[i] ./ i2m,
        rt_stator[i, :] * chord_stator[i] ./ i2m;
        label="$i",
    )
end
savefig(savepath * "../figures/statorsecitonsanitycheck.pdf")

plot3d(; aspectratio=:equal)
for is in 1:length(interp_rotor_sections[:, 1, 1])
    plot3d!(
        interp_rotor_sections[is, :, 1],
        interp_rotor_sections[is, :, 2],
        interp_rotor_sections[is, :, 3];
        color=:blue,
        label="",
    )
end
for is in 1:length(rotor_sections[:, 1, 1])
    plot3d!(
        rotor_sections[is, :, 1],
        rotor_sections[is, :, 2],
        rotor_sections[is, :, 3];
        color=:gray,
        label="",
    )
end
savefig(savepath * "../figures/3dsanitycheck.pdf")
