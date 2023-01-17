
plot(; aspectratio=:equal)

plot!(dxo, dro; color=mycolors[1], label="Duct Surface")
plot!(dxi, dri; color=mycolors[1], label="")

plot!(hx, hr; color=mycolors[2], label="Hub Surface")

#Plot rotor quarter chord
plot!(
    [rotor_c4_pos; rotor_c4_pos],
    [blade_elements.radial_positions[1]; blade_elements.radial_positions[end]];
    linewidth=1.5,
    linestyle=:dash,
    color=mycolors[3],
    label="Rotor Quarter Chord Line",
)

savefig("examples/lilium_like_geometry_side.pdf")

#myyticks = [0.0; Rhub; P.Rtip; maximum(ductrmesh)]
#yticklabels = [@sprintf("%1.2f", x) for x in myyticks]
#yticklabels[1] = "0"
#plot!(;
#    xticks=(
#        [
#            0.0
#            P.xlocation
#            maximum(ductxmesh);
#            maximum(hubx)
#        ],
#        [
#            "0"
#            @sprintf("%1.2f", P.xlocation)
#            @sprintf("%1.2f", maximum(ductxmesh))
#            @sprintf("%1.2f", maximum(hubx))
#        ],
#    ),
#    yticks=(myyticks, yticklabels),
#    tick_direction=:out,
#    grid=:true,
#)

##plot section cut line+arrows
#plot!(
#    [rotor_le[1]; rotor_te[1]],
#    [Rhub + (Rtip - Rhub) * 0.01; Rhub + (Rtip - Rhub) * 0.01];
#    linestyle=:dash,
#    color=:black,
#    linewidth=1.5,
#    label="",
#)
#plot!(
#    [
#        rotor_le[1]
#        rotor_le[1] - abs(minimum(hubxsym) - rotor_le[1]) * 0.1
#        rotor_le[1] - abs(minimum(hubxsym) - rotor_le[1]) * 0.1
#    ],
#    [Rhub + (Rtip - Rhub) * 0.01; Rhub + (Rtip - Rhub) * 0.01; Rhub + (Rtip - Rhub) * 0.15];
#    color=:black,
#    arrow=true,
#    linewidth=1.5,
#    label="",
#)
#plot!(
#    [
#        rotor_te[1]
#        rotor_te[1] + (maximum(hubxsym) - rotor_te[1]) * 0.1
#        rotor_te[1] + (maximum(hubxsym) - rotor_te[1]) * 0.1
#    ],
#    [Rhub + (Rtip - Rhub) * 0.01; Rhub + (Rtip - Rhub) * 0.01; Rhub + (Rtip - Rhub) * 0.15];
#    color=:black,
#    arrow=true,
#    linewidth=1.5,
#    label="",
#)

# annotate!(
#     rotor_te[1] + (maximum(hubxsym) - rotor_te[1]) * 0.25,
#     Rhub + (Rtip - Rhub) * 0.1,
#     text(L"\mathrm{\cref{fig:initialblade}}", 10),
# )
# annotate!(
#     rotor_le[1] - abs(minimum(hubxsym) - rotor_le[1]) * 0.3,
#     Rhub + (Rtip - Rhub) * 0.1,
#     text(L"\mathrm{\cref{fig:initialblade}}", 10),
# )

#######################################################
######     Starting Point Geometry FRONT VIEW     #####
#######################################################

#th = range(0.0, 2 * pi; length=360 * 2)
#_, leidx = findmin(ductxmesh)

#plot(;
#    aspectratio=:equal,
#    tick_direction=:out,
#    xticks=(
#        [0.0; maximum(hubr); Rtip],
#        [
#            "0"
#            @sprintf("%1.2f", maximum(hubr))
#            @sprintf("%1.2f", Rtip)
#        ],
#    ),
#    yticks=(
#        [0.0; maximum(hubr); minimum(ductrmesh); ductrmesh[leidx]; maximum(ductrmesh)],
#        [
#            "0"
#            @sprintf("%1.2f", maximum(hubr))
#            @sprintf("%1.2f", minimum(ductrmesh))
#            @sprintf("%1.2f", ductrmesh[leidx])
#            @sprintf("%1.2f", maximum(ductrmesh))
#        ],
#    ),
#    grid=true,
#)

##Get rotor LE and TE
#rotor_le = (section_chords .* P.Rtip .* 0.25) .* cosd.(section_twists)
#rotor_te = (section_chords .* P.Rtip .* 0.75) .* cosd.(section_twists)

##assemble blade
#qc = [dim_rad_stash zeros(length(dim_rad_stash))]
#le = [dim_rad_stash rotor_le]
#te = [dim_rad_stash -rotor_te]

##plot blades
#for i in 1:num_blades
#    angle = 2 * pi / num_blades * i
#    R(angle) = [cos(angle) -sin(angle); sin(angle) cos(angle)]
#    qcrot = R(angle) * qc'
#    leangle = 2.0 * asin(rotor_le[end] / (2 * Rtip))
#    teangle = 2.0 * asin(rotor_te[end] / (2 * Rtip))
#    lerot = R(angle + leangle) * qc'
#    terot = R(angle - teangle) * qc'

#    if i == 1
#        #quarter chord
#        plot!(
#            qcrot[1, :],
#            qcrot[2, :];
#            linestyle=:dash,
#            color=mycolors[3],
#            label="Rotor Quarter Chord Lines",
#        )

#        #le
#        plot!(
#            lerot[1, :],
#            lerot[2, :];
#            color=mycolors[3],
#            label="Rotor Blade Profiles (including twist)",
#        )

#        #te
#        plot!(terot[1, :], terot[2, :]; color=mycolors[3], label="")
#    else
#        plot!(qcrot[1, :], qcrot[2, :]; linestyle=:dash, color=mycolors[3], label="")
#        plot!(lerot[1, :], lerot[2, :]; color=mycolors[3], label="")
#        plot!(terot[1, :], terot[2, :]; color=mycolors[3], label="")
#    end
#end

##plot circle for duct
#x = maximum(ductrmesh) .* cos.(th)
#y = maximum(ductrmesh) .* sin.(th)
#plot!(x, y; color=mycolors[1], label="Duct Profile") #fill=(0, 1.0, mycolors[2]))
#x = minimum(ductrmesh) .* cos.(th)
#y = minimum(ductrmesh) .* sin.(th)
#xtip = Rtip .* cos.(th)
#ytip = Rtip .* sin.(th)
##nullify where blades are
#for i in 1:num_blades
#    for j in 1:length(th)
#        Rduct = minimum(ductrmesh)
#        cle = rotor_le[end] * Rduct / Rtip
#        cte = rotor_te[end] * Rduct / Rtip
#        angle = 2.0 * pi / num_blades * i
#        nangle = 2.0 * pi / num_blades * (i + 1)
#        leangle = 2.0 * asin(cle / (2 * Rduct))
#        teangle = 2.0 * asin(cte / (2 * Rduct))
#        if (th[j] > angle - teangle && th[j] < angle + leangle)
#            y[j] = NaN
#        end
#        if (th[j] > 0.0 && th[j] < leangle)
#            y[j] = NaN
#        end
#        if th[j] < nangle - teangle && th[j] > angle + leangle
#            ytip[j] = NaN
#        end
#        if th[j] < 2.0 * pi / num_blades - teangle && th[j] > leangle
#            ytip[j] = NaN
#        end
#    end
#end
#plot!(x, y; color=mycolors[1], label="")
#plot!(xtip, ytip; color=mycolors[3], label="")

##plot duct leading edge
#x = ductrmesh[leidx] .* cos.(th)
#y = ductrmesh[leidx] .* sin.(th)
#plot!(x, y; color=mycolors[1], linestyle=:dot, label="Duct Leading Edge")

##plot circle for hub
#x = Rhub .* cos.(th)
#y = Rhub .* sin.(th)

#plot!(x, y; z_order=:front, color=mycolors[2], label="Hub Profile") #fill=(0, 1.0, mycolors[2]))

#savefig("./latex/figures/initial_geometry_front.tikz")

