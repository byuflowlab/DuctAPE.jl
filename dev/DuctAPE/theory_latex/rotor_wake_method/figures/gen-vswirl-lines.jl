hc = [
    0.0 0.0
    0.000220493763 0.00327158696
    0.000818265835 0.0062371837
    0.00176513905 0.00907455198
    0.00305903051 0.0118233347
    0.00469737826 0.0144944247
    0.00667320006 0.01708664
    0.00898412429 0.0195938926
    0.0116183488 0.0220112894
    0.014565343 0.0243293047
    0.0178080108 0.0265385918
    0.0213311911 0.0286334325
    0.0251127146 0.0306049269
    0.0291345045 0.0324492492
    0.0333757736 0.0341628119
    0.0378172994 0.0357434638
    0.042441953 0.0371914469
    0.0472304933 0.0385068245
    0.0521632098 0.0396893248
    0.0572345369 0.0407429487
    0.0624375977 0.0416737311
    0.0677433684 0.0424805768
    0.0731483325 0.0431673825
    0.0786491632 0.0437394157
    0.0842350498 0.044200331
    0.0898977965 0.044553142
    0.0956536308 0.0448064096
    0.101440892 0.0449600294
    0.107219525 0.0450017527
    0.11342375 0.0449719056
    0.120000005 0.0449555703
    0.126668498 0.0449563377
    0.133435443 0.0449485704
    0.140258059 0.0449409299
    0.147098973 0.04493507
    0.153945565 0.0449282154
    0.160799354 0.0449208207
    0.167654142 0.0449139029
    0.174511507 0.0449073091
    0.181371018 0.0449007712
    0.188229889 0.0448941104
    0.195088565 0.0448872671
    0.20194827 0.0448803529
    0.208807036 0.0448735543
    0.21566385 0.0448670611
    0.22251825 0.0448606126
    0.229360789 0.0448528677
    0.236195073 0.0448435619
    0.242983848 0.0448389612
    0.249588013 0.0448331758
    0.25564155 0.0447736457
    0.26116699 0.0446058325
    0.266400695 0.0443055779
    0.271526247 0.0438646637
    0.276578337 0.0432802737
    0.281544089 0.0425501764
    0.286387861 0.0416728891
    0.291097105 0.040644031
    0.295679897 0.0394582525
    0.300161272 0.0381119996
    0.304465979 0.0366399921
    0.30637899 0.0359279998
]

using FLOWMath
x = range(0, 1, 500)

scale = maximum(hc[:, 1])

hc1 = hc ./ scale
hc1[:, 2] ./= 1.5
hc1finer = FLOWMath.akima(hc1[:, 1], hc1[:, 2], x)
hc1fine = [x hc1finer]

hc5 = hc1fine .* 5.0

f = open("scaled_dfdc_hub_coordinates.dat", "w")
for (x, z) in zip(eachrow(hc5[:, 1]), eachrow(hc5[:, 2]))
    write(f, "$(x[1]) $(z[1])\n")
end
write(f, "$(hc5[end,1]) 0.0\n")
close(f)

using Plots
using DuctAPE
const dt = DuctAPE

plot(hc5[:, 1], hc5[:, 2]; aspectratio=1, label="")

b1 = findfirst(x -> x > 1.3, x * 5)
lower_wall = [hc5[b1:end, :]' [6.0; hc5[end, 2]]]
upper_wall = copy(lower_wall)
upper_wall[2, :] .= 5.0

rshift = sort([
    [0.0; 0.22; 0.44] .+ 0.04
    [0.0; 0.22; 0.44] .+ 0.04 * 2
])

grid = similar(lower_wall, 2, size(lower_wall, 2), length(rshift) + 2) .= 0.0
grid[1, :, 1:end] .= lower_wall[1, :]
grid[2, :, 1] = lower_wall[2, :]
grid[2, :, end] = upper_wall[2, :]
for (i, r) in enumerate(rshift)
    grid[2, :, i + 1] = lower_wall[2, :] .+ r
end

dt.relax_grid!(grid)
zs = grid[1, :, 1]
rs = grid[2, :, 2:(end - 1)]

s = [0.0; 0.22; 0.44] .+ 0.04
g = 2
v1 = 1.35
v2 = 1.42
v3 = 2.0
v4 = 2.05
v5 = 2.5

for (i, r) in enumerate(eachcol(rs))
    if i%2==0
        b = findfirst(x -> x > 2.52, zs)
        vi = [zs r][(b + g * i):end, :]
    else
        vi = [zs r][(1 + g * i):end, :]
    end
    vi = [[vi[1, 1] v1]; vi]
    plot!(vi[:, 1], vi[:, 2]; label="", color=2)

    f = open("swirl-velocity-horseshoe$(i).dat", "w")
    for (x, z) in zip(eachrow(vi[:, 1]), eachrow(vi[:, 2]))
        write(f, "$(x[1]) $(z[1])\n")
    end
    close(f)
    f = open("swirl-velocity-vert1$(i).dat", "w")
    write(f, "$(vi[1,1]) $(v2)\n")
    write(f, "$(vi[1,1]) $(v3)\n")
    close(f)
    f = open("swirl-velocity-vert2$(i).dat", "w")
    write(f, "$(vi[1,1]) $(v4)\n")
    write(f, "$(vi[1,1]) $(v5)\n")
    close(f)

end
plot!()

#TODO: need to get another line matching the hub profile but futher up to be the dotted line showing the alingment of the contours to integrate over.  need to chop it at the starting horseshoe, and then at the point you want to put the second contour, maybe at 4ish
s = [1.5 - maximum(hc5[:, 2])]
b1 = findfirst(x -> x > 1.05, x * 5)
b2 = findfirst(x -> x > 4.0, x * 5)
for i in 1
    vi = copy(hc5[b1:b2, :])
    vi[:, 2] .+= s[i]
    plot!(vi[:, 1], vi[:, 2]; label="", color=4, linestyle=:dash)
    f = open("swirl-velocity-radialpos.dat", "w")
    for (x, z) in zip(eachrow(vi[:, 1]), eachrow(vi[:, 2]))
        write(f, "$(x[1]) $(z[1])\n")
    end
    close(f)
end
plot!()

