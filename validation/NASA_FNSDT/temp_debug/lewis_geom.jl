using DuctAPE
const dt = DuctAPE

project_dir = dirname(dirname(dirname(dirname(@__FILE__))))
if project_dir == ""
    project_dir = "."
end

savepath = project_dir * "/examples/dfdc_partial/dfdc_gen_files/"

# - Duct - #
include(project_dir * "/test/data/naca_662-015.jl")

ductcoordinates = [x_duct r_duct]
drefcoordinates = dt.repanel_airfoil(ductcoordinates; N=80, normalize=false)
#check that leading edge is at zero
# le, _ = findmin(drefcoordinates[:,1])
# drefcoordinates[:,1] .-= le

# f = open(savepath * "lewis_refined_duct.jl", "w")
# write(f, "duct_coordinates = [\n")
# for xr in eachrow(refcoordinates)
#     write(f, "$(xr[1]) $(xr[2])\n")
# end
# write(f, "]")
# close(f)

# f = open(savepath * "lewis_inf_duct.jl", "w")
# write(f, "duct_coordinates = [\n")
# for xr in eachrow(refcoordinates)
#     write(f, "$(xr[1]) $(1e1+xr[2])\n")
# end
# write(f, "]")
# close(f)

# - Hub - #
include(project_dir * "/test/data/bodyofrevolutioncoords.jl")

hubcoordinates = [x_hub r_hub]
hrefcoordinates = dt.repanel_revolution(hubcoordinates; N=40, normalize=false)

# f = open(savepath * "lewis_refined_hub.jl", "w")
# write(f, "hub_coordinates = [\n")
# for xr in eachrow(refcoordinates)
#     write(f, "$(xr[1]) $(xr[2])\n")
# end
# write(f, "]")
# close(f)

# f = open(savepath * "lewis_inf_hub.jl", "w")
# write(f, "hub_coordinates = [\n")
# for xr in eachrow(refcoordinates)
#     write(f, "$(xr[1]) 1.0E-8\n")
# end
# write(f, "]")
# close(f)

for i in 1:10
    run(`cp $(savepath)fl $(savepath)fh$(i)`)
    f = open(savepath * "fh$(i)", "a")
    for xr in eachrow(hrefcoordinates)
        write(f, "$(xr[1]) $(xr[2])\n")
    end
    write(f, "999.0 999.0\n")

    ddimcoordinates = [drefcoordinates[:, 1] drefcoordinates[:, 2] .+ (i)]
    for xr in eachrow(ddimcoordinates)
        write(f, "$(xr[1]) $(xr[2])\n")
    end
    write(f, "ENDGEOM")
    close(f)
end

for i in 1:10
    run(`cp $(savepath)fl $(savepath)fd$(i)`)
    f = open(savepath * "fd$(i)", "a")

    hdimcoordinates = [hrefcoordinates[:, 1]  hrefcoordinates[:, 2] / (2 * i)]
    for xr in eachrow(hdimcoordinates)
        write(f, "$(xr[1]) $(xr[2])\n")
    end
    write(f, "999.0 999.0\n")

    for xr in eachrow(drefcoordinates)
        write(f, "$(xr[1]) $(xr[2])\n")
    end
    write(f, "ENDGEOM")
    close(f)
end

run(`cp $(savepath)fl $(savepath)fl0`)
f = open(savepath * "fl0", "a")
for xr in eachrow(hrefcoordinates)
    write(f, "$(xr[1]) $(xr[2])\n")
end
write(f, "999.0 999.0\n")
for xr in eachrow(drefcoordinates)
    write(f, "$(xr[1]) $(xr[2])\n")
end
write(f, "ENDGEOM")
close(f)
