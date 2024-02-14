project_dir = dirname(dirname(@__FILE__))
if project_dir == ""
    project_dir = "."
end

using FLOWMath
include(project_dir * "/plots_default.jl")
include(project_dir * "/dev_debug_archive/simple_airfoil.jl")

# read in naca 4412 polar file
a = Float64[]
l = Float64[]
d = Float64[]

open(project_dir * "/test/data/naca4412.dat") do f
    # skip header
    info = readline(f)
    Re = parse(Float64, readline(f))
    Mach = parse(Float64, readline(f))
    for line in eachline(f)
        parts = split(line)
        push!(a, parse(Float64, parts[1]))
        push!(l, parse(Float64, parts[2]))
        push!(d, parse(Float64, parts[3]))
    end
end
# TODO

#spline data super fine
alphafine = range(-5.0, 20.0; step=0.01)
clfine = FLOWMath.akima(a * 180.0 / pi, l, alphafine)
cdfine = FLOWMath.akima(a * 180.0 / pi, d, alphafine)

## -- Extract Lift Parameters -- ##

clmax, clmaxidx = findmax(clfine)
clmin = minimum(clfine)
alpha0idx = findfirst(x -> x >= 0.0, clfine)
alpha0 = alphafine[alpha0idx]
dclda =
    (clfine[alpha0idx + 10] - clfine[alpha0idx]) /
    (alphafine[alpha0idx + 10] - alphafine[alpha0idx])

dclda_stall =
    (clfine[clmaxidx + 10] - clfine[clmaxidx]) /
    (alphafine[clmaxidx + 10] - alphafine[clmaxidx])

blend_hardness = 5

lift_params = (; clmax, clmin, alpha0, dclda, dclda_stall, blend_hardness)

cdo, cdoidx = findmin(cdfine)
clo = clfine[cdoidx]
b = 4e-3 #arbiratry guess, doesn't really matter probably
f = -0.5 # somewhere between -0.3 and -0.5 for Re<100k
Re_ref = 5e5
Re = 5e5

drag_params = (; cdo, clo, b, Re_ref, f)

generate_airfoil_data(
    lift_params, drag_params, Re; filename=project_dir * "dev_debug_archive/simple_naca4412.dat"
)

# check that things look reasonable...

# Read in simplified data
as = []
ls = []
ds = []
open(project_dir * "/dev_debug_archive/simple_naca4412.dat") do f
    for line in eachline(f)
        parts = split(line)
        push!(as, parse(Float64, parts[1]))
        push!(ls, parse(Float64, parts[2]))
        push!(ds, parse(Float64, parts[3]))
    end
end

# plot them both

pcl = plot(; xlabel="angle of attack", ylabel="cl")
plot!(pcl, a, l; label="NACA4412")
plot!(pcl, as, ls; label="Simplified")
savefig(project_dir * "simplepolarcomp-cl.pdf")

pcd = plot(; xlabel="angle of attack", ylabel="cl")
plot!(pcd, a, d; label="NACA4412")
plot!(pcd, as, ds; label="Simplified")
savefig(project_dir * "simplepolarcomp-cd.pdf")
