project_dir = dirname(dirname(dirname(dirname(@__FILE__))))
if project_dir == ""
    project_dir = "."
end

using Infiltrator
using FLOWMath
using Xfoil
include(project_dir * "/plots_default.jl")
include(project_dir * "/test/data/APC_5x5.75_geom.jl")

function estimate_re(RPM, Vinf, Rtip, chords, radial_positions; rho=1.225, mu=1.81e-5)

    #rotation rate in rev per sec
    n = RPM / 60.0
    #rad/s
    omega = n * 2.0 * pi
    #anglular velocity m/s
    vtheta = omega * radial_positions

    #total velocity
    vtot = sqrt.(vtheta .^ 2 .+ Vinf .^ 2)

    re = rho * vtot .* chords / mu

    #rotor diameter
    D = 2.0 * Rtip

    return re, Vinf ./ (n * D)
end

RPM = 4500
Vinf = range(2.0, 20.0)
rerange = zeros(length(radial_positions), length(Vinf))
J = zeros(length(Vinf))
for i in 1:length(Vinf)
    rerange[:, i], J[i] = estimate_re(RPM, Vinf[i], Rtip, chords, radial_positions)
end

function generate_inviscid_polars(geometry_files, alphas=range(-15.0, 25.0; step=0.5))
    for i in 1:length(geometry_files)
        println("Aifoil: $(geometry_files[i])")
        #read in geometry file
        x = []
        z = []
        open(geometry_files[i] * ".coord") do f
            for line in eachline(f)
                parts = split(line)
                push!(x, parse(Float64, parts[1]))
                push!(z, parse(Float64, parts[2]))
            end
        end

        #run inviscid xfoil
        println("Running Inviscid")
        cl, _ = Xfoil.alpha_sweep(reverse(x), reverse(z), alphas)

        println("Writing Lift Polar")
        # open write file
        fw = open(geometry_files[i] * "_Lift.polar", "w")
        # place header
        write(fw, geometry_files[i] * " Inviscid\n")

        #write data if converged
        for p in 1:length(cl)
            #check convergence
            write(fw, "$(alphas[p]) $(cl[p])\n")
        end

        close(fw)
    end
    return nothing
end

generate_inviscid_polars(
    project_dir * "/test/data/APC-5x5.75_airfoils/" .* saved_files,
    range(-15.0, 25.0; step=0.5),
)

function generate_viscous_polars(
    geometry_files, rerange, alphas=range(-15.0, 25.0; step=0.5)
)
    for i in 1:length(geometry_files)
        println("Aifoil: $(geometry_files[i])")
        #read in geometry file
        x = []
        z = []
        open(geometry_files[i] * ".coord") do f
            for line in eachline(f)
                parts = split(line)
                push!(x, parse(Float64, parts[1]))
                push!(z, parse(Float64, parts[2]))
            end
        end

        #initialize polar values
        nr = length(rerange[1, :])

        #run for range of alphas and res
        for r in 1:5:nr
            println("running for Re = ", rerange[i, r])

            #run xfoil
            cl, cd, cdp, cm, conv = Xfoil.alpha_sweep(
                reverse(x),
                reverse(z),
                alphas,
                rerange[i, r];
                reinit=true,
                percussive_maintenance=true,
            )

            println("Writing Polar")
            # open write file
            fw = open(geometry_files[i] * "_Re_$(round(Int,rerange[i,r])).polar", "w")
            # place header
            write(fw, geometry_files[i] * " RAW\n")
            write(fw, "$(rerange[i,r])\n")
            write(fw, "0.0\n")

            #write data if converged
            for p in 1:length(cl)
                #check convergence
                if conv[p]
                    write(fw, "$(alphas[p]*pi/180.0) $(cl[p]) $(cd[p])\n")
                end
            end

            close(fw)
        end
    end
    return nothing
end

generate_viscous_polars(
    project_dir * "/test/data/APC-5x5.75_airfoils/" .* saved_files,
    rerange,
    range(-15.0, 25.0; step=0.5),
)
