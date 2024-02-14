# compare my polar with dfdc's polar to make sure they're in the same ballpark (probably will need to just use dfdc's polar)


project_dir = dirname(dirname(dirname(@__FILE__)))
if project_dir == ""
    project_dir = "."
end

include(project_dir * "/plots_default.jl")

#load dfdc lift/drag polar
include(project_dir*"/dev_debug_archive/dfdc_comp/DFDC_POLAR.jl")
#load my lift/drag polar
include(project_dir*"/dev_debug_archive/dfdc_comp/mypolar.jl")



pcl = plot(xlabel="Angle of Attack (deg)", ylabel=L"c_l")
pcd = plot(xlabel="Angle of Attack (deg)", ylabel=L"c_d")

plot!(pcl,dfdc_polar[:,1]*180.0/pi, dfdc_polar[:,2], label = "DFDC Polar")
plot!(pcl,mypolar[:,1]*180.0/pi, mypolar[:,2], label = "My Polar")

plot!(pcd,dfdc_polar[:,1]*180.0/pi, dfdc_polar[:,3], label = "DFDC Polar")
plot!(pcd,mypolar[:,1]*180.0/pi, mypolar[:,3], label = "My Polar")

savefig(pcl, project_dir*"/dev_debug_archive/dfdc_comp/clcomp.pdf")
savefig(pcd, project_dir*"/dev_debug_archive/dfdc_comp/cdcomp.pdf")
