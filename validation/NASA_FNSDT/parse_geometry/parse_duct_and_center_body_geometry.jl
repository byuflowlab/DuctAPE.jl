include("geometry_parsing_functions.jl") # includes inch to meter conversion, i2m

datapath = "geometry_parsing/NASA_geometry_files/"
savepath = "geometry_parsing/parsed_geometry_files/"

# - Initialize Plots - #
p2 = plot(; apsectratio=1, xlabel="x (m)", ylabel="r (m)", grid=false)

#---------------------------------#
#       Duct Outer Surface        #
#---------------------------------#
file = datapath * "nacelle.xr" #in inches

# parse file into x,r coordinates
x, r = parse2d(file)

# plot
plot!(p2, x * i2m, r * i2m; aspectratio=1, label="nacelle")

ductout = [x r] .* i2m

#---------------------------------#
#       Duct Inner Surface        #
#---------------------------------#
file = datapath * "casing.xr" #in inches

# parse file into x,r coordinates
x, r = parse2d(file)

# plot
plot!(p2, x .* i2m, r .* i2m; aspectratio=1, label="casing")

# - DuctTAPE needs these reversed - #
ductin = reverse([x r]; dims=1) .* i2m

# put the duct coordinates together
duct = [ductin[1:(end - 1), :]; ductout]

#---------------------------------#
#          Center Body            #
#---------------------------------#
file = datapath * "centerbody.xr" #in inches

# parse file into x,r coordinates
x, r = parse2d(file)

# plot
plot!(p2, x .* i2m, r .* i2m; aspectratio=1, label="centerbody")

hub = [x r] .* i2m

# - Close TE Geoemtry - #
hubclosed = copy(hub)
dx = x[end] - x[end - 1]
dr = r[end] - r[end - 1]
nr = r[end] / abs(dr)
ter = range(r[end], 0.0; step=dr)
tex = range(x[end], x[end] + dx * nr; step=dx)
hubclosed = [hubclosed[1:end-1,:]; [tex ter]*i2m]

#---------------------------------#
#              SAVE               #
#---------------------------------#
savefig(savepath * "../figures/extracted_duct_and_hub_geometry.pdf")

# Write data
write2d(
    savepath * "extracted_duct_and_centerbody_geometry_in_julia_arrays.jl",
    "duct_coordinates",
    duct,
    "w",
)
write2d(
    savepath * "extracted_duct_and_centerbody_geometry_in_julia_arrays.jl",
    "hub_coordinates_open_TE",
    hub,
    "a",
)
write2d(
    savepath * "extracted_duct_and_centerbody_geometry_in_julia_arrays.jl",
    "hub_coordinates_closed_TE",
    hubclosed,
    "a",
)
