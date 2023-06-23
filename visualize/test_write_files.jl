#---------------------------------#
#          General Setup          #
#---------------------------------#
# savepath="visualize/"
savepath=""
include(savepath*"write_openvsp.jl")

#---------------------------------#
#            Rotor .bem           #
#---------------------------------#
# - Rotor Data - #
# Header
B_r=5
center_r=[0.12; 0.0; 0.0]
# Blade Data
rs_r = [0.0254, 0.03668832268354582, 0.047977664536709166, 0.059267006389872506, 0.07055634824303583, 0.08184059424811163, 0.09313299361012747, 0.10442539297214334, 0.11571779233415924, 0.127]
D_r=2.0*maximum(rs_r)
rs_r ./= 0.5*D_r
chord_r = [0.089142, 0.079785, 0.0713, 0.063979, 0.057777, 0.052541, 0.048103, 0.044316, 0.041061, 0.038243]
chord_r ./= 0.5*D_r
twist_r = [1.2044866233863267, 1.0322226262144865, 0.9045168848460613, 0.8075987514828161, 0.732200527796661, 0.6721088416504963, 0.6230650896694556, 0.5821371187101886, 0.5471432672077023, 0.5165476454202417]*180.0/pi
num_sections_r = length(rs_r)

# - Write Rotor File - #
write_openvsp_bem(savepath=savepath, filename="rotor",
rotor_name="rotor",
num_sections=num_sections_r,
num_blades=B_r,
diameter=D_r,
center=center_r,
# Blade Data
radius=rs_r,
chord=chord_r,
twist=twist_r,
)

#---------------------------------#
#           Stator .bem           #
#---------------------------------#
# - Stator Data - #
# Header
B_s=9
center_s=[0.22; 0.0; 0.0]
# Blade Data
D_s=D_r
rs_s=rs_r
num_sections_s = length(rs_s)
chord_s=0.05*ones(num_sections_s)#/(0.5*D_s)
twist_s=90.0*ones(num_sections_s)

# - Write Stator File - #
write_openvsp_bem(
    # File Info
    savepath=savepath,
    filename="stator",
    # Header
    rotor_name="stator",
    num_sections=num_sections_s,
    num_blades=B_s,
    diameter=D_s,
    center=center_s,
    # Blade Data
    radius=rs_s,
    chord=chord_s,
    twist=twist_s,
   )

#---------------------------------#
#            Duct .dat            #
#---------------------------------#
# - Duct Coordinates - #
# TODO: need to shift duct coordinates up to actual position to see if that works for placing things in OpenVSP
duct_coordinates = [0.304466 0.158439; 0.294972 0.158441; 0.28113 0.158423; 0.266505 0.158365; 0.251898 0.158254; 0.237332 0.158088; 0.222751 0.157864; 0.208123 0.157586; 0.193399 0.157258; 0.178507 0.156897; 0.16349 0.156523; 0.148679 0.156177; 0.134222 0.155902; 0.12 0.155721; 0.106044 0.155585; 0.092531 0.155498; 0.079836 0.155546; 0.067995 0.155792; 0.057025 0.156294; 0.046983 0.157103; 0.037937 0.158256; 0.029956 0.159771; 0.02311 0.161648; 0.017419 0.163862; 0.012842 0.166404; 0.009324 0.169289; 0.006854 0.172546; 0.005484 0.176154; 0.005242 0.180005; 0.006112 0.184067; 0.00809 0.188086; 0.011135 0.192004; 0.015227 0.19579; 0.020339 0.199393; 0.026403 0.202735; 0.033312 0.205736; 0.040949 0.208332; 0.049193 0.210487; 0.057935 0.212174; 0.067113 0.21339; 0.076647 0.214136; 0.086499 0.214421; 0.09661 0.214255; 0.10695 0.213649; 0.117508 0.212618; 0.12838 0.211153; 0.139859 0.209267; 0.151644 0.207051; 0.163586 0.204547; 0.175647 0.201771; 0.187807 0.198746; 0.20002 0.19549; 0.212269 0.192017; 0.224549 0.188335; 0.236794 0.18447; 0.249026 0.180416; 0.261206 0.176188; 0.273301 0.171796; 0.28524 0.16727; 0.29644 0.162842; 0.304542 0.159526]

_, dle = findmin(duct_coordinates[:,1])
dir, _ = findmin(duct_coordinates[:,2])
duct_coordinates[:,2] .+= maximum(rs_r) - dir

ductlower = duct_coordinates[1:dle,:]
ductupper = duct_coordinates[dle:end,:]

# - Write Duct File - #
write_openvsp_dat(savepath=savepath, filename="duct.dat", airfoilname="Duct", upper_coordinates=ductupper, lower_coordinates=ductlower)

#---------------------------------#
#             Hub .dat            #
#---------------------------------#
# - Hub Coordinates - #
hubupper = [0.0 0.0; 0.000586 0.005293; 0.002179 0.010047; 0.004736 0.014551; 0.008231 0.018825; 0.012632 0.022848; 0.01788 0.026585; 0.023901 0.030001; 0.030604 0.033068; 0.0379 0.035771; 0.045705 0.038107; 0.053933 0.040075; 0.06254 0.04169; 0.071451 0.042966; 0.08063 0.043916; 0.090039 0.044561; 0.09968 0.044922; 0.109361 0.044999; 0.12 0.044952; 0.135773 0.04495; 0.151899 0.04493; 0.16806 0.044913; 0.184232 0.044898; 0.200407 0.044882; 0.21658 0.044866; 0.232723 0.044847; 0.248578 0.044839; 0.262095 0.044564; 0.274184 0.043576; 0.285768 0.041795; 0.296701 0.039168; 0.306379 0.035928]

hublower = hubupper
hublower[:,1] .*-1.0

# - Write Hub File - #
write_openvsp_dat(savepath=savepath, filename="hub.dat", airfoilname="Duct", upper_coordinates=hubupper, lower_coordinates=hublower)

