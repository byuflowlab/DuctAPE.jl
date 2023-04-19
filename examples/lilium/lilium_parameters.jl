#=

7-seater Lilium Jet Propulsor Paramters

Estimated from publicly available data (span, length, and weight) and public images

Authors: Judd Mehr,

=#

## -- Major Paramters From Public Images -- ##

#Blade Geometry: estimated from public images on lillium website
Rhub = 0.06
Rtip = 0.26
croot_rotor_measured = 0.0475
ctip_rotor_measured = 0.03
rotor_te_pos = 0.26
rotor_c4_pos = rotor_te_pos - 3.0 * ctip_rotor_measured / 4.0

rotor_root_twist_guess = 60.0
rotor_chord_guess = croot_rotor_measured / sind(rotor_root_twist_guess)
rotor_tip_twist_guess = asind(ctip_rotor_measured / rotor_chord_guess)

#Stator Geometry: estimated from public images on lillium website
#for now assume that stator has zero twist...
stator_le_x = 0.3
stator_root_chord = 0.155
stator_c4_pos = stator_le_x + stator_root_chord / 4.0
stator_tip_chord = 0.07
stator_twist_guess = 0.0

#Duct Geometry: estimated from public images on lillium website
duct_OD = 0.35 #meters
duct_chord = 0.75 #meters
duct_te = [0.75; 0.08]
duct_rle = 0.015 / 0.75
duct_le = [0.0; 0.135]

#Hub Geometry: estimated from public images on lillium website
hub_le = 0.05
hub_te = 0.725
hub_chord = 0.675
cyl_le = rotor_te_pos - croot_rotor_measured
cyl_te = stator_le_x + stator_root_chord / 2.0

# Misc Parameters: From published data and public images
max_weight = 31000.0 #Newtons
nprops = 36 #number of propulsors (according to image on website)
n_rotor_blades = 27 #according to front view image on google
n_stator_blades = 8 #according to front view image on google
