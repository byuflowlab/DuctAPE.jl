# echo 'converting wake to body aic'
# convert -delay 10 -loop 0 w2baic*.png w2baic.gif
# echo 'converting rotor to body aic'
# convert -delay 10 -loop 0 r2baic*.png r2baic.gif
# echo 'converting body to rotor vhat-z'
# convert -delay 20 -loop 0 b2rvhatz*.png b2rvhatz.gif
# echo 'converting body to rotor vhat-r'
# convert -delay 20 -loop 0 b2rvhatr*.png b2rvhatr.gif
# echo 'converting body to body vhat-z'
# convert -delay 10 -loop 0 b2bvhatz*.png b2bvhatz.gif
# echo 'converting body to body vhat-r'
# convert -delay 10 -loop 0 b2bvhatr*.png b2bvhatr.gif
# echo 'converting rotor to body vhat-z'
# convert -delay 10 -loop 0 r2bvhatz*.png r2bvhatz.gif
# echo 'converting rotor to body vhat-r'
# convert -delay 10 -loop 0 r2bvhatr*.png r2bvhatr.gif
# echo 'converting wake to body vhat-z'
# convert -delay 10 -loop 0 w2bvhatz*.png w2bvhatz.gif
# echo 'converting wake to body vhat-r'
# convert -delay 10 -loop 0 w2bvhatr*.png w2bvhatr.gif
# echo 'converting body to wake vhat-z'
# convert -delay 3 -loop 0 b2wvhatz*.png b2wvhatz.gif
# echo 'converting body to wake vhat-r'
# convert -delay 3 -loop 0 b2wvhatr*.png b2wvhatr.gif
# echo 'converting rotor to wake vhat-z'
# convert -delay 3 -loop 0 r2wvhatz*.png r2wvhatz.gif
# echo 'converting rotor to wake vhat-r'
# convert -delay 3 -loop 0 r2wvhatr*.png r2wvhatr.gif
echo 'converting wake to wake vhat-z'
convert -delay 3 -loop 0 w2wvhatz*.png w2wvhatz.gif
echo 'converting wake to wake vhat-r'
convert -delay 3 -loop 0 w2wvhatr*.png w2wvhatr.gif
echo 'complete.'
