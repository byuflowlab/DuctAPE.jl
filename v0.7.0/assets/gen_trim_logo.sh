julia --project=../../ gen_logo.jl
magick mogrify -trim logo.png
