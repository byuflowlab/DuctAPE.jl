julia --project=../../ gen_logo.jl
magick logo.svg -trim logo.svg
