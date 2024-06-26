using Documenter
using DuctAPE

# - LaTeX Stuff - #
mathengine = Documenter.MathJax(
    Dict(:TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        # :Macros => Dict(:ket => ["|#1\\rangle", 1], :bra => ["\\langle#1|", 1]),
    ))
)

# - Make the Docs - #
makedocs(;
    modules=[DuctAPE, DuctAPE.C4Blade],
    format=Documenter.HTML(;
        repolink="https://github.com/byuflowlab/DuctAPE.jl/blob/{commit}{path}#L{line}",
        edit_link="main",
        mathengine=mathengine,
    ),
    pages=[
        "Home" => "index.md",
        "DuctAPE" => [
            "Getting Started" => "DuctAPE/tutorial.md",
            "Advanced Usage" => [
                "Options" => "DuctAPE/advanced_usage/option.md",
                "Preallocation" => "DuctAPE/advanced_usage/precompilation.md",
                "Multi-Point Analysis" => "DuctAPE/advanced_usage/multi_point.md",
                # "Multi-Rotor Analysis" => "DuctAPE/advanced_usage/multi_rotor.md",
                "Outputs" => "DuctAPE/advanced_usage/outputs.md",
                # "Manual Geometry" => "DuctAPE/advanced_usage/manual_repaneling.md",
            ],
            "API" => [
                "Public API Reference" => "DuctAPE/api/public_api.md",
                "Private API Reference" => [
                    "Prelims" => "DuctAPE/api/private_prelims.md",
                    "Preprocess" => "DuctAPE/api/private_preprocess.md",
                    "Process" => "DuctAPE/api/private_process.md",
                    "Postprocess" => "DuctAPE/api/private_postprocess.md",
                    "Utilities" => "DuctAPE/api/private_utilities.md",
                ],
                "API Index" => "DuctAPE/api/api_index.md",
            ],
            "Theory" => "DuctAPE/theory.md",
        ],
        "C\$^4\$Blade" => [
            "Intro" => "C4Blade/intro.md",
            "Blade Element Types" => [
                "C4Blade/airfoil_types/DFDC.md",
                "C4Blade/airfoil_types/CCBlade.md",
                "C4Blade/airfoil_types/actuator_disk.md",
                "C4Blade/airfoil_types/cascade.md",
            ],
            "Polar Modification" => "C4Blade/corrections.md",
            "API Reference" => "C4Blade/api.md",
        ],
    ],
    sitename="DuctAPE.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    checkdocs=:exports,
)

deploydocs(; repo="https://github.com/byuflowlab/DuctAPE.jl.git", devbranch="main")
