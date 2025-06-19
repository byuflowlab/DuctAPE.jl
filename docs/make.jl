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
        size_threshold_ignore=["DuctAPE/tutorial.md"],
    ),
    pages=[
        "Home" => "index.md",
        "DuctAPE" => [
            "Quick Start" => "DuctAPE/tutorial.md",
            "API Reference" => [
                "Analysis" => "DuctAPE/api/analysis.md",
                "Required Inputs" => "DuctAPE/api/inputs.md",
                "Options" => "DuctAPE/api/options.md",
                "Precompiled Caches" => "DuctAPE/api/caches.md",
                "Outputs" => "DuctAPE/api/outputs.md",
                "API Index" => "DuctAPE/api/api_index.md",
            ],
            "Visualization" => "DuctAPE/visualization.md",
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
            # "Polar Modification" => "C4Blade/corrections.md",
            "API Index" => "C4Blade/api.md",
        ],
    ],
    sitename="DuctAPE.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    # checkdocs=:exports,
    checkdocs=:none,
)

deploydocs(; repo="https://github.com/byuflowlab/DuctAPE.jl.git", devbranch="main")
