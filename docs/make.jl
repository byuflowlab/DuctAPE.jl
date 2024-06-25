using Documenter
using DuctAPE

makedocs(;
    modules=[DuctAPE, DuctAPE.C4Blade],
    format=Documenter.HTML(;
        repolink="https://github.com/byuflowlab/DuctAPE.jl/blob/{commit}{path}#L{line}",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "DuctAPE" => [
            "Getting Started" => "DuctAPE/tutorial.md",
            "Examples" => [
                "Options" => "DuctAPE/examples/option.md",
                "Outputs" => "DuctAPE/examples/outputs.md",
                "Multi-Point Analysis" => "DuctAPE/examples/multi_point.md",
                "Preallocation" => "DuctAPE/examples/precompilation.md",
                "Manual Geometry" => "DuctAPE/examples/manual_repaneling.md",
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
            "Examples" => [
                           "C4Blade/examples/DFDC.md",
                           "C4Blade/examples/CCBlade.md",
                           "C4Blade/examples/actuator_disk.md",
                           "C4Blade/examples/cascade.md",
                           "C4Blade/examples/corrections.md",
                          ],
            "API Reference" => "C4Blade/api.md",
        ],
    ],
    sitename="DuctAPE.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    checkdocs=:exports,
)

deploydocs(; repo="https://github.com/byuflowlab/DuctAPE.jl.git", devbranch="main")
