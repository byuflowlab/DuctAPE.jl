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
            "Examples" => "DuctAPE/examples.md",
            "API Reference" => "DuctAPE/public_api.md",
            "API Reference" => "DuctAPE/private_api.md",
            "Index" => "DuctAPE/api_index.md",
            "Theory" => "DuctAPE/theory.md",
        ],
        "C\$^4\$Blade" => [
            "Intro" => "C4Blade/intro.md",
            "Examples" => "C4Blade/examples.md",
            "API Reference" => "C4Blade/api.md",
        ],
    ],
    sitename="DuctAPE.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    checkdocs=:exports,
)

deploydocs(; repo="https://github.com/byuflowlab/DuctAPE.jl.git", devbranch="main")
