using Documenter
using DuctAPE

makedocs(;
    modules=[DuctAPE],
    format=Documenter.HTML(;
        repolink="https://github.com/byuflowlab/DuctAPE.jl/blob/{commit}{path}#L{line}",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "tutorial.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
        "Theory" => "theory.md",
    ],
    sitename="DuctAPE.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    checkdocs=:exports,
)

deploydocs(; repo="https://github.com/byuflowlab/DuctAPE.jl.git", devbranch="main")
