using DuctTAPE
using Documenter

DocMeta.setdocmeta!(DuctTAPE, :DocTestSetup, :(using DuctTAPE); recursive=true)

makedocs(;
    modules=[DuctTAPE],
    authors="Judd Mehr",
    repo="https://github.com/byuflowlab/DuctTAPE.jl/blob/{commit}{path}#{line}",
    sitename="DuctTAPE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://flow.byu.edu/DuctTAPE.jl",
        assets=String[],
    ),
    pages=[
        "Intro" => "index.md",
        "Quick Start" => "tutorial.md",
        "Guided Examples" => "examples.md",
        "API Reference" => "reference.md",
        "Theory" => "theory.md",
    ],
)

deploydocs(; repo="github.com/byuflowlab/DuctTAPE.jl", devbranch="main")
