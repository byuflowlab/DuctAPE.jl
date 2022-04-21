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
        canonical="https://byuflowlab.github.io/DuctTAPE.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/DuctTAPE.jl",
    devbranch="main",
)
