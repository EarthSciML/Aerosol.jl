using AerosolMTK
using Documenter

DocMeta.setdocmeta!(AerosolMTK, :DocTestSetup, :(using AerosolMTK); recursive=true)

makedocs(;
    modules=[AerosolMTK],
    authors="EarthSciML authors and contributors",
    repo="https://github.com/EarthSciML/AerosolMTK.jl/blob/{commit}{path}#{line}",
    sitename="AerosolMTK.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EarthSciML.github.io/AerosolMTK.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EarthSciML/AerosolMTK.jl",
)
