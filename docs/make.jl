using Aerosol
using Documenter

DocMeta.setdocmeta!(Aerosol, :DocTestSetup, :(using Aerosol); recursive=true)

makedocs(;
    modules=[Aerosol],
    authors="EarthSciML authors and contributors",
    repo="https://github.com/EarthSciML/Aerosol.jl/blob/{commit}{path}#{line}",
    sitename="Aerosol.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EarthSciML.github.io/Aerosol.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EarthSciML/Aerosol.jl",
)
