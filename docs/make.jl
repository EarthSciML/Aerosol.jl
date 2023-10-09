using Aerosol
using Documenter

DocMeta.setdocmeta!(Aerosol, :DocTestSetup, :(using Aerosol); recursive=true)

makedocs(;
    modules=[Aerosol],
    authors="EarthSciML authors and contributors",
    sitename="Aerosol.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EarthSciML.github.io/Aerosol.jl",
        assets=String[],
        repolink="https://github.com/EarthSciML/Aerosol.jl/blob/{commit}{path}#{line}",
        #size_threshold=10000000,
    ),
    pages=[
        "Home" => "index.md",
        "Thermodynamics" => [
            "Isorropia" => "isorropia.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/EarthSciML/Aerosol.jl",
)
