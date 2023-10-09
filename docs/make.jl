using Aerosol
using Documenter

DocMeta.setdocmeta!(Aerosol, :DocTestSetup, :(using Aerosol); recursive=true)

makedocs(;
    modules=[Aerosol],
    authors="EarthSciML authors and contributors",
    sitename="Aerosol.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aerosol.earthsci.dev",
        assets=String[],
        #size_threshold=10000000,
    ),
    repo = Remotes.GitHub("EarthSciML", "Aerosol.jl"),
    pages=[
        "Home" => "index.md",
        "Thermodynamics" => [
            "Isorropia" => [
                "Overview" => "isorropia/overview.md",
                "Examples" => "isorropia/examples.md",
                "Implementation details" => "isorropia/implementation.md",
            ],
        ],
        "API" => "api.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(;
    repo="github.com/EarthSciML/Aerosol.jl",
)
