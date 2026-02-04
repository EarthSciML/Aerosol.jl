using Aerosol
using Documenter

DocMeta.setdocmeta!(Aerosol, :DocTestSetup, :(using Aerosol); recursive = true)

makedocs(;
    modules = [Aerosol],
    authors = "EarthSciML authors and contributors",
    sitename = "Aerosol.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://aerosol.earthsci.dev",
        assets = String[]        #size_threshold=10000000,
    ),
    repo = Remotes.GitHub("EarthSciML", "Aerosol.jl"),
    pages = [
        "Home" => "index.md",
        "Population Dynamics" => "dynamics.md",
        "Single Particle Dynamics" => "single_particle_dynamics.md",
        "Size Distribution" => "size_distribution.md",
        "VBS" => "VBS.md",
        "Organic Aerosol" => "organic_aerosol.md",
        "Mass Transfer" => "mass_transfer.md",
        "Nucleation" => "nucleation.md",
        "Aqueous Chemistry" => "aqueous_chemistry.md",
        "Thermodynamics" => [
            "ISORROPIA" => "isorropia.md",
            "Seinfeld & Pandis Ch. 10" => "seinfeld_pandis_ch10.md"
        ],
        "API" => "api.md"
    ],
    warnonly = [:missing_docs]
)

deploydocs(; repo = "github.com/EarthSciML/Aerosol.jl")
