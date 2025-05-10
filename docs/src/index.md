```@meta
CurrentModule = Aerosol
```

# Aerosol.jl: Symbolic equation-based aerosol modeling

Documentation for [Aerosol](https://github.com/EarthSciML/Aerosol.jl).


## Installation

To install Aerosol.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("Aerosol")
```

## Features

Currently, we have a mostly-working version of the [Isorropia aerosol thermodynamics model](@ref "isorropia").

## Contributing

...coming soon

## Reproducibility

```@raw html
<details><summary>The documentation of this EarthSciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/EarthSciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/EarthSciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
