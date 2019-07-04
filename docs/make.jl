using Documenter, CatmullRom

makedocs(
    modules = [CatmullRom],
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    pages = Any[
        "Overview" => "index.md",
        "Examples" => "Examples.md",
        "Circle" => "Circle.md",
        "Piriform" => "Piriform.md",
        "Spherical Spiral" => "SphericalSpiral.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/CatmullRom.jl.git",
    target = "build"
)
