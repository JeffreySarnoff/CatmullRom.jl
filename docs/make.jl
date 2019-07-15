using Documenter, CatmullRom

makedocs(
    modules = [CatmullRom],
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = Any[
        "Overview" => "index.md",
        "Seeing The Contiguity" => "perspective.md",
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
