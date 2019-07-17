using Documenter

#    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),

makedocs(
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    pages = Any[
        "Overview" => "index.md",
        "Perspecive" => "perspective.md",
        "Three Functions" => "threefunctions.md",
        "Examples" => "Examples.md",
        "Points along a path" => "pointsalongapath.md",
        "A few hints" => "hints.md",
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/CatmullRom.jl.git",
    target = "build"
)
