using Documenter, CatmullRom

makedocs(
    modules = [CatmullRom],
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = Any[
        "Overview" => "index.md",
        "Perspecive" => "perspective.md",
        "Three Functions" => "threefunctions.md",
        "Examples" => "Examples.md",
        "A few hints" => "hints.md",
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/CatmullRom.jl.git",
    target = "build"
)
