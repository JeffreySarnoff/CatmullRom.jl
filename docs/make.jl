using Documenter, CatmullRom

makedocs(
    modules = [CatmullRom],
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    pages = Any[
        "Overview" => "index.md",
        "Guide" => "userguide.md",
        "Notes" => "notes.md",
        "Refs" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/CatmullRom.jl.git",
    target = "build"
)
