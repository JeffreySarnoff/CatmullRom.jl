using Documenter, CatmullRom

makedocs(
    modules = [CatmullRom],
    sitename = "CatmullRom.jl",
    authors = "Jeffrey Sarnoff",
    pages = Any[
        "Overview" => "index.md",
        "Guide" => "guide.md",
        "Notes" => "notes.md",
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JeffreySarnoff/CatmullRom.jl.git",
    target = "build"
)
