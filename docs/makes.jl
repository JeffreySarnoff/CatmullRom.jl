using Documenter, CatmullRom

makedocs(
    sitename = "CatmullRom",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Jeffrey Sarnoff",
    pages = [
        "Introduction" => "index.md",
        "Overview" => "overview.md",
        "Quick Start Guide" => "quickstart.md",
        "Examples" => "examples.md",
        "Style Guide" => "style.md",
    ],
    assets = [
        "assets/jump-logo-with-text.svg",
        "assets/numfocus-logo.png"
    ]
)

deploydocs(
    repo   = "github.com/JeffreySarnoff/CatmullRom.jl.git",
)
