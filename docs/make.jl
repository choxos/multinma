using Documenter
using multinma

makedocs(
    sitename = "multinma.jl",
    modules  = [multinma],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://choxos.github.io/multinma/",
    ),
    pages    = [
        "Home" => "index.md",
        "User Guide" => "guide.md",
        "Datasets" => "datasets.md",
        "API Reference" => "api.md",
    ],
    warnonly = true,
)

deploydocs(
    repo = "github.com/choxos/multinma.git",
    devbranch = "main",
    dirname = "",
    push_preview = true,
)
