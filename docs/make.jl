using Documenter
using PanelLocalProjections

DocMeta.setdocmeta!(PanelLocalProjections, :DocTestSetup, :(using PanelLocalProjections);
    recursive = true)

makedocs(;
    modules = [PanelLocalProjections],
    sitename = "PanelLocalProjections.jl",
    authors = "Andrea Recine",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages = ["Home" => "index.md",
        "Guide" => ["guide/example.md", "guide/data.md", "guide/formula.md",
            "guide/inference.md"],
        "Methods" => ["methods/absorption.md", "methods/spj.md"],
        "Replications" => ["replication/mss2026.md"],
        "Audit" => ["audit/mss_r.md"],
        "API" => "api.md"],
    checkdocs = :exports,
    # docs/paper, docs/references and docs/notes are local, gitignored reference
    # material; only the pages listed above are built.
    pagesonly = true
)

deploydocs(; repo = "github.com/andrerecio/PanelLocalProjections.jl", devbranch = "main",
    push_preview = false)
