using Documenter
using PanelLocalProjections

DocMeta.setdocmeta!(PanelLocalProjections, :DocTestSetup, :(using PanelLocalProjections);
    recursive = true)

# The cross-check page is crosscheck/README.md, the map of reference commands to
# `panellp` calls, copied at build time (one source; the copy is gitignored).
let readme = read(joinpath(@__DIR__, "..", "crosscheck", "README.md"), String)
    page = replace(readme,
        r"^# crosscheck/\n"m =>
            "# Cross-checks\n\nThis page is " *
            "`crosscheck/README.md` of the repository.\n")
    write(joinpath(@__DIR__, "src", "crosscheck.md"), page)
end

makedocs(;
    modules = [PanelLocalProjections],
    sitename = "PanelLocalProjections.jl",
    authors = "Andrea Recine",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    pages = ["Home" => "index.md",
        "Guide" => ["guide/example.md", "guide/data.md", "guide/formula.md",
            "guide/inference.md"],
        "Methods" => ["methods/absorption.md", "methods/spj.md", "methods/cdols.md",
            "methods/tlahr.md"],
        "Replications" => ["replication/mss2026.md", "replication/ugarte2026.md",
            "replication/firm_level_wip.md"],
        "Cross-checks" => ["crosscheck.md", "deviations.md"],
        "Audit" => ["audit/xtlp.md", "audit/mss_r.md", "audit/ugarte.md",
            "audit/as_code.md", "audit/tidymacro.md", "audit/ramey.md"],
        "API" => "api.md"],
    checkdocs = :exports,
    # docs/paper, docs/references and docs/notes are local, gitignored reference
    # material; only the pages listed above are built.
    pagesonly = true
)

deploydocs(; repo = "github.com/andrerecio/PanelLocalProjections.jl", devbranch = "main",
    push_preview = false)
