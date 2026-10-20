#!/usr/bin/env julia
"""
Format the package with the SciML style, from the local `dev` environment.

    julia --project=dev dev/format.jl            # reformat in place
    julia --project=dev dev/format.jl --check    # exit 1 if any file would change

Only the repository's own Julia sources are touched. `docs/paper`, `docs/references` and
`docs/notes` are read-only, gitignored reference material and are never formatted.
Adapted from LocalProjections.jl (andrerecio fork), `dev/format.jl`.
"""

using JuliaFormatter

const ROOT = dirname(@__DIR__)
const TARGETS = ["src", "ext", "test", "dev", "docs/make.jl"]

function main(args)
    check = "--check" in args
    allclean = true
    for t in TARGETS
        path = joinpath(ROOT, t)
        ispath(path) || continue
        clean = format(path; overwrite = !check)
        allclean &= clean
        println(clean ? "  ok          " : (check ? "  needs format " : "  reformatted  "), t)
    end
    println(allclean ? "\nAll files formatted." :
            (check ? "\nRun `julia --project=dev dev/format.jl` to fix." :
             "\nSome files were reformatted; review `git diff`."))
    return allclean
end

if abspath(PROGRAM_FILE) == @__FILE__
    clean = main(ARGS)
    "--check" in ARGS && !clean && exit(1)
end
