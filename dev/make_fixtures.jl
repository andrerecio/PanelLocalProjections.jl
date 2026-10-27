#!/usr/bin/env julia
"""
Generate (or check) the committed synthetic fixture panels in `test/fixtures/synthetic/`.

    julia --project=dev dev/make_fixtures.jl            # (re)write every fixture
    julia --project=dev dev/make_fixtures.jl --check    # regenerate in memory and compare

The data-generating processes live in `test/support/PLPDGPs.jl`, which the test suite
shares. Rewriting a fixture invalidates every committed Stata/R output computed on it, so
`--check` is part of each phase gate.
"""

include(joinpath(@__DIR__, "..", "test", "support", "PLPDGPs.jl"))
using .PLPDGPs

const OUTDIR = normpath(joinpath(@__DIR__, "..", "test", "fixtures", "synthetic"))

function _same(a, b)
    (ismissing(a) && ismissing(b)) ||
        (!ismissing(a) && !ismissing(b) &&
         (a == b || (a isa AbstractFloat &&
           b isa AbstractFloat && isapprox(a, b; rtol = 1e-14, atol = 1e-12))))
end

function check()
    ok = true
    for name in keys(FIXTURES)
        path = joinpath(OUTDIR, string(name, ".csv"))
        if !isfile(path)
            println("  missing   ", name)
            ok = false
            continue
        end
        disk = read_csv(path)
        mem = make_fixture(name)
        good = length(mem) == length(disk) &&
               all(k == string(dk) for ((k, _), dk) in zip(mem, keys(disk))) &&
               all(all(_same.(v, disk[Symbol(k)])) for (k, v) in mem)
        println(good ? "  ok        " : "  MISMATCH  ", name)
        ok &= good
    end
    return ok
end

function main(args)
    if "--check" in args
        check() || exit(1)
    else
        mkpath(OUTDIR)
        for name in keys(FIXTURES)
            path = write_csv(joinpath(OUTDIR, string(name, ".csv")), make_fixture(name))
            println("  wrote ", relpath(path, dirname(@__DIR__)))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
