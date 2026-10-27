# Build the US state panel of the panel LP-IV example (docs/src/examples/us_states.md) from
# public sources, into docs/data/us_states/ (gitignored: the data are not committed).
#
#     julia --project examples/us_states/build.jl
#
# Sources (downloaded once into docs/data/us_states/raw/; needs `unzip` on the PATH):
# - BEA, State quarterly personal income (SQINC, bulk download apps.bea.gov/regional/zip/):
#   SQINC1 (personal income, line 1; quarterly population and per-capita income are not
#   published for the sample's early years) and SQINC5N (earnings by place of work, line 35,
#   and manufacturing earnings, line 500). US government work, public domain.
# - FRED, GS1: 1-year Treasury constant-maturity rate, monthly (Board of Governors, H.15).
# - Jarociński & Karadi (2020, AEJ: Macroeconomics), "Deconstructing Monetary Policy
#   Surprises — The Role of Information Shocks", updated monetary policy surprises
#   (github.com/marekjarocinski/jkshocks_update_fed, `shocks_fed_jk_m.csv`, column MP_pm),
#   licensed CC BY 4.0.
#
# Output: docs/data/us_states/processed/state_panel.csv, one row per state (50 + DC) and
# quarter, 1999Q1–2019Q4:
#   fips, state, year, quarter, t (1999Q1 = 1, consecutive), lpi (100 × log personal income,
#   nominal, seasonally adjusted), manuf (manufacturing share of earnings, 1998 average), s
#   (manuf standardised across states), gs1 (quarterly average, percent), mp (MP_pm summed
#   over the quarter's months, percentage points).

using Downloads, Statistics
include(joinpath(@__DIR__, "csv.jl"))

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs", "data", "us_states"))
const RAW = joinpath(ROOT, "raw")
const OUT = joinpath(ROOT, "processed")
const FIRST, LAST = 1999, 2019            # sample years
const BASE = 1998                         # exposure base year (before the sample)

const SOURCES = (
    sqinc = ("https://apps.bea.gov/regional/zip/SQINC.zip", "SQINC.zip"),
    gs1 = ("https://fred.stlouisfed.org/graph/fredgraph.csv?id=GS1", "GS1.csv"),
    jk = ("https://raw.githubusercontent.com/marekjarocinski/jkshocks_update_fed/main/shocks_fed_jk_m.csv",
        "shocks_fed_jk_m.csv"))

function fetch(url, file)
    path = joinpath(RAW, file)
    isfile(path) || Downloads.download(url, path)
    return path
end

value(s) = (v = tryparse(Float64, strip(s)); v === nothing ? NaN : v)   # (D), (T), (NA) → NaN

"""
    bea_lines(path, linecodes) -> Dict((fips, linecode) => Dict((year, quarter) => value))

The series of the given line codes for the 50 states and DC (FIPS 01000–56000, multiples
of 1000) in a BEA regional CSV; suppressed or missing cells are NaN.
"""
function bea_lines(path, linecodes)
    lines = readlines(path)
    hdr = csvfields(lines[1])
    qcols = [(k, parse(Int, h[1:4]), parse(Int, h[7:7]))
             for (k, h) in enumerate(hdr) if occursin(r"^\d{4}:Q\d$", h)]
    out = Dict{Tuple{Int, Int}, Dict{Tuple{Int, Int}, Float64}}()
    names = Dict{Int, String}()
    for l in lines[2:end]
        f = csvfields(l)
        length(f) == length(hdr) || continue            # footnotes
        fips = tryparse(Int, strip(f[1]))
        fips === nothing && continue
        (1000 <= fips <= 56000 && fips % 1000 == 0) || continue
        lc = parse(Int, strip(f[5]))
        lc in linecodes || continue
        names[fips] = replace(strip(f[2]), r"\s*\*$" => "")
        out[(fips, lc)] = Dict((y, q) => value(f[k]) for (k, y, q) in qcols)
    end
    return out, names
end

function main()
    mkpath(RAW)
    mkpath(OUT)
    zip = fetch(SOURCES.sqinc...)
    gs1 = fetch(SOURCES.gs1...)
    jk = fetch(SOURCES.jk...)
    sq = joinpath(RAW, "sqinc")
    mkpath(sq)
    want = ["SQINC1__ALL_AREAS_1948_*.csv", "SQINC5N__ALL_AREAS_1998_*.csv"]
    run(`unzip -o -q $zip $want -d $sq`)
    file(pat) = joinpath(sq, only(filter(f -> occursin(pat, f), readdir(sq))))
    pi1, names = bea_lines(file(r"^SQINC1__ALL_AREAS"), (1,))
    e5n, _ = bea_lines(file(r"^SQINC5N__ALL_AREAS"), (35, 500))
    states = sort!(collect(keys(names)))
    @assert length(states) == 51

    # exposure: manufacturing share of earnings by place of work, average of the base year's
    # quarters with both values
    manuf = Dict{Int, Float64}()
    for f in states
        shares = Float64[]
        for q in 1:4
            m, e = e5n[(f, 500)][(BASE, q)], e5n[(f, 35)][(BASE, q)]
            isfinite(m) && isfinite(e) && push!(shares, m / e)
        end
        manuf[f] = mean(shares)
    end
    μ, σ = mean(values(manuf)), std(collect(values(manuf)))

    # quarterly aggregates: GS1 averaged over the quarter's months, MP_pm summed
    g = Dict{Tuple{Int, Int}, Vector{Float64}}()
    for l in readlines(gs1)[2:end]
        f = csvfields(l)
        y, m = parse(Int, f[1][1:4]), parse(Int, f[1][6:7])
        push!(get!(g, (y, cld(m, 3)), Float64[]), value(f[2]))
    end
    lines = readlines(jk)
    hdr = csvfields(lines[1])
    jm = findfirst(==("MP_pm"), hdr)
    mpq = Dict{Tuple{Int, Int}, Float64}()
    for l in lines[2:end]
        f = csvfields(l)
        y, m = parse(Int, f[1]), parse(Int, f[2])
        key = (y, cld(m, 3))
        mpq[key] = get(mpq, key, 0.0) + value(f[jm])
    end

    path = joinpath(OUT, "state_panel.csv")
    open(path, "w") do io
        println(io, "fips,state,year,quarter,t,lpi,manuf,s,gs1,mp")
        for f in states, y in FIRST:LAST, q in 1:4
            gq = g[(y, q)]
            length(gq) == 3 || error("GS1: $(length(gq)) months in $y Q$q")
            t = 4 * (y - FIRST) + q
            lpi = 100 * log(pi1[(f, 1)][(y, q)])
            isfinite(lpi) || error("no personal income for $(names[f]) in $y Q$q")
            println(io, join((f, "\"$(names[f])\"", y, q, t, repr(lpi), repr(manuf[f]),
                    repr((manuf[f] - μ) / σ), repr(mean(gq)), repr(mpq[(y, q)])), ","))
        end
    end
    println("wrote ", relpath(path, normpath(joinpath(@__DIR__, "..", ".."))), ": ",
        length(states), " states × ", 4 * (LAST - FIRST + 1), " quarters")
end

main()
