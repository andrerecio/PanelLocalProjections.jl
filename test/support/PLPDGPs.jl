"""
    PLPDGPs

Data-generating processes and plain-CSV I/O for the committed synthetic fixtures in
`test/fixtures/synthetic/`. The module is included both by the test suite
(`@testmodule PLPSupport`) and by `dev/make_fixtures.jl`. It deliberately has no
dependency on the package under test.

Every fixture is a deterministic function of a `StableRNG` seed. The arithmetic uses only
`+ - * /`, `sqrt`, `floor`, `randn` and `rand`, with every reduction an explicit
left-to-right loop, so regenerated values agree bit for bit across platforms, Julia
versions and compiler flags (the tests still compare numerically). Changing a DGP
invalidates every committed Stata/R output computed on that fixture; add a new fixture
instead.
"""
module PLPDGPs

using Dates: Date, Month, firstdayofquarter
using StableRNGs: StableRNG

export FIXTURES, make_fixture, write_csv, read_csv

const BURNIN = 200

# Columns are kept as an ordered list of name => vector pairs; `missing` marks an absent
# value (written as an empty CSV field).
const Columns = Vector{Pair{String, Vector}}

# Left-to-right mean. `sum` may reassociate (SIMD/pairwise), and its result then depends
# on compiler flags such as `--check-bounds=yes` under `Pkg.test`.
function _seqmean(v::AbstractVector{Float64})
    acc = 0.0
    for x in v
        acc += x
    end
    return acc / length(v)
end

"""
Panel DGP shared by `fx_balanced` and `fx_unbalanced` (simulated on the full N × T grid):

    X_t   = 0.4 X_{t-1} + eX_t                     aggregate shock (same for all units)
    Z_t   ~ N(0, 1)                                 omitted aggregate shock
    λ_t   = 0.3 eλ_t                                time effect
    x_it  = 0.2 μ_i + 0.6 x_{i,t-1} + ex_it         unit-specific shock, correlated with μ_i
    w_it  = 0.5 w_{i,t-1} + ew_it                   control
    sv_it = s_i + 0.3 ζ_it                          time-varying exposure
    y_it  = μ_i + λ_t + 0.5 y_{i,t-1} - 0.6 x_it + 0.8 s_i X_t + 0.3 w_it + γ_i Z_t + e_it

with μ_i ~ N(0,1), γ_i = 0.5 + 0.5 N(0,1), s_i = 1 + N(0,1). Unit groupings:
`grp = mod(id-1, 5) + 1` (5 groups) and `reg = mod(id-1, 6) + 1` (6 regions). Each unit
belongs to exactly one group and one region.
"""
function _panel_grid(seed::Integer, N::Integer, T::Integer)
    rng = StableRNG(seed)
    TT = BURNIN + T
    X = zeros(TT)
    Z = zeros(TT)
    λ = zeros(TT)
    for t in 1:TT
        eX = randn(rng)
        X[t] = (t == 1 ? 0.0 : 0.4 * X[t - 1]) + eX
        Z[t] = randn(rng)
        λ[t] = 0.3 * randn(rng)
    end
    n = N * T
    id = Vector{Int}(undef, n)
    tt = Vector{Int}(undef, n)
    y = Vector{Float64}(undef, n)
    x = similar(y)
    Xc = similar(y)
    s = similar(y)
    sv = similar(y)
    w = similar(y)
    grp = Vector{Int}(undef, n)
    reg = Vector{Int}(undef, n)
    r = 0
    for i in 1:N
        μ = randn(rng)
        γ = 0.5 + 0.5 * randn(rng)
        si = 1.0 + randn(rng)
        xp = 0.0
        wp = 0.0
        yp = μ / 0.5     # stationary mean of y given μ_i (λ, X, Z have mean zero)
        for t in 1:TT
            xt = 0.2 * μ + 0.6 * xp + randn(rng)
            wt = 0.5 * wp + randn(rng)
            ζ = randn(rng)
            yt = μ + λ[t] + 0.5 * yp - 0.6 * xt + 0.8 * si * X[t] + 0.3 * wt + γ * Z[t] +
                 randn(rng)
            if t > BURNIN
                r += 1
                id[r] = i
                tt[r] = t - BURNIN
                y[r] = yt
                x[r] = xt
                Xc[r] = X[t]
                s[r] = si
                sv[r] = si + 0.3 * ζ
                w[r] = wt
                grp[r] = mod(i - 1, 5) + 1
                reg[r] = mod(i - 1, 6) + 1
            end
            xp, wp, yp = xt, wt, yt
        end
    end
    return Pair{String, Vector}[
        "id" => id, "t" => tt, "y" => y, "x" => x, "X" => Xc, "s" => s,
        "sv" => sv, "w" => w, "grp" => grp, "reg" => reg]
end

function _select(cols::Columns, keep::AbstractVector{Bool})
    Columns([k => v[keep] for (k, v) in cols])
end

function _col(cols::Columns, name::AbstractString)
    for (k, v) in cols
        k == name && return v
    end
    throw(KeyError(name))
end

"`fx_balanced`: N = 30, T = 60, the full grid of `_panel_grid`."
fx_balanced() = _panel_grid(20261020, 30, 60)

"""
`fx_unbalanced`: N = 40, t ∈ 1..70 from `_panel_grid`, then

- staggered entry/exit: unit i is observed on `entry_i..exit_i`, with entry_i ∈ 1..15 and
  exit_i ∈ 56..70;
- period 35 is absent for every unit;
- about 4% of the remaining rows are removed at random (internal holes);
- unit 38 keeps only t = 30, 31; units 39 and 40 keep only t = 20 and t = 50 (singletons);
- `w` is missing on about 1% of the remaining rows.
"""
function fx_unbalanced()
    N, T = 40, 70
    cols = _panel_grid(20261021, N, T)
    rng = StableRNG(7)
    entry = [1 + floor(Int, 15 * rand(rng)) for _ in 1:N]
    exit_ = [56 + floor(Int, 15 * rand(rng)) for _ in 1:N]
    id = _col(cols, "id")
    t = _col(cols, "t")
    n = length(id)
    keep = trues(n)
    for r in 1:n
        i, ti = id[r], t[r]
        keep[r] = entry[i] <= ti <= exit_[i] && ti != 35
        u = rand(rng)   # one draw per row, whatever the outcome, keeps the stream aligned
        keep[r] && u < 0.04 && (keep[r] = false)
        i == 38 && (keep[r] = ti in (30, 31))
        i == 39 && (keep[r] = ti == 20)
        i == 40 && (keep[r] = ti == 50)
    end
    out = _select(cols, keep)
    w = Vector{Union{Missing, Float64}}(_col(out, "w"))
    for r in eachindex(w)
        rand(rng) < 0.01 && (w[r] = missing)
    end
    return Columns([k == "w" ? (k => w) : (k => v) for (k, v) in out])
end

"""
`fx_dated`: N = 6 units, quarterly `Date` column (first day of the quarter), 2001Q1–2010Q4;
unit 2 has no observation in 2005Q3. `y_it = 0.5 y_{i,t-1} + 0.4 x_it + e_it`, x iid N(0,1).
"""
function fx_dated()
    rng = StableRNG(20261022)
    N, T = 6, 40
    ids = Int[]
    dates = Date[]
    ys = Float64[]
    xs = Float64[]
    for i in 1:N
        yp = 0.0
        for k in 1:(BURNIN + T)
            x = randn(rng)
            yk = 0.5 * yp + 0.4 * x + randn(rng)
            if k > BURNIN
                q = k - BURNIN - 1
                d = firstdayofquarter(Date(2001, 1, 1) + Month(3 * q))
                if !(i == 2 && d == Date(2005, 7, 1))
                    push!(ids, i)
                    push!(dates, d)
                    push!(ys, yk)
                    push!(xs, x)
                end
            end
            yp = yk
        end
    end
    return Pair{String, Vector}["id" => ids, "date" => dates, "y" => ys, "x" => xs]
end

"""
`fx_cdols`: the Ugarte-Ruiz (2026, eqs. 1–2) DGP with correlated effects, N = 100, T = 50:

    s_it = 0.8 s_{i,t-1} + v_it,               v ~ N(0, 1)
    y_it = 0.8 y_{i,t-1} + s_it + u_i + e_it,  e ~ N(0, 1)
    u_i  = α (1 - ρ_y)(1 + s̄_i - s̄),          α = 10

where s̄_i and s̄ are means over the retained window t = 1..T. Both processes run a
200-period burn-in; y starts at u_i / (1 - ρ_y), its mean given u_i. This DGP is our own
documented version, not a reproduction of the paper's `xtarsim` draws.
"""
function fx_cdols()
    rng = StableRNG(20261023)
    N, T = 100, 50
    ρy, ρs, β, α = 0.8, 0.8, 1.0, 10.0
    TT = BURNIN + T
    S = zeros(TT, N)
    for i in 1:N
        sp = 0.0
        for k in 1:TT
            S[k, i] = ρs * sp + randn(rng)
            sp = S[k, i]
        end
    end
    sbar_i = [_seqmean(@view S[(BURNIN + 1):TT, i]) for i in 1:N]
    sbar = _seqmean(sbar_i)
    u = [α * (1 - ρy) * (1 + sbar_i[i] - sbar) for i in 1:N]
    ids = Int[]
    ts = Int[]
    ys = Float64[]
    ss = Float64[]
    for i in 1:N
        yp = u[i] / (1 - ρy)
        for k in 1:TT
            yk = ρy * yp + β * S[k, i] + u[i] + randn(rng)
            if k > BURNIN
                push!(ids, i)
                push!(ts, k - BURNIN)
                push!(ys, yk)
                push!(ss, S[k, i])
            end
            yp = yk
        end
    end
    return Pair{String, Vector}["id" => ids, "t" => ts, "y" => ys, "s" => ss]
end

"""
`fx_tlahr`: the Almuzara–Sancibrián (2026, §2) stylized design plus a unit effect,
N = 200, T = 40:

    y_it = μ_i + 0.7 y_{i,t-1} + β_i X_t + Z_t + u_it

with X_t, Z_t ~ N(0, 1) iid, β_i ~ N(1, 1), exposure s_i = β_i, μ_i ~ N(0, 1), and
u_it ~ N(0, σ_u²), σ_u² = N(1 - R̄²)/R̄² for R̄² = 0.9 (A&S fn 10). `sv = s + 0.5 ζ_it`
is a time-varying exposure.
"""
function fx_tlahr(; seed = 20261024, N = 200, T = 40)
    rng = StableRNG(seed)
    ρ, R2 = 0.7, 0.9
    σu = sqrt(N * (1 - R2) / R2)
    TT = BURNIN + T
    X = [randn(rng) for _ in 1:TT]
    Z = [randn(rng) for _ in 1:TT]
    ids = Int[]
    ts = Int[]
    ys = Float64[]
    Xs = Float64[]
    ss = Float64[]
    svs = Float64[]
    for i in 1:N
        β = 1.0 + randn(rng)
        μ = randn(rng)
        yp = μ / (1 - ρ)
        for k in 1:TT
            ζ = randn(rng)
            yk = μ + ρ * yp + β * X[k] + Z[k] + σu * randn(rng)
            if k > BURNIN
                push!(ids, i)
                push!(ts, k - BURNIN)
                push!(ys, yk)
                push!(Xs, X[k])
                push!(ss, β)
                push!(svs, β + 0.5 * ζ)
            end
            yp = yk
        end
    end
    return Pair{String, Vector}[
        "id" => ids, "t" => ts, "y" => ys, "X" => Xs, "s" => ss, "sv" => svs]
end

"`fx_tlahr_unbal`: `fx_tlahr` with about 10% of the rows removed at random."
function fx_tlahr_unbal()
    cols = fx_tlahr()
    rng = StableRNG(11)
    keep = [rand(rng) >= 0.10 for _ in 1:length(_col(cols, "id"))]
    return _select(cols, keep)
end

"""
`fx_lpiv`: a panel LP-IV design after Almuzara–Sancibrián (2026, §3.4, eq. 14), N = 60,
T = 50:

    X_t   ~ N(0, 1)                                structural shock (not observed)
    Z_t   ~ N(0, 1)                                other aggregate shock
    ffr_t = 0.5 ffr_{t-1} + X_t + 0.6 Z_t          endogenous aggregate (b₀ = 1)
    mp_t  = X_t + 0.5 ν_t,  ν_t ~ N(0, 1)          proxy (a₀ = 1, measurement error ν)
    w_it  = 0.5 w_{i,t-1} + ew_it                  exogenous control
    y_it  = μ_i + 0.6 y_{i,t-1} + β_i X_t + γ_i Z_t + 0.3 w_it + u_it

with s_i = β_i = 1 + N(0, 1), γ_i = 0.5 + 0.5 N(0, 1), μ_i ~ N(0, 1), u_it ~ N(0, 1);
`sv = s + 0.3 ζ_it` is a time-varying exposure; `grp = mod(id-1, 3) + 1` with dummies `g1`,
`g2`, `g3` (several endogenous regressors: g_k & ffr instrumented by g_k & mp). `X` is the
structural shock, kept for reference.
"""
function fx_lpiv(; seed = 20261025, N = 60, T = 50)
    rng = StableRNG(seed)
    TT = BURNIN + T
    X = zeros(TT)
    Z = zeros(TT)
    ffr = zeros(TT)
    mp = zeros(TT)
    for t in 1:TT
        X[t] = randn(rng)
        Z[t] = randn(rng)
        ffr[t] = (t == 1 ? 0.0 : 0.5 * ffr[t - 1]) + X[t] + 0.6 * Z[t]
        mp[t] = X[t] + 0.5 * randn(rng)
    end
    cols = Dict(k => Float64[] for k in ("y", "X", "ffr", "mp", "s", "sv", "w"))
    ids = Int[]
    ts = Int[]
    grp = Int[]
    for i in 1:N
        si = 1.0 + randn(rng)
        γ = 0.5 + 0.5 * randn(rng)
        μ = randn(rng)
        g = mod(i - 1, 3) + 1
        yp = μ / 0.4
        wp = 0.0
        for k in 1:TT
            wk = 0.5 * wp + randn(rng)
            ζ = randn(rng)
            yk = μ + 0.6 * yp + si * X[k] + γ * Z[k] + 0.3 * wk + randn(rng)
            if k > BURNIN
                push!(ids, i)
                push!(ts, k - BURNIN)
                push!(grp, g)
                push!(cols["y"], yk)
                push!(cols["X"], X[k])
                push!(cols["ffr"], ffr[k])
                push!(cols["mp"], mp[k])
                push!(cols["s"], si)
                push!(cols["sv"], si + 0.3 * ζ)
                push!(cols["w"], wk)
            end
            yp, wp = yk, wk
        end
    end
    dummy(k) = [Int(g == k) for g in grp]
    return Pair{String, Vector}[
        "id" => ids, "t" => ts, "y" => cols["y"], "X" => cols["X"], "ffr" => cols["ffr"],
        "mp" => cols["mp"], "s" => cols["s"], "sv" => cols["sv"], "w" => cols["w"],
        "grp" => grp, "g1" => dummy(1), "g2" => dummy(2), "g3" => dummy(3)]
end

"""
`fx_lpiv_unbal`: `fx_lpiv` with staggered entry (unit i enters at t ∈ 1..10), period 25
absent for every unit, and about 8% of the remaining rows removed at random.
"""
function fx_lpiv_unbal()
    cols = fx_lpiv()
    rng = StableRNG(12)
    entry = [1 + floor(Int, 10 * rand(rng)) for _ in 1:60]
    id = _col(cols, "id")
    t = _col(cols, "t")
    keep = trues(length(id))
    for r in eachindex(id)
        u = rand(rng)   # one draw per row keeps the stream aligned
        keep[r] = t[r] >= entry[id[r]] && t[r] != 25 && u >= 0.08
    end
    return _select(cols, keep)
end

"""
    mss_toy(rng; N = 30, T = 60, ρ = 0.8, β = -0.6, δ = 0.2, η = 0.2)

One draw of the Mei–Sheng–Shi (2026) main-text Monte Carlo design (their `simul_LP.R`,
lagY = 0, as used for Figures 1–3 and `toy_-0.6.xlsx`):

    x_it = δ + ρ x_{i,t-1} + u_it,   x_i0 ~ N(0, 1) (no burn-in)
    y_it = α_i + β x_it + ε_it,      α_i = η √T x̄_i + χ

with u, ε ~ N(0, 1), x̄_i the mean of x_i1..x_iT, and χ ~ N(0, 1) one draw per replication
(as in their code). The true response of y_{t+h} to x_t is β ρ^h. Not a committed fixture:
the Monte Carlo tests draw from it.
"""
function mss_toy(rng; N = 30, T = 60, ρ = 0.8, β = -0.6, δ = 0.2, η = 0.2)
    x = Matrix{Float64}(undef, T, N)
    for i in 1:N
        xp = randn(rng)
        for t in 1:T
            xp = δ + ρ * xp + randn(rng)
            x[t, i] = xp
        end
    end
    χ = randn(rng)
    id = Int[]
    tt = Int[]
    ys = Float64[]
    xs = Float64[]
    for i in 1:N
        α = η * sqrt(T) * _seqmean(view(x, :, i)) + χ
        for t in 1:T
            push!(id, i)
            push!(tt, t)
            push!(xs, x[t, i])
            push!(ys, α + β * x[t, i] + randn(rng))
        end
    end
    return (id = id, t = tt, y = ys, x = xs)
end

"Name → generator for every committed synthetic fixture."
const FIXTURES = (fx_balanced = fx_balanced, fx_unbalanced = fx_unbalanced,
    fx_dated = fx_dated, fx_cdols = fx_cdols, fx_tlahr = fx_tlahr,
    fx_tlahr_unbal = fx_tlahr_unbal, fx_lpiv = fx_lpiv, fx_lpiv_unbal = fx_lpiv_unbal)

make_fixture(name::Symbol) = getfield(FIXTURES, name)()

_fmt(::Missing) = ""
_fmt(v::Integer) = string(v)
_fmt(v::Date) = string(v)
_fmt(v::AbstractFloat) = isnan(v) ? "" : repr(Float64(v))

"""
    write_csv(path, cols)

Write `cols` as a comma-separated file with a header row. Floats use the shortest
round-trip representation (`repr`), and `missing`/`NaN` become empty fields.
"""
function write_csv(path::AbstractString, cols::Columns)
    n = length(last(first(cols)))
    open(path, "w") do io
        println(io, join(first.(cols), ","))
        for r in 1:n
            println(io, join((_fmt(v[r]) for (_, v) in cols), ","))
        end
    end
    return path
end

function _parsecol(fields::Vector{SubString{String}})
    if all(f -> !isempty(f) && tryparse(Int, f) !== nothing, fields)
        return parse.(Int, fields)
    elseif all(f -> !isempty(f) && occursin(r"^\d{4}-\d{2}-\d{2}$", f), fields)
        return Date.(fields)
    end
    vals = [isempty(f) ? missing : parse(Float64, f) for f in fields]
    return any(ismissing, vals) ? Vector{Union{Missing, Float64}}(vals) :
           Vector{Float64}(vals)
end

"""
    read_csv(path) -> NamedTuple

Read a file written by `write_csv` into a column table (a NamedTuple of vectors).
Integer columns become `Vector{Int}`, ISO dates `Vector{Date}`, and the rest `Float64`;
any empty field turns a column into `Vector{Union{Missing,Float64}}`.
"""
function read_csv(path::AbstractString)
    lines = readlines(path)
    names = split(lines[1], ",")
    rows = [split(l, ","; keepempty = true) for l in lines[2:end] if !isempty(l)]
    cols = [_parsecol([row[j] for row in rows]) for j in eachindex(names)]
    return NamedTuple{Tuple(Symbol.(names))}(Tuple(cols))
end

end # module
