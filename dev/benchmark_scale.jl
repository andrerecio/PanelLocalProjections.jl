#!/usr/bin/env julia
"""
Compustat-scale benchmark of `panellp` (Phase 1 exit criterion: ≤ 5 min, ≤ 2 GB RSS).

    julia --project=. -t 1 -e 'push!(LOAD_PATH, "dev"); include("dev/benchmark_scale.jl")'

Synthetic firm panel: N = 15,500 firms over T = 124 quarters, 8 sectors, staggered entry
and exit (each firm active for a contiguous window of at least 20 quarters), firm and
sector × quarter fixed effects, and an exposure × aggregate shock. Horizons 0–20; one-way
and two-way clustered and Driscoll–Kraay variances.
"""

using PanelLocalProjections, StatsModels, StableRNGs, Printf

function firm_panel(; N = 15_500, T = 124, S = 8, seed = 4242)
    rng = StableRNG(seed)
    mp = [randn(rng) for _ in 1:T]                 # aggregate shock
    sectorshock = [randn(rng) for _ in 1:S, _ in 1:T]
    firm = Int[]
    quarter = Int[]
    sector = Int[]
    y = Float64[]
    x = Float64[]
    lev = Float64[]
    mps = Float64[]
    for i in 1:N
        len = 20 + floor(Int, (T - 20) * rand(rng))
        start = 1 + floor(Int, (T - len + 1) * rand(rng))
        s = 1 + floor(Int, S * rand(rng))
        μ = randn(rng)
        β = 1.0 + 0.5 * randn(rng)
        yp = μ / 0.4
        xp = 0.0
        levp = 0.5 * randn(rng)
        for t in start:(start + len - 1)
            xt = 0.5 * xp + randn(rng)
            levt = 0.9 * levp + 0.3 * randn(rng)
            yt = μ + 0.6 * yp + 0.3 * xt + β * levt * mp[t] + sectorshock[s, t] + randn(rng)
            push!(firm, i)
            push!(quarter, t)
            push!(sector, s)
            push!(y, yt)
            push!(x, xt)
            push!(lev, levt)
            push!(mps, mp[t])
            yp, xp, levp = yt, xt, levt
        end
    end
    return (; firm, quarter, sector, y, x, lev, mp = mps)
end

println("generating the panel…")
data = firm_panel()
println("rows: ", length(data.firm))
f = @formula(leads(y) ~
             lev & mp + lev + x + lags(y, 1) + fe(firm) + fe(sector) & fe(quarter))

GC.gc()
t_est = @elapsed r = panellp(f, data; id = :firm, time = :quarter, horizon = 20,
    shock = "lev & mp")
t_vcov = @elapsed begin
    v1 = vcov(Clustered(:unit), r)
    v2 = vcov(Clustered(:unit, :time), r)
    v3 = vcov(DriscollKraaySE(), r)
end
s = summarize(r, v2)
@printf("estimation: %.1f s   variances (3 specs): %.2f s   max RSS: %.0f MB\n",
    t_est, t_vcov, Sys.maxrss() / 2^20)
@printf("obs by horizon: %s\n", join(nobspath(r), ", "))
@printf("sweeps by horizon: %s\n", join((fit.iterations for fit in r.fits), ", "))
@printf("IRF (h = 0, 4, 8, 20): %s\n",
    join((@sprintf("%.4f", s.coef[h + 1]) for h in (0, 4, 8, 20)), ", "))
@printf("stored result: %.1f MB\n", Base.summarysize(r) / 2^20)
