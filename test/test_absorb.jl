# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule LSDV begin
    using LinearAlgebra, PanelLocalProjections, StableRNGs
    const PLP = PanelLocalProjections
    # Within-transformed design of one horizon by explicit dummies (the brute-force oracle):
    # projection on an orthonormal basis of their span (`D \ A` is not reliable when the
    # dummies have several redundant columns, e.g. disconnected blocks).
    function within(r, fit, cols::Vector{<:AbstractVector})
        rows = findall(fit.esample)
        D = reduce(hcat, [Float64.(c[rows] .== permutedims(sort(unique(c[rows])))) for c in cols])
        rk = rank(D)
        Q = Matrix(qr(D, ColumnNorm()).Q)[:, 1:rk]
        return rows, A -> A - Q * (Q' * A), rk
    end
    # Cluster meat Σ_g (Σ_{r∈g} u_r)(Σ_{r∈g} u_r)' of scores U by the codes of each row.
    function meat(U, codes)
        S = reduce(vcat, [sum(U[codes .== g, :]; dims = 1) for g in unique(codes)])
        return S' * S
    end
    inter(a, b) = string.(a, "×", b)
    # OLS on within-transformed data: coefficients, bread, scores.
    function sandwich(Xt, yt)
        β = Xt \ yt
        return (; β, B = inv(Xt' * Xt), U = Xt .* (yt - Xt * β))
    end
    # A rotating panel that mixes slowly under alternating projections: unit i is observed
    # in L consecutive periods starting at ⌈i/2⌉ (two units enter each period).
    function chainpanel(; N = 120, L = 3, seed = 7)
        rng = StableRNG(seed)
        id = [i for i in 1:N for _ in 1:L]
        t = [cld(i, 2) + k for i in 1:N for k in 0:(L - 1)]
        n = length(id)
        a, g = randn(rng, N), randn(rng, maximum(t))
        x = randn(rng, n) .+ a[id] .+ g[t]
        w = randn(rng, n) .+ 0.5 .* a[id]
        y = 0.5 .* x .- 0.2 .* w .+ a[id] .+ g[t] .+ randn(rng, n)
        return (; id, t, x, w, y, X = g[t] .+ 0.0)
    end
end

@testitem "two-way absorption equals LSDV" tags = [:core, :absorb] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    for name in (:fx_balanced, :fx_unbalanced), fes in (:one, :two)
        u = fixture(name)
        f = fes === :two ?
            @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t)) :
            @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id))
        r = panellp(f, u; id = :id, time = :t, horizon = 3)
        cv = vcov(Clustered(:unit; ssc = :none), r)
        p = PLP.PanelIndex(u.id, u.t)
        for h in 0:3
            fit = r.fits[h + 1]
            rows, M, rk = LSDV.within(r, fit, fes === :two ? [u.id, u.t] : [u.id])
            X = hcat(Float64.(u.x), PLP._lag(p, Float64.(u.x), 1), PLP._lag(p, Float64.(u.y), 1),
                PLP._nanfloat(u.w))[rows, :]
            y = PLP._lead(p, Float64.(u.y), h)[rows]
            s = LSDV.sandwich(M(X), M(y))
            @test fit.coef ≈ s.β rtol = 1e-10
            # the variance is first order in the demeaning error, the coefficients second
            @test cv.V[h + 1] ≈ s.B * LSDV.meat(s.U, u.id[rows]) * s.B rtol = 1e-9
            @test fit.dof.df_a == rk
            @test fit.converged
        end
    end
end

@testitem "absorbed degrees of freedom on disconnected blocks" tags = [:core, :absorb] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using StatsModels
    b = fixture(:fx_balanced)
    # units 1–15 at t ≤ 30, units 16–30 at t > 30: two connected components
    d = map(v -> v[.!(((b.id .<= 15) .& (b.t .> 30)) .| ((b.id .> 15) .& (b.t .<= 30)))], b)
    r = panellp(@formula(leads(y) ~ x + w + fe(id) + fe(t)), d; id = :id, time = :t, horizon = 2)
    for h in 0:2
        fit = r.fits[h + 1]
        rows, M, rk = LSDV.within(r, fit, [d.id, d.t])
        @test fit.dof.fecomponents[1, 2] == 2 && fit.dof.df_a == rk
        @test fit.coef ≈ M(hcat(d.x, d.w)[rows, :]) \ M(PanelLocalProjections._lead(
            PanelLocalProjections.PanelIndex(d.id, d.t), Float64.(d.y), h)[rows]) rtol = 1e-10
    end
end

@testitem "absorption is invariant to units and levels" tags = [:core, :absorb] setup = [LSDV] begin
    using StatsModels
    # slowly mixing two-way design: thousands of sweeps
    c = LSDV.chainpanel()
    f = @formula(y ~ x + w + fe(id) + fe(t))
    r0 = panellp(f, c; id = :id, time = :t, horizon = 0)
    rows, M, _ = LSDV.within(r0, r0.fits[1], [c.id, c.t])
    s = LSDV.sandwich(M(hcat(c.x, c.w)[rows, :]), M(c.y[rows]))
    se(r) = sqrt.(vcov(Clustered(:unit; ssc = :none), r).V[1][[1, 4]])
    @test coef(r0, 0) ≈ s.β rtol = 1e-12
    @test se(r0) ≈ sqrt.((s.B * LSDV.meat(s.U, c.id[rows]) * s.B)[[1, 4]]) rtol = 1e-8
    for k in (1e-9, 1e-6, 1e6)       # rescaled data: same estimates, same sweeps
        rk = panellp(f, map(v -> v isa Vector{Float64} ? k .* v : v, c); id = :id, time = :t,
            horizon = 0)
        @test coef(rk, 0) ≈ coef(r0, 0) rtol = 1e-11
        @test se(rk) ≈ se(r0) rtol = 1e-10
        @test rk.fits[1].iterations == r0.fits[1].iterations && rk.fits[1].converged
    end
    for shift in (1e4, 1e6)          # levels absorbed by the fixed effects
        rs = panellp(f, merge(c, (; x = c.x .+ shift, w = c.w .- shift)); id = :id, time = :t,
            horizon = 0)
        @test coef(rs, 0) ≈ coef(r0, 0) rtol = 1e-10
        @test se(rs) ≈ se(r0) rtol = 1e-10
    end
    # an aggregate shock is absorbed by the time effects at any scale
    for k in (1.0, 1e-6, 1e-9)
        @test_throws ArgumentError panellp(@formula(y ~ X + x + fe(id) + fe(t)),
            merge(c, (; X = k .* c.X)); id = :id, time = :t, horizon = 0)
    end
end

@testitem "absorption methods" tags = [:core, :absorb] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    f = @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t))
    ra = panellp(f, u; id = :id, time = :t, horizon = 2)
    rp = panellp(f, u; id = :id, time = :t, horizon = 2,
        absorb = AlternatingProjections(; accelerate = false))
    @test all(isapprox(ra.fits[k].coef, rp.fits[k].coef; rtol = 1e-10) for k in 1:3)
    va, vp = vcov(Clustered(:unit, :time), ra), vcov(Clustered(:unit, :time), rp)
    @test all(isapprox(va.V[k], vp.V[k]; rtol = 1e-9) for k in 1:3)
    @test ra.fits[2].iterations < rp.fits[2].iterations
    # OneSweep: unit means, then time means of the unit-demeaned data, once (LP_panel_all.r)
    rows = findall(ra.fits[1].esample)
    v = randn(length(rows))
    function onesweep(v, id, t)
        w = copy(v)
        for g in unique(id); m = id .== g; w[m] .-= sum(w[m]) / count(m); end
        for g in unique(t); m = t .== g; w[m] .-= sum(w[m]) / count(m); end
        w
    end
    fs = PLP._fe_sample([Int32.(u.id), Int32.(u.t)], rows)
    w = copy(v)
    PLP._absorb!(w, fs, OneSweep(), PLP._absorb_buffers(length(w), fs))
    @test w ≈ onesweep(v, u.id[rows], u.t[rows]) rtol = 1e-13
    # one dimension: exact group demeaning
    fs1 = PLP._fe_sample([Int32.(u.id)], rows)
    w1 = copy(v)
    PLP._absorb!(w1, fs1, AlternatingProjections(), PLP._absorb_buffers(length(w1), fs1))
    ids = u.id[rows]
    @test w1 ≈ v .- [sum(v[ids .== ids[k]]) / count(==(ids[k]), ids) for k in eachindex(v)] rtol = 1e-13
end

@testitem "singletons and degrees of freedom" tags = [:core, :absorb] begin
    const PLP = PanelLocalProjections
    # unit 3 is alone; dropping it leaves period 3 with one row, which is dropped next
    a = Int32[1, 1, 2, 2, 3, 1]
    b = Int32[1, 2, 1, 2, 3, 3]
    mask = trues(6)
    @test PLP._drop_singletons!(mask, [a, b]) == 2
    @test mask == BitVector([1, 1, 1, 1, 0, 0])
    # two components: {units 1,2} × {periods 1,2} and {unit 3} × {period 3}
    fs = PLP._fe_sample([Int32[1, 1, 2, 2, 3, 3], Int32[1, 2, 1, 2, 3, 3]], collect(1:6))
    m = PLP._fe_components(fs)
    @test m[1, 2] == 2
    @test PLP._df_fe(length.(fs.counts), m, [1, 2], false) == 3 + 3 - 2
    @test PLP._isnested(Int32[1, 1, 2, 2], 2, Int32[5, 5, 6, 6])
    @test !PLP._isnested(Int32[1, 1, 2, 2], 2, Int32[5, 6, 6, 6])
end

@testitem "absorbed and collinear regressors" tags = [:core, :ols] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    b = fixture(:fx_balanced)
    bw = merge(b, (; w2 = 2 .* b.w, s2 = 3 .* b.s))
    # s is time-invariant: absorbed by the unit effects, omitted (NaN), not a fake regressor
    r = @test_logs (:warn, r"omitted, absorbed.*: s$") panellp(
        @formula(leads(y) ~ x + s + fe(id)), b; id = :id, time = :t, horizon = 1)
    j = findfirst(==("s"), coefnames(r))
    @test isnan(coef(r, 0)[j]) && !r.fits[1].keep[j] && r.fits[1].absorbed[j]
    @test r.fits[1].dof.k == 1
    @test occursin("note: omitted", sprint(show, MIME"text/plain"(), r))
    # every call reports its own omissions (no once-per-session log limit)
    @test_logs (:warn, r"absorbed.*: s2$") panellp(@formula(leads(y) ~ x + s2 + fe(id)), bw;
        id = :id, time = :t, horizon = 1)
    # an aggregate shock is absorbed by the time effects: an error, not a zero coefficient
    err = try
        panellp(@formula(leads(y) ~ X + fe(id) + fe(t)), b; id = :id, time = :t, horizon = 1)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("every horizon", err.msg) &&
          occursin("time effects", err.msg)
    # exact collinearity among regressors
    r2 = @test_logs (:warn, r"omitted, collinear.*: w2$") panellp(
        @formula(leads(y) ~ x + w + w2 + fe(id)), bw; id = :id, time = :t, horizon = 0)
    @test count(r2.fits[1].keep) == 2 && !any(r2.fits[1].absorbed)
    # equal after equilibration, so QR keeps the first-listed column and drops the shock
    @test_throws ArgumentError panellp(@formula(leads(y) ~ w2 + w + fe(id)), bw; id = :id,
        time = :t, horizon = 0, shock = "w")
end

@testitem "a degenerate horizon is singular, not fatal" tags = [:core, :ols] begin
    using StatsModels, StableRNGs
    # staggered treatment dummy, first switched on at t = 6, 7 or 8 (T = 10): with leads up
    # to h, the sample ends at t = 10 − h, where D ≡ 0 for h ≥ 5
    rng = StableRNG(3)
    N, T = 30, 10
    τ = rand(rng, 6:8, N)
    d = (id = repeat(1:N, inner = T), t = repeat(1:T, N), y = randn(rng, N * T),
        x = randn(rng, N * T))
    d = merge(d, (; D = Float64.(d.t .>= τ[d.id])))
    f = @formula(leads(y) ~ D + x + fe(id) + fe(t))
    r = @test_logs (:warn, r"shock is absorbed.*h = 5, 6") panellp(f, d; id = :id, time = :t,
        horizon = 6)
    @test [f.status for f in r.fits] == [fill(:ok, 5); :singular; :singular]
    @test all(isnan, coef(r, 5)) && !any(isnan, coef(r, 4)) && nobs(r, 5) == 150
    @test all(isnan, summarize(r, Clustered(:unit)).se[6:7])
    # one row per unit at the last horizon, singletons kept: the shock demeans to zero there
    d1 = (id = repeat(1:N, inner = 6), t = repeat(1:6, N), y = randn(rng, 6N), x = randn(rng, 6N))
    r1 = @test_logs (:warn, r"h = 5") panellp(@formula(leads(y) ~ x + fe(id)), d1; id = :id,
        time = :t, horizon = 5, singletons = :keep)
    @test r1.fits[6].status === :singular && r1.fits[5].status === :ok
end
