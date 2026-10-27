# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "SPJ equals 2θ̂ − ½(θ̂_a + θ̂_b) on explicit halves" tags = [:core, :spj] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id))
    for est_fes in (f, @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t)))
        # h = 1 on the balanced panel: rows t = 2..59 for every unit (58, even)
        rs = panellp(est_fes, b; id = :id, time = :t, horizon = 1, estimator = SplitPanelJackknife())
        rf = panellp(est_fes, b; id = :id, time = :t, horizon = 1)
        ra = panellp(est_fes, b; id = :id, time = :t, horizon = 1, subset = 2 .<= b.t .<= 30)
        rb = panellp(est_fes, b; id = :id, time = :t, horizon = 1, subset = 31 .<= b.t .<= 59)
        @test coef(rs, 1) ≈ 2 .* coef(rf, 1) .- (coef(ra, 1) .+ coef(rb, 1)) ./ 2 rtol = 1e-10
        # with an even number of periods per unit the two split rules coincide
        rfl = panellp(est_fes, b; id = :id, time = :t, horizon = 1,
            estimator = SplitPanelJackknife(; split = :floor))
        @test coef(rfl, 1) ≈ coef(rs, 1) rtol = 1e-12
    end
end

@testitem "SPJ scores equal the paper's d*" tags = [:core, :spj] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra
    const PLP = PanelLocalProjections
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + w + fe(id))
    r = panellp(f, b; id = :id, time = :t, horizon = 1, estimator = SplitPanelJackknife())
    fit = r.fits[2]
    rows = findall(fit.esample)              # t = 1..59, 59 rows per unit: odd
    p = PLP.PanelIndex(b.id, b.t)
    X = hcat(Float64.(b.x), Float64.(b.w))[rows, :]
    y = PLP._lead(p, Float64.(b.y), 1)[rows]
    id, t = b.id[rows], b.t[rows]
    demean(A, g) = A .- reduce(vcat, [sum(A[g .== g[k], :]; dims = 1) ./ count(==(g[k]), g) for k in eachindex(g)])
    Xt, yt = demean(X, id), demean(reshape(y, :, 1), id)
    # half a = the first ⌈59/2⌉ = 30 periods. For odd T_i the code's d = 2X̃ − X̃_half is
    # W − 2W̄ + W̄_own (xtlp notes); build the variance from that d by hand
    ina = t .<= 30
    meanover(A, sel) = reduce(vcat, [sum(A[(id .== id[k]) .& sel, :]; dims = 1) ./ count((id .== id[k]) .& sel) for k in eachindex(id)])
    dcode = 2 .* Xt .- ((X .- meanover(X, ina)) .* ina .+ (X .- meanover(X, .!ina)) .* .!ina)
    e = vec(yt) - Xt * coef(r, 1)
    B = inv(Xt' * Xt)
    S = reduce(vcat, [sum((dcode .* e)[id .== g, :]; dims = 1) for g in unique(id)])
    @test vcov(Clustered(:unit; ssc = :none), r).V[2] ≈ B * (S' * S) * B rtol = 1e-10
    # on an even split (h = 2: rows t = 1..58) the package's scores are the paper's d*
    # (Assumption 2; appendix B.2.3 with time effects): W − W̄_{i,other half}, minus its
    # cross-sectional mean at t with time effects. Robust and time-clustered variances
    # tell d from d* row by row (unit clustering cannot: they differ by a unit constant).
    for twoway in (false, true)
        f2 = twoway ? @formula(leads(y) ~ x + w + fe(id) + fe(t)) : f
        r2 = panellp(f2, b; id = :id, time = :t, horizon = 2, estimator = SplitPanelJackknife(),
            absorb = AlternatingProjections(; tol = 1e-14))
        rows2 = findall(r2.fits[3].esample)
        W = hcat(Float64.(b.x), Float64.(b.w))[rows2, :]
        y2 = PLP._lead(p, Float64.(b.y), 2)[rows2]
        id2, t2 = b.id[rows2], b.t[rows2]
        ina2 = t2 .<= 29
        um(A) = reduce(vcat, [sum(A[id2 .== id2[k], :]; dims = 1) ./ count(==(id2[k]), id2) for k in eachindex(id2)])
        tm(A) = reduce(vcat, [sum(A[t2 .== t2[k], :]; dims = 1) ./ count(==(t2[k]), t2) for k in eachindex(t2)])
        wt(A) = twoway ? A .- um(A) .- tm(A) .+ sum(A; dims = 1) ./ size(A, 1) : A .- um(A)
        other = reduce(vcat, [sum(W[(id2 .== id2[k]) .& (ina2 .!= ina2[k]), :]; dims = 1) ./
                              count((id2 .== id2[k]) .& (ina2 .!= ina2[k])) for k in eachindex(id2)])
        dstar = twoway ? (W .- other) .- tm(W .- other) : W .- other
        Wt = wt(W)
        e2 = vec(wt(reshape(y2, :, 1))) - Wt * coef(r2, 2)
        B2 = inv(Wt' * Wt)
        U = dstar .* e2
        St = reduce(vcat, [sum(U[t2 .== tt, :]; dims = 1) for tt in unique(t2)])
        @test vcov(Robust(; ssc = :none), r2).V[3] ≈ B2 * (U' * U) * B2 rtol = 1e-10
        @test vcov(Clustered(:time; ssc = :none), r2).V[3] ≈ B2 * (St' * St) * B2 rtol = 1e-10
    end
end

@testitem "SPJ split rules, singletons and degenerate halves" tags = [:core, :spj] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    const PLP = PanelLocalProjections
    # per-unit split by count: ⌈n/2⌉ or ⌊n/2⌋ first rows (in time order) in half 1
    p = PLP.PanelIndex([1, 1, 1, 2, 2, 2, 2, 1], [3, 1, 2, 1, 2, 3, 4, 4])
    rows = [1, 2, 3, 4, 5, 6, 7]                     # unit 1: t = 3, 1, 2; unit 2: t = 1..4
    @test PLP._spj_halves(p, rows, :ceil) == Int8[2, 1, 1, 1, 1, 2, 2]
    @test PLP._spj_halves(p, rows, :floor) == Int8[2, 1, 2, 1, 1, 2, 2]
    u = fixture(:fx_unbalanced)
    f = @formula(leads(y) ~ x + lags(x, 1) + w + fe(id))
    rc = panellp(f, u; id = :id, time = :t, horizon = 2, estimator = SplitPanelJackknife())
    rfl = panellp(f, u; id = :id, time = :t, horizon = 2, estimator = SplitPanelJackknife(; split = :floor))
    @test !(coef(rc, 0) ≈ coef(rfl, 0))
    # unit fixed effects only: singleton rows demean to zero, so keeping them changes nothing
    rk = panellp(f, u; id = :id, time = :t, horizon = 2, singletons = :keep, estimator = SplitPanelJackknife())
    @test all(coef(rk, h) ≈ coef(rc, h) for h in 0:2)
    @test nobs(rk, 0) > nobs(rc, 0)
    # SPJ needs a unit fixed effect
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + fe(t)), u; id = :id, time = :t,
        horizon = 1, estimator = SplitPanelJackknife())
    # halves with fewer rows than regressors at every horizon: nothing can be estimated
    tiny = (id = [1, 1, 1, 2, 2, 2], t = [1, 2, 3, 1, 2, 3], y = randn(6), x = randn(6), w = randn(6))
    err = try
        panellp(@formula(leads(y) ~ x + w + fe(id)), tiny; id = :id, time = :t, horizon = 0,
            estimator = SplitPanelJackknife())
    catch e
        e
    end
    @test err isa ArgumentError && occursin("too few periods", err.msg)
    @test_throws ArgumentError SplitPanelJackknife(; split = :middle)
end

@testitem "a regressor absorbed within one half makes the horizon singular" tags = [:core, :spj] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    b = fixture(:fx_balanced)
    # w2 is frozen at unit constants (non-dyadic: demeaning leaves round-off, not zeros)
    # from t = 31 on: identified in the full sample, absorbed by the unit effects inside
    # half b while that half starts at t = 31 (h = 0, 1: 60 and 59 rows per unit), not at
    # h = 2 (58 rows: half b starts at t = 30)
    frozen(v) = [t >= 31 ? 0.1 * i + 1 / 3 : x for (i, t, x) in zip(b.id, b.t, v)]
    bw = merge(b, (; w2 = frozen(b.w), x2 = frozen(b.x)))
    f = @formula(leads(y) ~ x + w2 + fe(id))
    r = @test_logs (:warn, r"rank-deficient at h = 0, 1 \(w2 absorbed") panellp(f, bw;
        id = :id, time = :t, horizon = 2, estimator = SplitPanelJackknife())
    @test [f.status for f in r.fits] == [:spj_half_singular, :spj_half_singular, :ok]
    @test all(isnan, coef(r, 0)) && !any(isnan, coef(r, 2))
    @test occursin("w2 absorbed", sprint(show, MIME"text/plain"(), r))
    # the same for the shock; with no estimable horizon, an error
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x2 + w + fe(id)), bw; id = :id,
        time = :t, horizon = 1, estimator = SplitPanelJackknife())
    # the fixed-effects fit itself is fine
    @test !any(isnan, coef(panellp(f, bw; id = :id, time = :t, horizon = 1), 1))
end

@testitem "SPJ equals the hand-split formula on an unbalanced panel" tags = [:core, :spj] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    p = PLP.PanelIndex(u.id, u.t)
    # rows of `es` in the first half of each unit (its first ⌈n_i/2⌉ or ⌊n_i/2⌋ in time order)
    function firsthalf(es, split)
        ina = falses(length(es))
        for i in unique(u.id)
            rows = findall(es .& (u.id .== i))
            rows = rows[sortperm(u.t[rows])]
            ina[rows[1:(split === :ceil ? cld(length(rows), 2) : fld(length(rows), 2))]] .= true
        end
        return ina
    end
    # projection off the span of dummies for the given codes on `rows` (any connectivity)
    function annihilate(rows, codes...)
        D = reduce(hcat, [Float64.(c[rows] .== permutedims(unique(c[rows]))) for c in codes])
        Q = Matrix(qr(D, ColumnNorm()).Q)[:, 1:rank(D)]
        return A -> A - Q * (Q' * A)
    end
    # without lags, units with two or three rows enter, with single-row halves
    for split in (:ceil, :floor), twoway in (false, true), rhs in ("x + lags(x, 1) + w", "x + w")
        f = eval(Meta.parse("@formula(leads(y) ~ $rhs + fe(id)" * (twoway ? " + fe(t))" : ")")))
        rs = panellp(f, u; id = :id, time = :t, horizon = 2, estimator = SplitPanelJackknife(; split))
        rf = panellp(f, u; id = :id, time = :t, horizon = 2)
        for h in 0:2
            es = esample(rs, h)
            ina = firsthalf(es, split)
            ra = panellp(f, u; id = :id, time = :t, horizon = h:h, subset = es .& ina, singletons = :keep)
            rb = panellp(f, u; id = :id, time = :t, horizon = h:h, subset = es .& .!ina, singletons = :keep)
            @test es == esample(rf, h)
            @test coef(rs, h) ≈ 2 .* coef(rf, h) .- (coef(ra, h) .+ coef(rb, h)) ./ 2 rtol = 1e-10
        end
        # the variance: full-sample bread, SPJ residuals, scores d = 2X̃ − X̃_half with every
        # fixed effect (time effects too) re-estimated in each half, by explicit dummies
        h = 1
        rows = findall(esample(rs, h))
        cols = rhs == "x + w" ? [Float64.(u.x), PLP._nanfloat(u.w)] :
               [Float64.(u.x), PLP._lag(p, Float64.(u.x), 1), PLP._nanfloat(u.w)]
        W, y = reduce(hcat, cols)[rows, :], PLP._lead(p, Float64.(u.y), h)[rows]
        codes = twoway ? (u.id, u.t) : (u.id,)
        M = annihilate(rows, codes...)
        Wt, yt = M(W), M(y)
        Wh = similar(Wt)
        ina = firsthalf(esample(rs, h), split)[rows]
        for half in (ina, .!ina)
            sel = findall(half)
            Wh[sel, :] = annihilate(rows[sel], codes...)(W[sel, :])
        end
        e = yt - Wt * coef(rs, h)
        B = inv(Wt' * Wt)
        U = (2 .* Wt .- Wh) .* e
        S = reduce(vcat, [sum(U[u.id[rows] .== g, :]; dims = 1) for g in unique(u.id[rows])])
        @test vcov(Clustered(:unit; ssc = :none), rs).V[h + 1] ≈ B * (S' * S) * B rtol = 1e-9
        @test vcov(Robust(; ssc = :none), rs).V[h + 1] ≈ B * (U' * U) * B rtol = 1e-9
    end
end
