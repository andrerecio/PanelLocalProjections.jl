# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "variances from score sums equal brute-force sandwiches" tags = [:core, :vcov] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    r = panellp(@formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t)), u;
        id = :id, time = :t, horizon = 3, clusters = (:reg,))
    p = PLP.PanelIndex(u.id, u.t)
    h = 2
    fit = r.fits[h + 1]
    rows, M, _ = LSDV.within(r, fit, [u.id, u.t])
    X = M(hcat(Float64.(u.x), PLP._lag(p, Float64.(u.x), 1), PLP._lag(p, Float64.(u.y), 1),
        PLP._nanfloat(u.w))[rows, :])
    e = M(PLP._lead(p, Float64.(u.y), h)[rows]) - X * fit.coef
    B = inv(X' * X)
    U = X .* e
    sums(codes) = (g = unique(codes); reduce(vcat, [sum(U[codes .== c, :]; dims = 1) for c in g]))
    meat(codes) = (S = sums(codes); S' * S)
    id, t, reg = u.id[rows], u.t[rows], u.reg[rows]
    oracle = Dict(
        Robust(; ssc = :none) => B * (U' * U) * B,
        Clustered(:unit; ssc = :none) => B * meat(id) * B,
        Clustered(:time; ssc = :none) => B * meat(t) * B,
        Clustered(:reg; ssc = :none) => B * meat(reg) * B,
        Clustered(:unit, :time; ssc = :none, psd = false) => B * (meat(id) + meat(t) - U' * U) * B,
        Clustered(:reg, :time; ssc = :none, psd = false) =>
            B * (meat(reg) + meat(t) - meat(string.(reg, "_", t))) * B)
    for (spec, V) in oracle
        @test vcov(spec, r).V[h + 1] ≈ V rtol = 1e-11
    end
    # Driscoll–Kraay with calendar lags: absent periods (here t = 35) count as zero sums
    m = 3
    G = zeros(70, size(U, 2))
    for (k, tt) in enumerate(t)
        G[tt, :] .+= U[k, :]
    end
    S = G' * G
    for j in 1:m
        Γ = G[(j + 1):end, :]' * G[1:(end - j), :]
        S .+= (1 - j / (m + 1)) .* (Γ + Γ')
    end
    @test vcov(DriscollKraaySE(m; ssc = :none), r).V[h + 1] ≈ B * S * B rtol = 1e-11
    pos = vcov(DriscollKraaySE(m; ssc = :none, timelags = :position), r).V[h + 1]
    @test !(pos ≈ B * S * B)
    # small-sample factors
    n, k = fit.dof.nobs, fit.dof.k
    @test vcov(Robust(), r).V[h + 1] ≈ n / (n - k - fit.dof.df_a) .* (B * (U' * U) * B) rtol = 1e-11
    G_id = length(unique(id))
    T_h = length(unique(t))
    # unit effects nested in the unit clusters: (n−1)/(n − K − T) · G/(G−1) with one component
    @test vcov(Clustered(:unit), r).V[h + 1] ≈
          (n - 1) / (n - k - T_h) * G_id / (G_id - 1) .* (B * meat(id) * B) rtol = 1e-11
    @test vcov(Clustered(:unit), r).df[h + 1][1] == G_id - 1
end

@testitem "Driscoll–Kraay bandwidth rules" tags = [:core, :vcov] begin
    const PLP = PanelLocalProjections
    nw, as, mss = DriscollKraaySE(:nw94), DriscollKraaySE(:as2026), DriscollKraaySE(:mss2026)
    for T in 25:100
        x = 4 * (T / 100)^(2 / 9)
        # A&S's kernel bandwidth M = ⌈x⌉ means M − 1 lags: xtscc's ⌊x⌋ unless x is an integer
        @test PLP._dk_lags(as, T) == (isinteger(x) ? floor(Int, x) - 1 : PLP._dk_lags(nw, T))
        @test PLP._dk_lags(mss, T) == floor(Int, T^(1 / 4))
    end
    @test PLP._dk_lags(DriscollKraaySE(50), 10) == 9      # capped at T − 1
    @test_throws ArgumentError DriscollKraaySE(:foo)
    @test_throws ArgumentError DriscollKraaySE(-1)
end

@testitem "PSD repair of two-way variances" tags = [:core, :vcov] begin
    using LinearAlgebra
    const PLP = PanelLocalProjections
    V = [1.0 0.9 0.3; 0.9 0.5 0.0; 0.3 0.0 -0.05]
    sd = [2.0, 0.5, 1e-3]
    for s in (nothing, sd)
        W, fixed = PLP._psd_clip(V, s)
        @test fixed
        @test minimum(eigvals(Symmetric(s === nothing ? W : W .* (s * s')))) > -1e-12
    end
    # the reghdfe/xtlp metric clips S V S, S = diag(sd)
    W, _ = PLP._psd_clip(V, sd)
    E = eigen(Symmetric(V .* (sd * sd')))
    @test W ≈ (E.vectors * Diagonal(max.(E.values, 0)) * E.vectors') ./ (sd * sd')
    @test PLP._psd_clip([2.0 0.1; 0.1 1.0], nothing) == ([2.0 0.1; 0.1 1.0], false)
    @test Clustered(:unit, :time).psd === :stata
    @test Clustered(:unit, :time; psd = false).psd === :none
    @test Clustered(:unit).psd === :none
    @test_throws ArgumentError Clustered(:unit, :unit)
    @test_throws ArgumentError Clustered(:a, :b, :c)
    # with an intercept (reghdfe reporting _cons): clip the full matrix, then replace the
    # slope block by the clip of the slope block alone
    Vi = [4.0 0.3 0.2; 0.3 1.0 0.9; 0.2 0.9 0.5]
    Vi[2, 3] = Vi[3, 2] = 1.2                  # slope block indefinite
    sdi = [1.0, 0.5, 2.0]
    W, fixed = PLP._psd_clip(Vi, sdi, 1)
    S = Vi .* (sdi * sdi')
    full = PLP._clip(S)[1]
    full[2:3, 2:3] = PLP._clip(S[2:3, 2:3])[1]
    @test fixed && W ≈ full ./ (sdi * sdi')
end

@testitem "PSD repair with an intercept follows reghdfe" tags = [:core, :vcov] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w)      # pooled: an intercept
    r = panellp(f, b; id = :id, time = :t, horizon = 3, clusters = (:reg,))
    rs = panellp(f, merge(b, (; x = b.x .+ 1000, w = b.w .- 1000)); id = :id, time = :t,
        horizon = 3, clusters = (:reg,))
    cv, cvs = vcov(Clustered(:reg, :time), r), vcov(Clustered(:reg, :time), rs)
    @test [i.psdfixed for i in cv.info] == [true, true, true, false]
    # the slopes' variance does not depend on the regressors' levels (the intercept's does)
    for k in 1:4
        @test diag(cvs.V[k])[2:5] ≈ diag(cv.V[k])[2:5] rtol = 1e-6
    end
end

@testitem "undefined and negative variances" tags = [:core, :vcov] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels, StableRNGs
    rng = StableRNG(5)
    N, T = 30, 8
    d = (id = repeat(1:N, inner = T), t = repeat(1:T, N), y = randn(rng, N * T),
        x = randn(rng, N * T))
    # h = 7 reads y at t + 7: only t = 1 is left, a single period
    r = panellp(@formula(leads(y) ~ x), d; id = :id, time = :t, horizon = 7)
    for spec in (Clustered(:unit, :time), Clustered(:time), DriscollKraaySE())
        cv = @test_logs (:warn, r"undefined at h = 7") vcov(spec, r)
        @test all(isnan, cv.V[8]) && all(isnan, cv.df[8]) && all(isfinite, cv.V[7])
    end
    # a two-way variance without PSD repair can be negative: NaN standard error
    b = fixture(:fx_balanced)
    rb = panellp(@formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id)), b; id = :id,
        time = :t, horizon = 3, clusters = (:reg,))
    cv = @test_logs (:warn, r"negative variances at h = 1") vcov(
        Clustered(:time, :reg; ssc = :none, psd = false), rb)
    @test cv.V[2][1, 1] < 0
    s = summarize(rb, cv)
    @test isnan(s.se[2]) && isnan(s.lower[2]) && all(isfinite, s.se[[1, 3, 4]])
    @test isnan(stderror(cv; term = "x")[2])
end

@testitem "cluster specifications" tags = [:core, :vcov] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    u = fixture(:fx_unbalanced)
    r = panellp(@formula(leads(y) ~ x + w + fe(id) + fe(grp) & fe(t)), u; id = :id, time = :t,
        horizon = 1, clusters = (:reg,))
    # the Mei–Sheng–Shi factor exists for unit and unit × time clustering only
    @test vcov(Clustered(:unit; ssc = :mss), r) isa PanelLPVcov
    @test vcov(Clustered(:id, :time; ssc = :mss, psd = :raw), r) isa PanelLPVcov
    @test_throws ArgumentError vcov(Clustered(:time; ssc = :mss), r)
    @test_throws ArgumentError vcov(Clustered(:reg; ssc = :mss), r)
    # two names of one dimension
    @test_throws ArgumentError vcov(Clustered(:unit, :id), r)
    @test_throws ArgumentError vcov(Clustered(:t, :time), r)
    # an intersection matches the fixed-effect interaction whatever the order of its columns
    @test vcov(Clustered((:t, :grp)), r).V == vcov(Clustered((:grp, :t)), r).V
    err = try
        vcov(Clustered((:reg, :t)), r)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("clusters = ((:reg, :t),)", err.msg)
end

@testitem "intersection and nested clusters equal brute-force sandwiches" tags = [:core, :vcov] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    r = panellp(@formula(leads(y) ~ x + w + lags(x, 1) + fe(id) + fe(grp) & fe(t)), u;
        id = :id, time = :t, horizon = 1, clusters = ((:reg, :t), :reg))
    fit = r.fits[2]
    gt = LSDV.inter(u.grp, u.t)
    rows, M, rk = LSDV.within(r, fit, [u.id, gt])
    p = PLP.PanelIndex(u.id, u.t)
    X = hcat(Float64.(u.x), PLP._nanfloat(u.w), PLP._lag(p, Float64.(u.x), 1))[rows, :]
    s = LSDV.sandwich(M(X), M(PLP._lead(p, Float64.(u.y), 1)[rows]))
    U, B = s.U, s.B
    id, t, reg, g = u.id[rows], u.t[rows], u.reg[rows], gt[rows]
    rt = LSDV.inter(reg, t)
    V(spec) = vcov(spec, r).V[2]
    meat(c) = LSDV.meat(U, c)
    @test fit.dof.df_a == rk
    @test V(Clustered((:grp, :t); ssc = :none)) ≈ B * meat(g) * B rtol = 1e-10
    @test V(Clustered((:reg, :t); ssc = :none)) ≈ B * meat(rt) * B rtol = 1e-10
    @test V(Clustered(:unit, (:reg, :t); ssc = :none, psd = false)) ≈
          B * (meat(id) + meat(rt) - U' * U) * B rtol = 1e-10
    # units are nested in regions: the two-way meat is the region meat (pair kind :first)
    @test V(Clustered(:unit, :reg; ssc = :none, psd = false)) ≈ B * meat(reg) * B rtol = 1e-10
    # two intersection dimensions, and an intersection with a column
    @test V(Clustered((:grp, :t), (:reg, :t); ssc = :none, psd = false)) ≈
          B * (meat(g) + meat(rt) - meat(LSDV.inter(g, rt))) * B rtol = 1e-10
    @test V(Clustered((:grp, :t), :reg; ssc = :none, psd = false)) ≈
          B * (meat(g) + meat(reg) - meat(LSDV.inter(g, reg))) * B rtol = 1e-10
end
