# t-LAHR (Almuzara & Sancibrián 2026): lag augmentation, the synthetic information
# criterion, CR0 by period, the HC2/Imbens–Kolesár refinement, the cosine comparator and
# the synthetic series.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule TLAHRSupport begin
    using LinearAlgebra, PanelLocalProjections, StableRNGs, Statistics, StatsModels
    # Brute-force t-LAHR on a within-transformed design Xw (columns = kept coefficients) and
    # response yw with period codes t: per shock j (column) with exposure e_j, the weights are
    # the per-period residuals of e_j on an intercept (`timefe`) and the other exposures;
    # X̄ the per-period weighted means of Xw, ω_t = Σ_{i∈t} w_i²; HC2 with d = 1 − diag(hat
    # of √ω X̄) (fn 24's weighted regression) and the Bell–McCaffrey df from the eigenvalues
    # of G'G (Imbens & Kolesár 2016).
    function oracle(Xw, yw, t, shockcols::Vector{Int}, exposures::Vector{Vector{Float64}},
            timefe::Bool)
        β = Xw \ yw
        e = yw - Xw * β
        B = inv(Xw' * Xw)
        ps = sort(unique(t))
        S = reduce(vcat, [sum(Xw[t .== p, :] .* e[t .== p]; dims = 1) for p in ps])
        out = []
        for (k, j) in enumerate(shockcols)
            w = similar(exposures[k])
            for p in ps
                rr = findall(==(p), t)
                Z = reduce(hcat, [timefe ? [ones(length(rr))] : Vector{Float64}[];
                                  [exposures[l][rr] for l in eachindex(exposures) if l != k]];
                           init = zeros(length(rr), 0))
                w[rr] = exposures[k][rr] - Z * (Z \ exposures[k][rr])
            end
            ω = [sum(abs2, w[t .== p]) for p in ps]
            Xb = reduce(vcat, [sum(Xw[t .== p, :] .* w[t .== p]; dims = 1) ./
                               sum(abs2, w[t .== p]) for p in ps])
            Xs = sqrt.(ω) .* Xb
            # columns the weights annihilate (other exposures of X) are round-off, ~1e-15
            # relative: a loose tolerance drops them
            H = Xs * pinv(Xs; rtol = 1e-10)
            d = 1 .- diag(H)
            a = Xs * pinv(Xs' * Xs; rtol = 1e-20)[:, j]
            G = (I - H) * Diagonal(a ./ sqrt.(d))
            λ = eigvals(Symmetric(G' * G))
            push!(out, (; j, β = β[j], v0 = (B * (S' * S) * B)[j, j],
                v2 = (B * (S' * (S ./ d)) * B)[j, j], ν = sum(λ)^2 / sum(abs2, λ), Xb))
        end
        return out
    end
    # CR2 clustered by period (Bell & McCaffrey 2002) for coefficient k of the pooled
    # regression of y_{t+h} on s·X and s·X_{t−1} without intercept, and its Bell–McCaffrey
    # degrees of freedom under homoskedastic errors: per period the block
    # (I − H_tt)^{−1/2} of the residual maker, G_{·t} = M_{·,t}(I − H_tt)^{−1/2}X_t B e_k
    # and ν = (Σλ)²/Σλ² for the eigenvalues λ of G'G.
    function cr2(d, h, k)
        p = PanelLocalProjections.PanelIndex(d.id, d.t)
        lag(v) = PanelLocalProjections._lag(p, Float64.(v), 1)
        Xall = hcat(d.s .* d.X, d.s .* lag(d.X))
        yall = PanelLocalProjections._lead(p, Float64.(d.y), h)
        rows = findall(r -> all(isfinite, Xall[r, :]) && isfinite(yall[r]), eachindex(yall))
        X, y, tt = Xall[rows, :], yall[rows], d.t[rows]
        B = inv(X' * X)
        β = B * (X' * y)
        u = y - X * β
        H = X * B * X'
        M = I - H
        ps = sort(unique(tt))
        meat = zeros(2, 2)
        G = zeros(length(rows), length(ps))
        for (q, pp) in enumerate(ps)
            rr = findall(==(pp), tt)
            E = eigen(Symmetric(I - H[rr, rr]))
            A = E.vectors * Diagonal(1 ./ sqrt.(E.values)) * E.vectors'
            sc = X[rr, :]' * (A * u[rr])
            meat .+= sc * sc'
            G[:, q] = M[:, rr] * (A * (X[rr, :] * B[:, k]))
        end
        λ = eigvals(Symmetric(G' * G))
        return (; β = β[k], se = sqrt((B * meat * B)[k, k]), ν = sum(λ)^2 / sum(abs2, λ))
    end
    # Almuzara & Sancibrián (2026, §2) stylized design: y_it = ρ y_i,t−1 + β_i X_t + Z_t +
    # u_it, β_i ~ N(1, 1), exposure s_i = β_i (so the response is ρ^h), X, Z ~ N(0, 1),
    # u ~ N(0, σ_u²) with σ_u² = c·N(1 − R̄²)/R̄²: c = 1 is their fn 10 (`:fn10`), c = 1 + ρ²
    # the exact eq. (7) value (`:exact`).
    function stylized(rng; N = 200, T = 30, ρ = 0.7, R2 = 0.9, cal = :fn10, burn = 100)
        σu = sqrt(N * (1 - R2) / R2 * (cal === :fn10 ? 1.0 : 1 + ρ^2))
        X, Z = randn(rng, burn + T), randn(rng, burn + T)
        id, tt, y, XX, s = Int[], Int[], Float64[], Float64[], Float64[]
        for i in 1:N
            β = 1 + randn(rng)
            yp = 0.0
            for k in 1:(burn + T)
                yk = ρ * yp + β * X[k] + Z[k] + σu * randn(rng)
                if k > burn
                    push!(id, i); push!(tt, k - burn); push!(y, yk); push!(XX, X[k])
                    push!(s, β)
                end
                yp = yk
            end
        end
        return (; id, t = tt, y, X = XX, s)
    end
    # Coverage of the 90% intervals of Table 1 (h = 0, 1, 2, 5) in one draw: unit (1W) and
    # two-way (2W, no PSD repair: a negative variance never covers, as their SE = 0) clusters,
    # Driscoll–Kraay with their bandwidth, plain t-LAHR (one lag, CR0, normal) and the
    # recommended one (ICLags, HC2 + Imbens–Kolesár); and the selected p̂.
    const TABLE1_H = (0, 1, 2, 5)
    function table1(d; ρ = 0.7)
        f = @formula(leads(y) ~ s & X + fe(id) + fe(t))
        kw = (; id = :id, time = :t, horizon = 5)
        w = panellp(f, d; kw...)
        la = panellp(f, d; kw..., estimator = LagAugmented())
        pl = panellp(f, d; kw..., estimator = LagAugmented(ConstantLags(1)))
        function cover(r, spec)
            s = summarize(r, spec; level = 0.9)
            return [s.lower[h + 1] <= ρ^h <= s.upper[h + 1] for h in TABLE1_H]
        end
        return (; oneway = cover(w, Clustered(:unit; ssc = :none)),
            twoway = cover(w, Clustered(:unit, :time; ssc = :none, psd = false)),
            dk = cover(w, DriscollKraaySE(:as2026; ssc = :none)),
            plain = cover(pl, TLAHR(; refine = false)),
            recommended = cover(la, TLAHR()), phat = la.lags.p)
    end
    # A stationary panel AR(2) with an aggregate shock times an exposure:
    # y_it = μ_i + 1.2 y_i,t−1 − 0.5 y_i,t−2 + s_i X_t + Z_t + u_it.
    function ar2panel(; N = 200, T = 60, seed = 2)
        rng = StableRNG(seed)
        burn = 100
        X = randn(rng, burn + T)
        Z = randn(rng, burn + T)
        id, tt, y, XX, s = Int[], Int[], Float64[], Float64[], Float64[]
        for i in 1:N
            si, μ = 1 + randn(rng), randn(rng)
            y1 = y2 = μ / 0.3
            for k in 1:(burn + T)
                yk = μ + 1.2 * y1 - 0.5 * y2 + si * X[k] + Z[k] + 3 * randn(rng)
                if k > burn
                    push!(id, i); push!(tt, k - burn); push!(y, yk); push!(XX, X[k])
                    push!(s, si)
                end
                y1, y2 = yk, y1
            end
        end
        return (; id, t = tt, y, X = XX, s)
    end
end

@testitem "LagAugmented: lag rules and the added columns" tags = [:core, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    d = fixture(:fx_tlahr)
    H = 6
    fit(f, rule; kw...) = panellp(f, d; id = :id, time = :t, horizon = H,
        estimator = LagAugmented(rule), kw...)
    f = @formula(leads(y) ~ s & X + fe(id) + fe(t))
    r = fit(f, ConstantLags(2))
    @test coefnames(r) == ["s & X", "L1.y", "L2.y", "L1.(s & X)", "L2.(s & X)"]
    @test lagorder(r) == fill(2, H + 1)
    @test r.lags.rule == "ConstantLags(2)" && r.lags.p == 2
    # the same regression as the explicit formula, on the same sample
    e = panellp(@formula(leads(y) ~ s & X + lags(y, 1:2) + lags(s & X, 1:2) + fe(id) + fe(t)),
        d; id = :id, time = :t, horizon = H, shock = "s & X")
    for h in 0:H
        @test coef(r, h)[1] ≈ coef(e, h)[findfirst(==("s & X"), coefnames(e))] rtol = 1e-12
        @test nobs(r, h) == nobs(e, h)
    end
    # CappedLags (the authors' code): no lag at h = 0; the unused columns are NaN, not
    # reported as omitted
    c = fit(f, CappedLags(2))
    @test lagorder(c) == min.(0:H, 2)
    @test all(isnan, coef(c, 0)[2:end]) && all(isfinite, coef(c, 2))
    @test !c.fits[1].active[2] && c.fits[3].active[5]
    @test isempty(PanelLocalProjections._fitnotes(c))
    @test nobs(c, 0) == 200 * 40 && nobs(c, 1) == 200 * 38
    # default pmax: ⌈(T − H)^{1/3}⌉ = ⌈34^{1/3}⌉ = 4
    @test fit(f, CappedLags()).lags.pmax == 4
    @test lagorder(fit(f, CappedLags(0))) == zeros(Int, H + 1)
    # default pmax with the largest horizon: T = 40, H = 13 gives ⌈27^{1/3}⌉ = 3
    for rule in (CappedLags(), ICLags())
        @test panellp(f, d; id = :id, time = :t, horizon = 13,
            estimator = LagAugmented(rule)).lags.pmax == 3
    end
    # lags already in the formula are not added twice
    r1 = fit(@formula(leads(y) ~ s & X + lags(y, 1) + fe(id) + fe(t)), ConstantLags(2);
        shock = "s & X")
    @test count(==("L1.y"), coefnames(r1)) == 1 && "L2.y" in coefnames(r1)
    # … whatever their name: s & L1.X is the lag of the shock s & L0.X, so the regression is
    # s & lags(X, 0:3) + lags(y, 1) at every horizon
    shocks = ["s & L0.X", "s & L1.X", "s & L2.X"]
    r3 = fit(@formula(leads(y) ~ s & lags(X, 0:2) + fe(id) + fe(t)), ConstantLags(1);
        shock = shocks)
    @test coefnames(r3) == [shocks; "L1.y"; "L1.(s & L2.X)"]
    e3 = panellp(@formula(leads(y) ~ s & lags(X, 0:3) + lags(y, 1) + fe(id) + fe(t)), d;
        id = :id, time = :t, horizon = H, shock = shocks)
    @test all(f -> f.status === :ok, r3.fits)
    for nm in shocks
        @test coefpath(r3; term = nm) ≈ coefpath(e3; term = nm) rtol = 1e-10
    end
    # an added lag that the formula's columns span (L1.D.y = L1.y − L2.y) gives way to them,
    # without a warning
    r4 = fit(@formula(ldiff(y) ~ s & X + lags(y, 1:2) + fe(id) + fe(t)), CappedLags(3);
        shock = "s & X")
    jy = [findfirst(==(nm), coefnames(r4)) for nm in ("L1.y", "L2.y")]
    jd = findfirst(==("L1.D.y"), coefnames(r4))
    @test all(f -> all(f.keep[jy]) && !f.keep[jd], r4.fits)
    @test isempty(PanelLocalProjections._fitnotes(r4))
    # an aggregate shock with time effects is absorbed, whatever the rule
    for rule in (ICLags(), ConstantLags(1), CappedLags(2))
        @test_throws r"absorbed by the fixed effects" fit(
            @formula(leads(y) ~ X + fe(id) + fe(t)), rule)
    end
    # the lagged outcome of each response
    names(g) = coefnames(fit(g, ConstantLags(1)))
    @test "L1.D.y" in names(@formula(ldiff(y) ~ s & X + fe(id) + fe(t)))
    @test "L1.D.y" in names(@formula(leads(firstdiff(y)) ~ s & X + fe(id) + fe(t)))
    @test "L1.y" in names(@formula(cumul(y) ~ s & X + fe(id) + fe(t)))
    @test "L1.y" in names(@formula(anchor(y, y) ~ s & X + fe(id) + fe(t)))
    # requirements: unit effects, or finer ones nested in the units
    @test_throws ArgumentError fit(@formula(leads(y) ~ s & X + fe(t)), ConstantLags(1))
    dh = merge(d, (; half = Int.(d.t .> 20)))
    @test panellp(@formula(leads(y) ~ s & X + fe(id) & fe(half) + fe(t)), dh; id = :id,
        time = :t, horizon = 1, estimator = LagAugmented(ConstantLags(1))) isa
          PanelLocalProjection
    @test_throws ArgumentError panellp(@formula(y ~ s & X + fe(id) + fe(t)), d; id = :id,
        time = :t, horizon = 0, estimator = LagAugmented(ConstantLags(1)))
    @test_throws ArgumentError ConstantLags(-1)
    @test_throws ArgumentError CappedLags(-1)
    @test_throws ArgumentError ICLags(; pmax = 0)
    @test_throws ArgumentError ICLags(; criterion = :hq)
    # common sample, and the display
    rc = fit(f, ConstantLags(1); sample = :common)
    @test allequal(nobspath(rc))
    @test occursin("lags by horizon: 1, 1", sprint(show, MIME"text/plain"(), rc))
end

@testitem "ICLags: the synthetic information criterion" tags = [:core, :tlahr] setup = [PLPSupport, LSDV, TLAHRSupport] begin
    using .PLPSupport: fixture
    using .TLAHRSupport: ar2panel
    using LinearAlgebra, Statistics, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_tlahr)
    f = @formula(leads(y) ~ s & X + fe(id) + fe(t))
    r = panellp(f, d; id = :id, time = :t, horizon = 6, estimator = LagAugmented())
    sel = r.lags
    @test sel.pmax == 4 && sel.p == 1 && lagorder(r) == ones(Int, 7)
    # the criterion: the h = 0 regression of y on the shock and p lags of y; its residuals
    # averaged per period with the weights ŝ_i = s_i − mean(s) (balanced, time effects);
    # log of their mean square plus 2p/T, T the periods of the p = 0 regression
    T0 = 40
    pidx = PLP.PanelIndex(d.id, d.t)
    for p in 0:sel.pmax
        g = p == 0 ? f :
            PLP.leads(Term(:y)) ~ Term(:s) & Term(:X) + PLP.lags(Term(:y), p) + fe(:id) + fe(:t)
        e = panellp(g, d; id = :id, time = :t, horizon = 0, shock = "s & X")
        fit = e.fits[1]
        rows, M, _ = LSDV.within(e, fit, [d.id, d.t])
        col(nm) = nm == "s & X" ? d.s .* d.X :
                  PLP._lag(pidx, Float64.(d.y), parse(Int, nm[2:(findfirst('.', nm) - 1)]))
        Xr = reduce(hcat, [col(nm) for nm in coefnames(e)])
        ξ = M(d.y[rows]) - M(Xr[rows, :]) * fit.coef
        ŝ = d.s[rows] .- mean(d.s)
        tt = d.t[rows]
        ξbar = [sum(ŝ[tt .== q] .* ξ[tt .== q]) / sum(abs2, ŝ[tt .== q]) for q in sort(unique(tt))]
        @test sel.periods[p + 1] == length(ξbar) == T0 - p
        @test sel.ic[p + 1] ≈ log(mean(abs2, ξbar)) + 2p / T0 rtol = 1e-10
    end
    # BIC, common sample, outcome-and-shock lags
    b = panellp(f, d; id = :id, time = :t, horizon = 6, estimator = LagAugmented(ICLags(; criterion = :bic)))
    @test b.lags.p == 1
    @test b.lags.ic[2] ≈ log(b.lags.sigma2[2]) + log(T0) / T0 rtol = 1e-12
    c = panellp(f, d; id = :id, time = :t, horizon = 6,
        estimator = LagAugmented(ICLags(; sample = :common, vars = :outcome_and_shock, pmax = 3)))
    @test allequal(c.lags.periods) && c.lags.periods[1] == T0 - 3
    @test c.lags.ic[3] ≈ log(c.lags.sigma2[3]) + 2 * 2 * 2 / (T0 - 3) rtol = 1e-12
    # vars = :outcome_and_shock: candidate p adds p lags of the shock too, and counts them
    o = panellp(f, d; id = :id, time = :t, horizon = 6,
        estimator = LagAugmented(ICLags(; vars = :outcome_and_shock)))
    for p in 1:(o.lags.pmax)
        g = PLP.leads(Term(:y)) ~ Term(:s) & Term(:X) + PLP.lags(Term(:y), p) +
            PLP.lags(Term(:s) & Term(:X), p) + fe(:id) + fe(:t)
        sy = synthetic(panellp(g, d; id = :id, time = :t, horizon = 0, shock = "s & X"), 0)
        @test o.lags.sigma2[p + 1] ≈ sum(sy.weight .* sy.resid .^ 2) / sum(sy.weight) rtol = 1e-10
        @test o.lags.ic[p + 1] ≈ log(o.lags.sigma2[p + 1]) + 2p * 2 / T0 rtol = 1e-12
    end
    # an unbalanced panel: the periods weighted by their synthetic weight ω_t (fn 24)
    du = fixture(:fx_tlahr_unbal)
    ru = panellp(f, du; id = :id, time = :t, horizon = 6, estimator = LagAugmented())
    for p in 0:(ru.lags.pmax)
        g = p == 0 ? f :
            PLP.leads(Term(:y)) ~ Term(:s) & Term(:X) + PLP.lags(Term(:y), p) + fe(:id) + fe(:t)
        sy = synthetic(panellp(g, du; id = :id, time = :t, horizon = 0, shock = "s & X"), 0)
        ok = sy.weight .> 0
        @test ru.lags.periods[p + 1] == count(ok)
        @test ru.lags.sigma2[p + 1] ≈
              sum(sy.weight[ok] .* sy.resid[ok] .^ 2) / sum(sy.weight[ok]) rtol = 1e-9
    end
    @test_throws r"lower pmax" panellp(f, d; id = :id, time = :t, horizon = 3,
        estimator = LagAugmented(ICLags(; pmax = 39, sample = :common)))
    # a persistent AR(2): the criteria pick two lags, and p(h) = max(1, min(h, 2))
    a = ar2panel()
    for crit in (:aic, :bic)
        ra = panellp(f, a; id = :id, time = :t, horizon = 4,
            estimator = LagAugmented(ICLags(; criterion = crit)))
        @test ra.lags.p == 2
        @test lagorder(ra) == [1, 1, 2, 2, 2]
    end
    # lags of the outcome already in the formula are the candidates' own: p̂, p(h) and the
    # coefficients do not change (from h = 2 on with lags(y, 1:2), which enters at h = 0, 1)
    ra = panellp(f, a; id = :id, time = :t, horizon = 4, estimator = LagAugmented())
    for (g, same) in ((@formula(leads(y) ~ s & X + lags(y, 1) + fe(id) + fe(t)), 1:5),
        (@formula(leads(y) ~ s & X + lags(y, 1:2) + fe(id) + fe(t)), 3:5))
        rg = panellp(g, a; id = :id, time = :t, horizon = 4, estimator = LagAugmented(),
            shock = "s & X")
        @test rg.lags.p == 2 && lagorder(rg) == lagorder(ra)
        @test coefpath(rg)[same] ≈ coefpath(ra)[same] rtol = 1e-10
    end
    rb = panellp(@formula(leads(y) ~ s & X + lags(y, 1) + fe(id) + fe(t)), a; id = :id,
        time = :t, horizon = 4, estimator = LagAugmented(), shock = "s & X")
    @test rb.lags.ic[1] == rb.lags.ic[2]         # one regression: L1.y is in the formula
    @test rb.lags.ic[3] ≈ log(rb.lags.sigma2[3]) + 2 / rb.lags.periods[1] rtol = 1e-12
    # the criterion needs a shock that multiplies an aggregate series
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + fe(id)), fixture(:fx_balanced);
        id = :id, time = :t, horizon = 2, estimator = LagAugmented())
    @test panellp(@formula(leads(y) ~ x + fe(id)), fixture(:fx_balanced); id = :id, time = :t,
        horizon = 2, estimator = LagAugmented(ConstantLags(1))) isa PanelLocalProjection
end

@testitem "t-LAHR equals the brute-force formulas" tags = [:core, :tlahr] setup = [PLPSupport, LSDV, TLAHRSupport] begin
    using .PLPSupport: fixture
    using .TLAHRSupport: oracle
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    for (name, timefe) in ((:fx_tlahr, true), (:fx_tlahr_unbal, true), (:fx_tlahr_unbal, false))
        d = fixture(name)
        p = PLP.PanelIndex(d.id, d.t)
        lag(v, l) = PLP._lag(p, Float64.(v), l)
        sX, svX = d.s .* d.X, d.sv .* d.X
        # one shock, and two shocks of the same aggregate series
        specs = timefe ?
            [(@formula(leads(y) ~ s & X + fe(id) + fe(t)), ["s & X"], [d.s], [sX]),
             (@formula(leads(y) ~ s & X + sv & X + fe(id) + fe(t)), ["s & X", "sv & X"],
              [d.s, d.sv], [sX, svX])] :
            [(@formula(leads(y) ~ X + fe(id)), ["X"], [ones(length(d.id))], [Float64.(d.X)])]
        for (f, shock, exps, cols) in specs
            r = panellp(f, d; id = :id, time = :t, horizon = 3,
                estimator = LagAugmented(ConstantLags(1)), shock)
            V0, V2 = vcov(TLAHR(; refine = false), r), vcov(TLAHR(), r)
            @test V0.V == vcov(Clustered(:time; ssc = :none), r).V
            @test all(isinf, filter(!isnan, reduce(vcat, V0.df)))
            for h in 0:3
                fit = r.fits[h + 1]
                rows, M, _ = LSDV.within(r, fit, timefe ? [d.id, d.t] : [d.id])
                # columns in coefnames order: shocks, L1.y, lags of the shocks
                Xr = hcat(cols..., lag(d.y, 1), [lag(c, 1) for c in cols]...)
                @test coefnames(r) == [shock; "L1.y"; ["L1." * PLP._wrap(s) for s in shock]]
                Xw = M(Xr[rows, :])
                yw = M(PLP._lead(p, Float64.(d.y), h)[rows])
                o = oracle(Xw, yw, d.t[rows], collect(eachindex(shock)),
                    [e[rows] for e in exps], timefe)
                for x in o
                    @test coef(r, h)[x.j] ≈ x.β rtol = 1e-10
                    # eq. (11): CR0 by period
                    @test vcov(V0, h)[x.j, x.j] ≈ x.v0 rtol = 1e-9
                    @test vcov(V2, h)[x.j, x.j] ≈ x.v2 rtol = 1e-9
                    @test V2.df[h + 1][x.j] ≈ x.ν rtol = 1e-8
                end
                # several shocks: the joint matrix is positive semi-definite
                @test minimum(eigvals(Symmetric(vcov(V2, h)))) > -1e-12 * maximum(diag(vcov(V2, h)))
            end
        end
    end
    # eq. (11) through Frisch–Waugh–Lovell: σ² = Σ_t (Σ_i x̃_it ξ̂_it)² / (Σ x̃²)²
    d = fixture(:fx_tlahr_unbal)
    r = panellp(@formula(leads(y) ~ s & X + fe(id) + fe(t)), d; id = :id, time = :t,
        horizon = 2, estimator = LagAugmented(ConstantLags(2)))
    p = PLP.PanelIndex(d.id, d.t)
    fit = r.fits[3]
    rows, M, _ = LSDV.within(r, fit, [d.id, d.t])
    sX = d.s .* d.X
    Xr = hcat(sX, PLP._lag(p, Float64.(d.y), 1), PLP._lag(p, Float64.(d.y), 2),
        PLP._lag(p, sX, 1), PLP._lag(p, sX, 2))[rows, :]
    Xw = M(Xr)
    ξ = M(PLP._lead(p, Float64.(d.y), 2)[rows]) - Xw * fit.coef
    x̃ = Xw[:, 1] - Xw[:, 2:end] * (Xw[:, 2:end] \ Xw[:, 1])
    tt = d.t[rows]
    num = sum(abs2, [sum(x̃[tt .== q] .* ξ[tt .== q]) for q in unique(tt)])
    @test vcov(vcov(TLAHR(; refine = false), r), 2)[1, 1] ≈ num / sum(abs2, x̃)^2 rtol = 1e-9
end

@testitem "Imbens–Kolesár degrees of freedom: closed form and numerics" tags = [:core, :tlahr] begin
    using LinearAlgebra, StableRNGs
    const PLP = PanelLocalProjections
    # the closed form equals (Σλ)²/Σλ² for the eigenvalues of G'G, for every contrast e_k
    function eigν(X, k)
        H = X * pinv(X)
        d = 1 .- diag(H)
        a = X * pinv(X' * X)[:, k]
        G = (I - H) * Diagonal(a ./ sqrt.(d))
        λ = eigvals(Symmetric(G' * G))
        return sum(λ)^2 / sum(abs2, λ)
    end
    rng = StableRNG(5)
    X = randn(rng, 30, 4)
    des = PLP._ikdesign(X)
    for k in 1:4
        @test PLP._iknu(des, k) ≈ eigν(X, k) rtol = 1e-10
    end
    # the contrast matters: e₁ (the authors' code, for every coefficient) is not e₂
    @test !(PLP._iknu(des, 2) ≈ PLP._iknu(des, 1))
    # a period of high leverage (d_t < 0.01)
    Xh = randn(rng, 25, 3)
    Xh[1, :] .*= 40
    dh = PLP._ikdesign(Xh)
    @test minimum(dh.d) < 0.01
    for k in 1:3
        @test PLP._iknu(dh, k) ≈ eigν(Xh, k) rtol = 1e-10
    end
    # no T × T matrix: memory linear in T
    dl = PLP._ikdesign(randn(rng, 8000, 3))
    PLP._iknu(dl, 1)
    @test (@allocated PLP._iknu(dl, 2)) < 2^20
    # an exactly collinear column is dropped by the rank tolerance
    Xc = hcat(X, X[:, 1] .+ X[:, 2])
    dc = PLP._ikdesign(Xc)
    @test size(dc.U, 2) == 4
    @test dc.d ≈ des.d rtol = 1e-10
    # an ill-conditioned design, κ(X̄'X̄) ≈ 1e11, against a BigFloat reference
    t = range(0, 1; length = 40)
    Xi = hcat(ones(40), t, t .+ 3e-6 .* randn(rng, 40))
    @test 1e10 < cond(Xi' * Xi) < 1e12
    dI = PLP._ikdesign(Xi)
    setprecision(BigFloat, 256) do
        Xb = BigFloat.(Xi)
        Q = inv(Xb' * Xb)
        Hb = Xb * Q * Xb'
        db = 1 .- diag(Hb)
        @test dI.d ≈ Float64.(db) rtol = 1e-8
        for k in 1:3
            a = Xb * Q[:, k]
            bb = a .^ 2 ./ db
            Mb = (I - Hb) .^ 2
            νb = sum(abs2, a)^2 / (bb' * Mb * bb)
            @test PLP._iknu(dI, k) ≈ Float64(νb) rtol = 1e-6
        end
    end
    # leverage one (as many periods as regressors) is undefined
    @test PLP._ikdesign(randn(rng, 3, 3)) === nothing
    @test PLP._ikdesign([1.0 NaN; 2.0 3.0; 4.0 5.0]) === nothing
end

@testitem "t-LAHR: exposures and weights" tags = [:core, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, Statistics, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_tlahr)
    build(f; data = d, shock = nothing) = PLP._build_design(f, data; id = :id, time = :t,
        frequency = nothing, subset = nothing, clusters = (), shock)
    # exposure detection, for every column: the aggregate factor, whatever its position and
    # form (its series per period), and the exposure
    Xp = [d.X[findfirst(==(q), d.t)] for q in 1:40]
    e = build(@formula(leads(y) ~ X & s + fe(id)))
    @test e.aggseries[1] == Xp && e.exposures[1] == d.s
    e = build(@formula(leads(y) ~ log(abs(X) + 1) & sv & s + fe(id)))
    @test e.aggseries[1] ≈ log.(abs.(Xp) .+ 1) && e.exposures[1] ≈ d.sv .* d.s
    e = build(@formula(leads(y) ~ X + fe(id)))
    @test e.aggseries[1] == Xp && all(==(1), e.exposures[1])
    e = build(@formula(leads(y) ~ s & sv + fe(id)))
    @test isempty(e.aggseries[1]) && isempty(e.exposures[1]) && e.owner == [0]
    # lags of the aggregate factor, and of the whole interaction (exposure and series)
    e = build(@formula(leads(y) ~ s & lags(X, 0:1) + fe(id) + fe(t)); shock = ["s & L0.X", "s & L1.X"])
    @test e.exposures == [d.s, d.s] && e.aggseries[1] == Xp
    @test isequal(e.aggseries[2], [NaN; Xp[1:39]])
    @test !PLP._sameagg(e.aggseries[1], e.aggseries[2])
    e = build(@formula(leads(y) ~ lags(s & X, 0:1) + fe(id) + fe(t)); shock = "L0.(s & X)")
    @test e.exposures[1] == d.s && isequal(e.aggseries[2], [NaN; Xp[1:39]])
    @test isequal(e.exposures[2], PLP._lag(PLP.PanelIndex(d.id, d.t), d.s, 1))
    @test e.owner == [1, 1]       # the lag's exposure is the shock's
    # the same aggregate series up to scale; the intercept and unit-level columns have none
    e = build(@formula(leads(y) ~ s & X + sv & X2 + sv); data = merge(d, (; X2 = 2 .* d.X)),
        shock = "s & X")
    @test e.coefnames == ["(Intercept)", "sv", "s & X", "sv & X2"]
    @test isempty(e.exposures[1]) && isempty(e.exposures[2])
    @test PLP._sameagg(e.aggseries[3], e.aggseries[4])
    # owners: a shock its own synthetic regression, another exposure × aggregate its own,
    # the others the first shock's
    @test e.owner == [3, 3, 3, 4]
    @test build(@formula(leads(y) ~ s & X + fe(id) + fe(t))).timefe == [2]
    # two exposures of X (disjoint groups): each shock's weights are residualised on the
    # other exposure, so in a balanced panel its synthetic regression has no column for it
    dd = merge(d, (; s1 = d.s .* (d.id .<= 100), s2 = d.s .* (d.id .> 100)))
    r = panellp(@formula(leads(y) ~ s1 & X + s2 & X + fe(id) + fe(t)), dd; id = :id, time = :t,
        horizon = 2, estimator = LagAugmented(ConstantLags(1)), shock = ["s1 & X", "s2 & X"])
    cv = vcov(TLAHR(), r)
    for h in 0:2
        fit = r.fits[h + 1]
        X1 = fit.synthetic[1].X
        @test norm(X1[:, 2]) < 1e-10 * norm(X1[:, 1])         # s2 & X under s1's weights
        @test norm(X1[:, 5]) < 1e-10 * norm(X1[:, 4])         # L1.(s2 & X)
        @test all(isfinite, cv.df[h + 1]) && all(>(1), cv.df[h + 1])
    end
    # raw weights (the authors' code) keep the other exposure's columns
    @test norm(r.fits[1].synthetic[1].Xraw[:, 2]) > 1e-3 * norm(r.fits[1].synthetic[1].Xraw[:, 1])
    # with time effects and one time-invariant exposure, raw and residualised weights agree
    # on a balanced panel (the leverage is scale-invariant)
    r = panellp(@formula(leads(y) ~ s & X + fe(id) + fe(t)), d; id = :id, time = :t,
        horizon = 2, estimator = LagAugmented(ConstantLags(1)))
    @test vcov(TLAHR(; weights = :raw), r).V ≈ vcov(TLAHR(), r).V rtol = 1e-10
    @test vcov(TLAHR(; weights = :raw), r).df ≈ vcov(TLAHR(), r).df rtol = 1e-10
end

@testitem "t-LAHR: zero-weight periods, owners, and CR2 by period" tags = [:core, :tlahr] setup = [PLPSupport, TLAHRSupport] begin
    using .PLPSupport: fixture
    using .TLAHRSupport: cr2
    using LinearAlgebra, Logging, StatsModels, Test
    const PLP = PanelLocalProjections
    d = fixture(:fx_tlahr)
    kw = (; id = :id, time = :t)
    # one unit in period 1, two group dummies as exposures of X: the per-period
    # residualisation never solves a singular system; the period identifies g1 & X only
    sel = [(i <= 50 && t <= 30 && (t > 1 || i == 1)) for (i, t) in zip(d.id, d.t)]
    dg = map(c -> c[sel], d)
    dg = merge(dg, (; g1 = Float64.(dg.id .<= 25), g2 = Float64.(dg.id .> 25)))
    for est in (WithinFE(), LagAugmented(CappedLags(2)))
        r = panellp(@formula(leads(y) ~ g1 & X + g2 & X + fe(id)), dg; kw..., horizon = 2,
            shock = ["g1 & X", "g2 & X"], estimator = est)
        @test [s.weight[1] for s in r.fits[1].synthetic] == [1.0, 0.0]
        @test all(isfinite, vcov(TLAHR(), r).df[1][1:2])
    end
    # residualised exposures that are zero in exact arithmetic are zero, not round-off: two
    # units and two exposures in period 1 with time effects; and one exposure equal to 0.1
    # for every unit of period 20
    sel = [(t > 1 || i <= 2) for (i, t) in zip(d.id, d.t)]
    d2 = map(c -> c[sel], d)
    d2 = merge(d2, (; s2 = sin.(d2.id) .+ 0.3))
    r = panellp(@formula(leads(y) ~ s & X + s2 & X + fe(id) + fe(t)), d2; kw..., horizon = 2,
        estimator = LagAugmented(), shock = ["s & X", "s2 & X"])
    @test synthetic(r, 0).weight[1] == 0 && synthetic(r, 0).x[1] == 0
    @test all(<(2), r.lags.sigma2) && r.lags.p == 1
    @test all(isfinite, vcov(TLAHR(), r).df[1])
    for v in (ifelse.(d.t .== 20, 0.1, d.s), d.s .* (d.t .!= 20))
        dz = merge(d, (; sz = v))
        rz = panellp(@formula(leads(y) ~ sz & X + fe(id) + fe(t)), dz; kw..., horizon = 2,
            estimator = LagAugmented(ConstantLags(1)))
        cz = vcov(TLAHR(), rz)
        for h in 0:2
            fz = rz.fits[h + 1]
            kidx = findall(fz.keep)
            sy = fz.synthetic[1]
            k20 = findfirst(==(20), synthetic(rz, h).period)
            @test sy.weight[k20] == 0 && all(iszero, sy.X[k20, kidx])
            @test isnan(synthetic(rz, h).y[k20])
            # the period has no leverage: its scores enter the HC2 variance unadjusted
            des = PLP._ikdesign(sqrt.(sy.weight) .* sy.X[:, kidx])
            @test des.d[k20] ≈ 1 atol = 1e-14
            G = fz.scores.S[2][fz.timepresent, kidx] * fz.bread[kidx, kidx]
            @test vcov(cz, h)[1, 1] ≈ sum(G[:, 1] .^ 2 ./ des.d) rtol = 1e-12
            @test isfinite(cz.df[h + 1][1])
        end
        rl = panellp(@formula(leads(y) ~ sz & X + fe(id) + fe(t)), dz; kw..., horizon = 2,
            estimator = LagAugmented())
        @test rl.lags.periods == [39, 38, 37, 36, 35]      # period 20 has no weight
    end
    # labels do not matter: the weights of s & X are residualised on the exposure of X (1),
    # a control or a shock
    for p in (d, fixture(:fx_tlahr_unbal))
        fl = @formula(leads(y) ~ X + s & X + lags(X, 1) + fe(id))
        est = LagAugmented(ConstantLags(1))
        l1 = panellp(fl, p; kw..., horizon = 4, estimator = est, shock = "s & X")
        l2 = panellp(fl, p; kw..., horizon = 4, estimator = est, shock = ["X", "s & X"])
        jl = findfirst(==("s & X"), coefnames(l1))
        @test coefnames(l1) == coefnames(l2)
        cl1, cl2 = vcov(TLAHR(), l1), vcov(TLAHR(), l2)
        for h in 0:4
            @test cl1.df[h + 1][jl] ≈ cl2.df[h + 1][jl] rtol = 1e-12
            @test vcov(cl1, h)[jl, jl] ≈ vcov(cl2, h)[jl, jl] rtol = 1e-12
        end
    end
    # owners: s2 & L1.X has the exposure of the shock s2 & X, and uses its synthetic
    # regression, as the same lag written lags(s2 & X, 1)
    dd = merge(d, (; s1 = d.s, s2 = cos.(d.id)))
    f = @formula(leads(y) ~ s1 & X + s2 & X + s1 & lags(X, 1) + s2 & lags(X, 1) + lags(y, 1) +
                 fe(id) + fe(t))
    g = @formula(leads(y) ~ s1 & X + s2 & X + lags(s1 & X, 1) + lags(s2 & X, 1) + lags(y, 1) +
                 fe(id) + fe(t))
    rf = panellp(f, dd; kw..., horizon = 2, shock = ["s1 & X", "s2 & X"])
    rg = panellp(g, dd; kw..., horizon = 2, shock = ["s1 & X", "s2 & X"])
    @test coefnames(rf) == ["L1.y", "s1 & X", "s2 & X", "s1 & L1.X", "s2 & L1.X"]
    @test rf.owner == [2, 2, 3, 2, 3]
    jf, jg = 5, findfirst(==("L1.(s2 & X)"), coefnames(rg))
    cf, cg = vcov(TLAHR(), rf), vcov(TLAHR(), rg)
    for h in 0:2
        @test cf.df[h + 1][jf] ≈ cg.df[h + 1][jg] rtol = 1e-10
    end
    # shocks of different aggregate series with one exposure, s & lags(X, 0:1): not
    # residualised on each other, so the regression s & X + lags(s & X, 1)
    r = panellp(@formula(leads(y) ~ s & lags(X, 0:1) + fe(id) + fe(t)), d; kw..., horizon = 2,
        shock = ["s & L0.X", "s & L1.X"])
    rb = panellp(@formula(leads(y) ~ s & X + lags(s & X, 1) + fe(id) + fe(t)), d; kw...,
        horizon = 2, shock = "s & X")
    cv, cb = vcov(TLAHR(), r), vcov(TLAHR(), rb)
    ja, jb = findfirst(==("s & L0.X"), coefnames(r)), findfirst(==("s & X"), coefnames(rb))
    for h in 0:2
        @test cv.df[h + 1][ja] ≈ cb.df[h + 1][jb] rtol = 1e-10
        @test vcov(cv, h)[ja, ja] ≈ vcov(cb, h)[jb, jb] rtol = 1e-10
        @test r.fits[h + 1].synthetic[1].weight ≈ rb.fits[h + 1].synthetic[1].weight rtol = 1e-12
    end
    # two exposures in disjoint groups: each coefficient's degrees of freedom are those of
    # the scalar exposure (the regressors are all aggregate series times the exposures)
    dd = merge(d, (; s1 = d.s .* (d.id .<= 100), s2 = d.s .* (d.id .> 100)))
    r2 = panellp(@formula(leads(y) ~ s1 & X + s2 & X + lags(s1 & X, 1) + lags(s2 & X, 1) +
                          fe(id) + fe(t)), dd; kw..., horizon = 3, shock = ["s1 & X", "s2 & X"])
    r1 = panellp(@formula(leads(y) ~ s & X + lags(s & X, 1) + fe(id) + fe(t)), d; kw...,
        horizon = 3, shock = "s & X")
    v2, v1 = vcov(TLAHR(), r2), vcov(TLAHR(), r1)
    j1 = findfirst(==("s & X"), coefnames(r1))
    for h in 0:3, nm in ("s1 & X", "s2 & X")
        @test v2.df[h + 1][findfirst(==(nm), coefnames(r2))] ≈ v1.df[h + 1][j1] rtol = 1e-10
    end
    # the joint variance: each coefficient's scores over √d of its owner's synthetic
    # regression, the covariances included
    r = panellp(@formula(leads(y) ~ s & X + sv & X + fe(id) + fe(t)), d; kw..., horizon = 2,
        estimator = LagAugmented(ConstantLags(1)), shock = ["s & X", "sv & X"])
    @test coefnames(r) == ["s & X", "sv & X", "L1.y", "L1.(s & X)", "L1.(sv & X)"]
    @test r.owner == [1, 2, 1, 1, 2]
    cv = vcov(TLAHR(), r)
    for h in 0:2
        fh = r.fits[h + 1]
        kidx = findall(fh.keep)
        G = fh.scores.S[2][fh.timepresent, kidx] * fh.bread[kidx, kidx]
        des = [PLP._ikdesign(sqrt.(sy.weight) .* sy.X[:, kidx]) for sy in fh.synthetic]
        Z = reduce(hcat, [G[:, c] ./ sqrt.(des[r.owner[j]].d) for (c, j) in enumerate(kidx)])
        @test vcov(cv, h) ≈ Z' * Z rtol = 1e-12
        @test cv.df[h + 1][5] ≈ PLP._iknu(des[2], 5) rtol = 1e-12
    end
    # every regressor an aggregate series times the exposure, on an unbalanced panel: the
    # HC2 variance with the synthetic leverage of √ω_t X̄_t is CR2 clustered by period, and
    # the degrees of freedom its Bell–McCaffrey ones; the authors' unweighted X̄_t is not
    small = map(c -> c[fixture(:fx_tlahr_unbal).id .<= 25], fixture(:fx_tlahr_unbal))
    r = panellp(@formula(leads(y) ~ 0 + s & X + s & lags(X, 1)), small; kw..., horizon = 2,
        shock = ["s & X", "s & L1.X"])
    cv, craw = vcov(TLAHR(), r), vcov(TLAHR(; weights = :raw), r)
    for h in 0:2, k in 1:2
        o = cr2(small, h, k)
        @test coef(r, h)[k] ≈ o.β rtol = 1e-12
        @test sqrt(vcov(cv, h)[k, k]) ≈ o.se rtol = 1e-12
        @test cv.df[h + 1][k] ≈ o.ν rtol = 1e-12
        @test !(craw.df[h + 1][k] ≈ o.ν)
    end
    # type-stable
    ds = PLP._build_design(@formula(leads(y) ~ s & X + sv & X + fe(id) + fe(t)), d; kw...,
        frequency = nothing, subset = nothing, clusters = (), shock = ["s & X", "sv & X"])
    y1 = PLP._lead(ds.pidx, Float64.(d.y), 1)
    w = PLP._within_fit(ds, y1, isfinite.(y1), AlternatingProjections())
    tc = view(ds.groups[2], w.rows)
    @test (@inferred PLP._shockweights(ds, 1, w.rows, w.keep)) isa NTuple{2, Vector{Float64}}
    @test (@inferred PLP._synthetic(ds, w, tc, PLP._present(tc, ds.ngroups[2]))) isa
          Vector{PLP.SyntheticSeries}
end

@testitem "t-LAHR: what it refuses, and undefined variances" tags = [:core, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    b = fixture(:fx_balanced)
    # the refinement needs shocks that multiply an aggregate series
    r = panellp(@formula(leads(y) ~ x + fe(id)), b; id = :id, time = :t, horizon = 2,
        estimator = LagAugmented(ConstantLags(1)))
    @test_throws ArgumentError vcov(TLAHR(), r)
    @test vcov(TLAHR(; refine = false), r).V == vcov(Clustered(:time; ssc = :none), r).V
    # not for the split-panel jackknife; the cosine variance is
    s = panellp(@formula(leads(y) ~ x + fe(id)), b; id = :id, time = :t, horizon = 2,
        estimator = SplitPanelJackknife())
    @test_throws ArgumentError vcov(TLAHR(; refine = false), s)
    @test vcov(CosineHAR(), s) isa PanelLPVcov
    @test_throws ArgumentError TLAHR(; weights = :demeaned)
    # with unit effects the synthetic columns sum to zero over the periods, so every
    # leverage is at most 1 − 1/T; without them, as many periods as regressors give
    # leverage one: NaN
    d = fixture(:fx_tlahr)
    small = map(c -> c[d.t .<= 6], d)
    r = panellp(@formula(leads(y) ~ s & X + fe(id) + fe(t)), small; id = :id, time = :t,
        horizon = 1, estimator = LagAugmented(ConstantLags(2)))
    @test all(isfinite, vcov(TLAHR(), r).df[2][r.fits[2].keep])
    tiny = map(c -> c[d.t .<= 4], d)
    r = panellp(@formula(leads(y) ~ X + lags(X, 1)), tiny; id = :id, time = :t, horizon = 0,
        shock = "X")
    @test count(r.fits[1].timepresent) == 3 && count(r.fits[1].keep) == 3
    cv = @test_logs (:warn, r"leverage one") vcov(TLAHR(), r)
    @test all(isnan, cv.V[1]) && all(isnan, cv.df[1])
end

@testitem "CosineHAR equals the brute-force EWC" tags = [:core, :vcov, :tlahr] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    # harreg's rule ⌊0.41 T^{2/3}⌋ in double precision (T = 1000 gives 40, not 41)
    @test PLP._llsw(1000) == 40 && PLP._llsw(30) == 3 && PLP._llsw(38) == 4
    u = fixture(:fx_unbalanced)       # period 35 is absent from every sample
    r = panellp(@formula(leads(y) ~ x + lags(y, 1) + fe(id) + fe(t)), u; id = :id, time = :t,
        horizon = 2)
    h = 1
    fit = r.fits[h + 1]
    rows, M, _ = LSDV.within(r, fit, [u.id, u.t])
    p = PLP.PanelIndex(u.id, u.t)
    X = M(hcat(Float64.(u.x), PLP._lag(p, Float64.(u.y), 1))[rows, :])
    e = M(PLP._lead(p, Float64.(u.y), h)[rows]) - X * fit.coef
    tt = u.t[rows]
    t0, t1 = extrema(tt)
    G = zeros(t1 - t0 + 1, 2)
    for (k, q) in enumerate(tt)
        G[q - t0 + 1, :] .+= X[k, :] .* e[k]
    end
    T = size(G, 1)
    for (spec, B) in ((CosineHAR(), PLP._llsw(T)), (CosineHAR(7), 7))
        Φ = [sqrt(2) * cos(π * j * (t - 0.5) / T) for j in 1:B, t in 1:T]
        Ω = (Φ * G)' * (Φ * G) / B
        Bread = inv(X' * X)
        cv = vcov(spec, r)
        @test vcov(cv, h) ≈ Bread * Ω * Bread rtol = 1e-10
        @test cv.df[h + 1] == fill(Float64(B), 2)
        @test cv.info[h + 1].lags == B && cv.info[h + 1].periods == T
    end
    # B must be below the number of periods
    @test_logs (:warn, r"fewer cosine terms than periods") vcov(CosineHAR(500), r)
    @test_throws ArgumentError CosineHAR(0)
    @test_throws ArgumentError CosineHAR(:nw)
end

@testitem "synthetic series" tags = [:core, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using Dates, Statistics, StatsModels
    d = fixture(:fx_tlahr)
    # s ≡ 1 with unit effects: the cross-sectional means of the within-transformed data
    r = panellp(@formula(leads(y) ~ X + lags(X, 1:2) + fe(id)), d; id = :id, time = :t,
        horizon = 2, shock = "X")
    s = synthetic(r, 1)
    tt = 3:39
    Xs = d.X[1:40]
    @test s.period == collect(tt)
    @test s.x ≈ Xs[tt] .- mean(Xs[tt]) rtol = 1e-10
    @test all(==(200.0), s.weight)
    ybar = [mean(d.y[d.t .== q + 1]) for q in tt]
    @test s.y ≈ ybar .- mean(ybar) rtol = 1e-9
    @test abs(sum(s.resid)) < 1e-10
    @test synthetic(r, 1; weights = :raw).x ≈ s.x
    # a lag of the shock: its column in the shock's synthetic regression
    sl = synthetic(r, 1; term = "L1.X")
    @test sl.x ≈ Xs[tt .- 1] .- mean(Xs[tt .- 1]) rtol = 1e-10
    @test sl.weight == s.weight && sl.y == s.y
    @test_throws ArgumentError synthetic(r, 1; weights = :other)
    # not for the split-panel jackknife, nor for a column outside the horizon's regression
    sp = panellp(@formula(leads(y) ~ s & X + fe(id)), d; id = :id, time = :t, horizon = 1,
        estimator = SplitPanelJackknife())
    @test_throws r"store no synthetic series" synthetic(sp, 1)
    c = panellp(@formula(leads(y) ~ X + fe(id)), d; id = :id, time = :t, horizon = 2,
        estimator = LagAugmented(CappedLags(1)))
    @test_throws r"not in the regression" synthetic(c, 0; term = "L1.X")
    # Date periods
    q = fixture(:fx_dated)
    qd = merge(q, (; X = [Dates.month(x) / 3 + Dates.year(x) / 100 for x in q.date]))
    rq = panellp(@formula(leads(y) ~ X + fe(id)), qd; id = :id, time = :date,
        frequency = :quarter, horizon = 1)
    @test first(synthetic(rq, 0).period) == Date(2001, 1, 1)
end
