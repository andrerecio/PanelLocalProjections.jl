# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

# Panel LP-IV (Phase 6a): two-stage least squares on the fixed-effects design.

@testmodule IVOracle begin
    using LinearAlgebra
    # 2SLS on within-transformed data: endogenous columns Xen, exogenous Xex, excluded
    # instruments Z (all already within-transformed). Columns of the result: [Xen Xex].
    function tsls(yt, Xen, Xex, Z)
        X = hcat(Xen, Xex)
        W = hcat(Z, Xex)
        Q = Matrix(qr(W).Q)[:, 1:size(W, 2)]
        Xh = Q * (Q' * X)
        β = (Xh' * Xh) \ (Xh' * yt)
        e = yt - X * β
        return (; β, B = inv(Xh' * Xh), U = Xh .* e, e)
    end
end

@testitem "2SLS equals the dummy-variable oracle, with every variance" tags = [:core, :lpiv] setup = [PLPSupport, LSDV, IVOracle] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    for (name, fes) in ((:fx_lpiv, :twoway), (:fx_lpiv_unbal, :twoway), (:fx_lpiv_unbal, :oneway))
        d = fixture(name)
        f = fes === :twoway ?
            @formula(leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id) + fe(t)) :
            @formula(leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id))
        r = panellp(f, d; id = :id, time = :t, horizon = 0:3)
        @test r.coefnames == ["s & ffr", "w", "L1.y"] && r.shocks == ["s & ffr"]
        @test r.endogenous == ["s & ffr"] && r.instruments == ["s & mp"]
        p = PLP.PanelIndex(d.id, d.t)
        for h in 0:3
            fit = r.fits[h + 1]
            rows, M, _ = LSDV.within(r, fit, fes === :twoway ? [d.id, d.t] : [d.id])
            o = IVOracle.tsls(M(PLP._lead(p, Float64.(d.y), h)[rows]),
                M(reshape(d.s .* d.ffr, :, 1)[rows, :]),
                M(hcat(d.w, PLP._lag(p, Float64.(d.y), 1))[rows, :]),
                M(reshape(d.s .* d.mp, :, 1)[rows, :]))
            @test fit.coef ≈ o.β rtol = 1e-10
            id, t = d.id[rows], d.t[rows]
            oracle = (Robust(; ssc = :none) => o.B * (o.U' * o.U) * o.B,
                Clustered(:unit; ssc = :none) => o.B * LSDV.meat(o.U, id) * o.B,
                Clustered(:time; ssc = :none) => o.B * LSDV.meat(o.U, t) * o.B,
                TLAHR(refine = false) => o.B * LSDV.meat(o.U, t) * o.B,
                Clustered(:unit, :time; ssc = :none, psd = false) =>
                    o.B * (LSDV.meat(o.U, id) + LSDV.meat(o.U, t) - o.U' * o.U) * o.B)
            for (spec, V) in oracle
                @test vcov(spec, r).V[h + 1] ≈ V rtol = 1e-9
            end
        end
    end
end

@testitem "IV identities: (x ~ x) is OLS; just-identified is (Z'X)⁻¹Z'y" tags = [:core, :lpiv] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_lpiv_unbal)
    kw = (; id = :id, time = :t, horizon = 0:4)
    ro = panellp(@formula(leads(y) ~ s & ffr + w + fe(id) + fe(t)), d; kw..., shock = "s & ffr")
    ri = panellp(@formula(leads(y) ~ (s & ffr ~ s & ffr) + w + fe(id) + fe(t)), d; kw...)
    pm = [findfirst(==(n), ro.coefnames) for n in ri.coefnames]
    for h in 0:4
        @test coef(ri, h) ≈ coef(ro, h)[pm] rtol = 1e-12
    end
    for spec in (Robust(), Clustered(:unit), Clustered(:time), Clustered(:unit, :time),
            DriscollKraaySE(), DriscollKraaySE(3; ssc = :none), CosineHAR(), TLAHR(refine = false))
        Vi, Vo = vcov(spec, ri), vcov(spec, ro)
        for k in eachindex(Vi.V)
            @test Vi.V[k] ≈ Vo.V[k][pm, pm] rtol = 1e-11
            @test Vi.df[k] == Vo.df[k][pm]
        end
    end
    # just-identified: β = (Z̃'X̃)⁻¹ Z̃'ỹ
    r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...)
    p = PLP.PanelIndex(d.id, d.t)
    fit = r.fits[3]
    rows, M, _ = LSDV.within(r, fit, [d.id, d.t])
    X = M(hcat(d.s .* d.ffr, d.w)[rows, :])
    Z = M(hcat(d.s .* d.mp, d.w)[rows, :])
    @test fit.coef ≈ (Z' * X) \ (Z' * M(PLP._lead(p, Float64.(d.y), 2)[rows])) rtol = 1e-10
end

@testitem "several endogenous regressors; collinear instruments left out" tags = [:core, :lpiv] setup = [PLPSupport, LSDV, IVOracle] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_lpiv)
    # group-specific responses to the instrumented aggregate (CFFS23-style)
    r = panellp(@formula(leads(y) ~ (g1 & ffr + g2 & ffr + g3 & ffr ~ g1 & mp + g2 & mp + g3 & mp) +
            w + fe(id)), d; id = :id, time = :t, horizon = 0:2)
    @test r.shocks == ["g1 & ffr", "g2 & ffr", "g3 & ffr"] && r.endogenous == r.shocks
    p = PLP.PanelIndex(d.id, d.t)
    G = hcat(d.g1, d.g2, d.g3)
    fit = r.fits[2]
    rows, M, _ = LSDV.within(r, fit, [d.id])
    o = IVOracle.tsls(M(PLP._lead(p, Float64.(d.y), 1)[rows]), M((G .* d.ffr)[rows, :]),
        M(reshape(Float64.(d.w), :, 1)[rows, :]), M((G .* d.mp)[rows, :]))
    @test fit.coef ≈ o.β rtol = 1e-10
    @test vcov(Clustered(:time; ssc = :none), r).V[2] ≈ o.B * LSDV.meat(o.U, d.t[rows]) * o.B rtol = 1e-9
    # a duplicated instrument (and one equal to an exogenous column) leaves the first stage,
    # with a note; the fit is the just-identified one. Of collinear instruments the later is
    # dropped (StatsModels lists main effects first: mp2, w2, s & mp).
    d2 = merge(d, (; mp2 = 2 .* d.s .* d.mp, w2 = 3 .* d.w))
    f2 = @formula(leads(y) ~ (s & ffr ~ s & mp + mp2 + w2) + w + fe(id))
    r2 = @test_logs (:warn, r"instruments left out of the first stage[^:]*: w2; s & mp") match_mode = :any panellp(
        f2, d2; id = :id, time = :t, horizon = 0:1)
    r1 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id)), d; id = :id, time = :t, horizon = 0:1)
    @test r2.instruments == ["mp2", "w2", "s & mp"]
    @test all(f -> f.instkeep == [true, false, false], r2.fits)
    for h in 0:1
        @test coef(r2, h) ≈ coef(r1, h) rtol = 1e-10
    end
    @test occursin("instrumented: s & ffr (2SLS); excluded instruments: mp2, w2, s & mp",
        sprint(show, MIME"text/plain"(), r2))
end

@testitem "IV specification errors" tags = [:core, :lpiv] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    d = fixture(:fx_lpiv)
    kw = (; id = :id, time = :t, horizon = 0)
    msg(f; k...) = try
        panellp(f, d; kw..., k...)
        ""
    catch err
        err isa ArgumentError ? sprint(showerror, err) : rethrow()
    end
    @test occursin("fe() terms cannot be inside", msg(@formula(leads(y) ~ (s & ffr ~ s & mp + fe(t)) + fe(id))))
    @test occursin("two instrumental-variable blocks", msg(@formula(leads(y) ~ (s & ffr ~ s & mp) + (w ~ sv) + fe(id))))
    @test occursin("is an exogenous regressor", msg(@formula(leads(y) ~ (s & ffr ~ s & mp + w) + w + fe(id))))
    @test occursin("underidentified", msg(@formula(leads(y) ~ (s & ffr + w ~ s & mp) + fe(id))))
    @test occursin("cannot be nested", msg(@formula(leads(y) ~ log(w ~ sv) + fe(id))))
    @test occursin("a constant cannot be", msg(@formula(leads(y) ~ (s & ffr ~ 1 + s & mp) + fe(id))))
    f = @formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id))
    @test occursin("split-panel jackknife", msg(f; estimator = SplitPanelJackknife()))
    @test occursin("no IV version", msg(@formula(ldiff(y) ~ (s & ffr ~ s & mp) + lags(s & ffr, 1)); estimator = CDOLS()))
    # an aggregate instrument with time effects is absorbed: underidentified everywhere
    @test occursin("underidentified at every horizon", msg(@formula(leads(y) ~ (s & ffr ~ mp) + w + fe(id) + fe(t))))
    r = panellp(f, d; kw...)
    @test_throws ArgumentError vcov(TLAHR(), r)
    @test vcov(TLAHR(refine = false), r).V[1] ≈ vcov(Clustered(:time; ssc = :none), r).V[1]
end

@testitem "over-identified 2SLS; collinear exogenous regressors; instrument tolerance" tags = [:core, :lpiv] setup = [PLPSupport, LSDV, IVOracle] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_lpiv_unbal)
    p = PLP.PanelIndex(d.id, d.t)
    kw = (; id = :id, time = :t, horizon = 0:3)
    # two instruments for one endogenous regressor: 2SLS, not IV with the first instrument
    r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp + sv & mp) + w + fe(id) + fe(t)), d; kw...)
    r1 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...)
    for h in 0:3
        fit = r.fits[h + 1]
        @test fit.instkeep == [true, true] && fit.dof.k == 2
        rows, M, _ = LSDV.within(r, fit, [d.id, d.t])
        o = IVOracle.tsls(M(PLP._lead(p, Float64.(d.y), h)[rows]),
            M(reshape(d.s .* d.ffr, :, 1)[rows, :]), M(reshape(Float64.(d.w), :, 1)[rows, :]),
            M(hcat(d.s .* d.mp, d.sv .* d.mp)[rows, :]))
        @test fit.coef ≈ o.β rtol = 1e-10
        @test !(coef(r1, h)[1] ≈ fit.coef[1])
        n, k, dfa = fit.dof.nobs, fit.dof.k, fit.dof.df_a
        # the default factor counts the second-stage regressors only (ivreg2, ivreghdfe, fixest)
        @test vcov(Robust(), r).V[h + 1] ≈ n / (n - k - dfa) .* (o.B * (o.U' * o.U) * o.B) rtol = 1e-9
        @test vcov(Clustered(:time; ssc = :none), r).V[h + 1] ≈ o.B * LSDV.meat(o.U, d.t[rows]) * o.B rtol = 1e-9
    end
    # of two collinear exogenous regressors one is omitted (which one is up to the pivoted
    # QR, as in OLS); the shock's estimate and variance are unchanged
    d3 = merge(d, (; w3 = 3 .* d.w))
    r3 = @test_logs (:warn, r"collinear with the other regressors") match_mode = :any panellp(
        @formula(leads(y) ~ (s & ffr ~ s & mp + sv & mp) + w + w3 + fe(id) + fe(t)), d3; kw...)
    for h in 0:3
        @test count(r3.fits[h + 1].keep) == 2
        @test coef(r3, h)[1] ≈ coef(r, h)[1] rtol = 1e-12
        @test vcov(Clustered(:time), r3).V[h + 1][1, 1] ≈ vcov(Clustered(:time), r).V[h + 1][1, 1] rtol = 1e-10
    end
    # an instrument close to, but not exactly, collinear with another stays in the first stage
    z2 = d.s .* d.mp .+ 1e-6 .* [sin(17.0 * i) for i in eachindex(d.mp)]
    rz = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp + z2) + w + fe(id) + fe(t)), merge(d, (; z2)); kw...)
    @test all(f -> f.instkeep == [true, true], rz.fits)
end

@testitem "underidentified at some horizons" tags = [:core, :lpiv] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    d0 = fixture(:fx_lpiv)
    # the aggregate mp (absorbed by the time effects) gets unit-level variation only in the
    # last two periods, which the leads of h ≥ 2 leave out
    z = d0.mp .+ ifelse.(d0.t .>= 49, d0.s .* d0.mp, 0.0)
    f = @formula(leads(y) ~ (s & ffr ~ z) + w + fe(id) + fe(t))
    r = @test_logs (:warn, r"underidentified at h = 2, 3") match_mode = :any panellp(f,
        merge(d0, (; z)); id = :id, time = :t, horizon = 0:3)
    @test [f.status for f in r.fits] == [:ok, :ok, :underidentified, :underidentified]
    @test all(isnan, coef(r, 2)) && all(isnan, coef(r, 3))
    @test occursin("underidentified at h = 2, 3", sprint(show, MIME"text/plain"(), r))
    for spec in (Robust(), Clustered(:unit, :time), DriscollKraaySE(2), TLAHR(refine = false))
        cv = vcov(spec, r)
        @test all(isfinite, cv.V[1]) && all(isnan, cv.V[3]) && all(isnan, cv.V[4])
    end
    s = summarize(r, Clustered(:time))
    @test isfinite(s.se[1]) && isnan(s.se[3])
end

@testitem "endogenous regressors the instruments do not predict are underidentified" tags = [:core, :lpiv] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    d0 = fixture(:fx_lpiv)
    n = length(d0.y)
    v = [sin(3.0 * i) for i in 1:n]
    x1 = d0.s .* d0.ffr
    # the pooled first stage: intercept, w, and the instruments s & mp, sv & mp
    M = hcat(ones(n), d0.w, d0.s .* d0.mp, d0.sv .* d0.mp)
    orth(u) = u .- M * (M \ u)
    kw = (; id = :id, time = :t, horizon = 0)
    msg(f, d; k...) = try
        panellp(f, d; kw..., k...)
        ""
    catch err
        err isa ArgumentError ? sprint(showerror, err) : rethrow()
    end
    # x2 orthogonal to the instruments and the exogenous columns: its fitted value is
    # round-off, which must not become an arbitrary coefficient (nor contaminate the shock)
    f = @formula(leads(y) ~ (s & ffr + x2 ~ s & mp + sv & mp) + w)
    @test occursin("underidentified at every horizon", msg(f, merge(d0, (; x2 = orth(v))); shock = "s & ffr"))
    # x2 whose fitted value is collinear with x1's: two instruments, one direction
    x2 = 2 .* (x1 .- orth(x1)) .+ orth(v)
    @test occursin("underidentified at every horizon", msg(f, merge(d0, (; x2)); shock = "s & ffr"))
    # a weak but genuine first stage is identified
    x2 = orth(v) .+ 1e-2 .* (d0.sv .* d0.mp)
    r = panellp(f, merge(d0, (; x2)); kw..., shock = "s & ffr")
    @test r.fits[1].status === :ok && all(isfinite, coef(r, 0))
    @test_throws ArgumentError synthetic(r, 0)
end

@testitem "lag-augmented LP-IV equals the hand-written regression" tags = [:core, :lpiv, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using LinearAlgebra, StatsModels
    d = fixture(:fx_lpiv_unbal)
    kw = (; id = :id, time = :t)
    f = @formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t))
    r = panellp(f, d; kw..., horizon = 0:3, estimator = LagAugmented(CappedLags(2)))
    @test lagorder(r) == [0, 1, 2, 2]
    @test r.coefnames == ["s & ffr", "w", "L1.y", "L2.y", "L1.(s & ffr)", "L2.(s & ffr)"]
    @test r.endogenous == ["s & ffr", "L1.(s & ffr)", "L2.(s & ffr)"]
    @test r.instruments == ["s & mp", "L1.(s & mp)", "L2.(s & mp)"] && r.instlag == [0, 1, 2]
    # horizon h uses the lags up to p(h) of the outcome, the endogenous shock and the instrument
    @test [f.instkeep for f in r.fits] == [[true, false, false], [true, true, false], [true, true, true], [true, true, true]]
    specs = (TLAHR(refine = false), Robust(), Clustered(:unit, :time), DriscollKraaySE(2), CosineHAR())
    for (h, p) in zip(0:3, lagorder(r))
        lh = p == 0 ? "" : " + lags(y, 1:$p)"
        endo = p == 0 ? "s & ffr" : "s & ffr + lags(s & ffr, 1:$p)"
        inst = p == 0 ? "s & mp" : "s & mp + lags(s & mp, 1:$p)"
        fh = eval(Meta.parse("@formula(leads(y) ~ ($endo ~ $inst) + w$lh + fe(id) + fe(t))"))
        rh = panellp(fh, d; kw..., horizon = h:h)
        pm = [findfirst(==(n), r.coefnames) for n in rh.coefnames]
        @test esample(r, h) == esample(rh, h)
        @test coef(r, h)[pm] ≈ coef(rh, h) rtol = 1e-12
        @test all(isnan, coef(r, h)[setdiff(eachindex(r.coefnames), pm)])
        for spec in specs
            @test vcov(spec, r).V[h + 1][pm, pm] ≈ vcov(spec, rh).V[1] rtol = 1e-11
        end
    end
    # t-LAHR for LP-IV is A&S's CR0 by period (supplement A.4); the refinement is not defined
    @test vcov(TLAHR(refine = false), r).V[3] ≈ vcov(Clustered(:time; ssc = :none), r).V[3]
    @test_throws ArgumentError vcov(TLAHR(), r)
    @test_throws ArgumentError synthetic(r, 1)
end

@testitem "ICLags on LP-IV is the criterion of the reduced form" tags = [:core, :lpiv, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    kw = (; id = :id, time = :t, horizon = 0:5)
    for ds in (:fx_lpiv, :fx_lpiv_unbal), rule in (ICLags(), ICLags(criterion = :bic, vars = :outcome_and_shock),
            ICLags(sample = :common, pmax = 3))
        d = fixture(ds)
        ri = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...,
            estimator = LagAugmented(rule))
        rr = panellp(@formula(leads(y) ~ s & mp + w + fe(id) + fe(t)), d; kw...,
            estimator = LagAugmented(rule), shock = "s & mp")
        @test ri.lags.p == rr.lags.p && lagorder(ri) == lagorder(rr)
        @test isequal(isfinite.(ri.lags.ic), isfinite.(rr.lags.ic))
        @test filter(isfinite, ri.lags.ic) ≈ filter(isfinite, rr.lags.ic) rtol = 1e-12
        @test ri.lags.periods == rr.lags.periods
    end
    # the reduced form puts the instruments in place of the endogenous regressors; the shock's
    # instrument (same exposure) is the shock, the other instruments are plain regressors
    d = fixture(:fx_lpiv)
    D = PanelLocalProjections._build_design(@formula(leads(y) ~ (s & ffr ~ s & mp + sv & mp) + w + fe(id)), d;
        id = :id, time = :t, frequency = nothing, subset = nothing, clusters = (), shock = nothing)
    R = PanelLocalProjections._reducedform(D)
    @test R.coefnames == ["s & mp", "sv & mp", "w"] && R.shockcols == [1] && isempty(R.endo)
    @test R.X[:, 3] == D.X[:, 2] && R.X[:, 1:2] == D.Z
end

@testitem "LP-IV Monte Carlo: t-LAHR coverage, IV against OLS" tags = [:mc, :lpiv, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: PLPDGPs
    using Logging, StableRNGs, Statistics, StatsModels
    # PLPDGPs.lpiv_mc: β(h) = 0.6^h; OLS on s & ffr is biased by about +0.4 at h = 0. With
    # T = 80 every estimator has a small-T bias at longer horizons (about −0.05 at h = 2–4,
    # shrinking with T, not with N), so unbiasedness is checked at h = 0, 1 only.
    R = 200
    hs = 0:4
    truth = 0.6 .^ hs
    rng = StableRNG(20261026)
    out = with_logger(NullLogger()) do
        map(1:R) do _
            d = PLPDGPs.lpiv_mc(rng)
            kw = (; id = :id, time = :t, horizon = hs, estimator = LagAugmented())
            ri = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + fe(id) + fe(t)), d; kw...)
            ro = panellp(@formula(leads(y) ~ s & ffr + fe(id) + fe(t)), d; kw..., shock = "s & ffr")
            covers(spec) = (s = summarize(ri, spec; level = 0.9);
                (s.lower .<= truth) .& (truth .<= s.upper))
            (; biv = coefpath(ri) .- truth, bols = coefpath(ro) .- truth,
                tlahr = covers(TLAHR(refine = false)), unit = covers(Clustered(:unit)))
        end
    end
    rate(k) = vec(sum(reduce(hcat, getfield.(out, k)); dims = 2)) ./ R
    med(k) = vec(median(reduce(hcat, getfield.(out, k)); dims = 2))
    # A&S Proposition 3: t-LAHR (supplement A.4) covers; unit clustering does not
    @test all(abs.(rate(:tlahr) .- 0.9) .<= 4 * sqrt(0.09 / R) + 0.02)
    @test all(rate(:unit) .<= 0.4)
    # the instrument removes the endogeneity bias of OLS
    @test all(abs.(med(:biv)[1:2]) .<= 0.03)
    @test med(:bols)[1] >= 0.3
end

@testitem "lag-augmented LP-IV: the shock's instrument weights ICLags; instrument lags follow p(h)" tags = [:core, :lpiv, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    d = fixture(:fx_lpiv_unbal)
    kw = (; id = :id, time = :t, horizon = 0:4)
    # two endogenous regressors, the shock s & ffr: the criterion is the reduced form's with
    # the instrument of the shock's exposure (s & mp), whatever the order of the instruments
    rr = panellp(@formula(leads(y) ~ s & mp + g1 & mp + fe(id) + fe(t)), d; kw...,
        estimator = LagAugmented(), shock = "s & mp")
    for f in (@formula(leads(y) ~ (s & ffr + g1 & ffr ~ s & mp + g1 & mp) + fe(id) + fe(t)),
        @formula(leads(y) ~ (s & ffr + g1 & ffr ~ g1 & mp + s & mp) + fe(id) + fe(t)))
        r = panellp(f, d; kw..., estimator = LagAugmented(), shock = "s & ffr")
        @test filter(isfinite, r.lags.ic) ≈ filter(isfinite, rr.lags.ic) rtol = 1e-12
    end
    # an endogenous non-shock whose instrument has no exposure structure does not matter
    r = panellp(@formula(leads(y) ~ (w + s & ffr ~ lags(w, 1) + s & mp) + fe(id) + fe(t)), d; kw...,
        estimator = LagAugmented(), shock = "s & ffr")
    @test r.lags.p >= 0 && all(f -> f.status === :ok, r.fits)
    # no instrument with the shock's exposure: the criterion is not defined (the fixed rules are)
    f = @formula(leads(y) ~ (s & ffr ~ sv & mp) + fe(id) + fe(t))
    @test_throws ArgumentError panellp(f, d; kw..., estimator = LagAugmented())
    @test lagorder(panellp(f, d; kw..., estimator = LagAugmented(CappedLags(2)))) == [0, 1, 2, 2, 2]
    # lag ℓ of the instruments enters only with lag ℓ of the endogenous shock: with
    # lags(s & ffr, 1) an exogenous control, no L1 instrument, and the fit at h depends on p(h)
    f = @formula(leads(y) ~ (s & ffr ~ s & mp) + lags(s & ffr, 1) + fe(id))
    r1 = panellp(f, d; kw..., estimator = LagAugmented(ConstantLags(1)))
    r2 = panellp(f, d; kw..., estimator = LagAugmented(CappedLags(2)))
    @test r2.instruments == ["s & mp", "L2.(s & mp)"]
    @test coef(r1, 1) ≈ filter(isfinite, coef(r2, 1)) rtol = 1e-12
end

@testitem "lag-augmented LP-IV: missing-value patterns, two instruments, notes" tags = [:core, :lpiv, :tlahr] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    d0 = fixture(:fx_lpiv_unbal)
    kw = (; id = :id, time = :t)
    est = LagAugmented(CappedLags(2))
    # hand-written regression of horizon h with p lags of y, of the endogenous shock and of
    # every instrument
    function hand(d, h, p, insts)
        lag(x) = p == 0 ? "" : " + lags($x, 1:$p)"
        inst = join((i * lag(i) for i in insts), " + ")
        f = eval(Meta.parse("@formula(leads(y) ~ (s & ffr$(lag("s & ffr")) ~ $inst) + w$(lag("y")) + fe(id) + fe(t))"))
        return panellp(f, d; kw..., horizon = h:h)
    end
    # an instrument observed from t = 7 on: the sample of each horizon needs the instrument
    # lags it uses; and an outcome observed from t = 4 on, on the balanced panel, where the
    # instrument lags are defined on the whole h = 0 sample but must not enter it (p(0) = 0)
    mpl = [t <= 6 ? NaN : v for (t, v) in zip(d0.t, d0.mp)]
    d1 = fixture(:fx_lpiv)
    yl = [t <= 3 ? NaN : v for (t, v) in zip(d1.t, d1.y)]
    for d in (merge(d0, (; mp = mpl)), merge(d1, (; y = yl)))
        r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...,
            horizon = 0:3, estimator = est)
        @test all(f -> f.status === :ok, r.fits)
        for (h, p) in zip(0:3, lagorder(r))
            rh = hand(d, h, p, ["s & mp"])
            pm = [findfirst(==(n), r.coefnames) for n in rh.coefnames]
            @test esample(r, h) == esample(rh, h)
            @test coef(r, h)[pm] ≈ coef(rh, h) rtol = 1e-12
            @test vcov(TLAHR(refine = false), r).V[h + 1][pm, pm] ≈ vcov(TLAHR(refine = false), rh).V[1] rtol = 1e-11
        end
        @test !any(n -> occursin("left out", n), PLP._fitnotes(r))
    end
    # two instruments: lag ℓ of each enters at p(h) ≥ ℓ
    r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp + sv & mp) + w + fe(id) + fe(t)), d0; kw...,
        horizon = 0:3, estimator = est)
    @test r.instruments == ["s & mp", "sv & mp", "L1.(s & mp)", "L2.(s & mp)", "L1.(sv & mp)", "L2.(sv & mp)"]
    @test [f.instkeep for f in r.fits] == [[1, 1, 0, 0, 0, 0], [1, 1, 1, 0, 1, 0], [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]
    for (h, p) in zip(0:3, lagorder(r))
        rh = hand(d0, h, p, ["s & mp", "sv & mp"])
        pm = [findfirst(==(n), r.coefnames) for n in rh.coefnames]
        @test coef(r, h)[pm] ≈ coef(rh, h) rtol = 1e-12
    end
    @test !any(n -> occursin("left out", n), PLP._fitnotes(r))
    # an instrument absorbed by the time effects is left out of the first stage, with a note,
    # and does not stop ICLags (it is not the shock's instrument in the reduced form)
    f = @formula(leads(y) ~ (s & ffr ~ s & mp + mp) + w + fe(id) + fe(t))
    ra = @test_logs (:warn, r"instruments left out of the first stage[^:]*: mp") match_mode = :any panellp(
        f, d0; kw..., horizon = 0:3, estimator = LagAugmented())
    rb = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d0; kw...,
        horizon = 0:3, estimator = LagAugmented())
    @test lagorder(ra) == lagorder(rb)
    for h in 0:3
        @test coef(ra, h) ≈ coef(rb, h) rtol = 1e-10
    end
end
