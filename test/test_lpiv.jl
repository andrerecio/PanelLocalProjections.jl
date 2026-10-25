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
    @test occursin("lag augmentation", msg(f; estimator = LagAugmented()))
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
