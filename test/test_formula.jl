# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "regressors equal calendar lookups" tags = [:core, :formula] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    f = @formula(leads(y) ~ x + lags(x, 1:2) + lags(firstdiff(y), 1) + lags(s & X, 1) + s & X)
    d = PLP._build_design(f, u; id = :id, time = :t, frequency = nothing, subset = nothing,
        clusters = (), shock = "x")
    @test d.coefnames == ["(Intercept)", "x", "L1.x", "L2.x", "L1.D.y", "L1.(s & X)", "s & X"]
    @test d.intercept && d.shockcols == [2]
    look = Dict((u.id[r], u.t[r]) => r for r in eachindex(u.id))
    at(v, r, k) = (j = get(look, (u.id[r], u.t[r] + k), 0); j == 0 ? NaN : Float64(v[j]))
    same(a, b) = (isnan(a) && isnan(b)) || a == b
    for r in eachindex(u.id)
        @test same(d.X[r, 3], at(u.x, r, -1))
        @test same(d.X[r, 4], at(u.x, r, -2))
        @test same(d.X[r, 5], at(u.y, r, -1) - at(u.y, r, -2))
        @test same(d.X[r, 6], at(u.s, r, -1) * at(u.X, r, -1))
        @test same(d.X[r, 7], u.s[r] * u.X[r])
    end
    ctx = PLP._EvalContext(d.pidx, d.columns)
    y = Float64.(u.y)
    for h in 0:3
        lhs = Dict(
            :leads => (PLP.LeadsResponse(PLP.VarTerm(:y)), r -> at(y, r, h)),
            :ldiff => (PLP.LDiffResponse(PLP.VarTerm(:y)), r -> at(y, r, h) - at(y, r, -1)),
            :anchor => (PLP.AnchorResponse(PLP.VarTerm(:y), PLP.VarTerm(:y)), r -> at(y, r, h) - y[r]),
            :cumul => (PLP.CumulResponse(PLP.VarTerm(:y)), r -> sum(at(y, r, k) for k in 0:h)),
            :diff => (PLP.LeadsResponse(PLP.FirstDiffTerm(PLP.VarTerm(:y))),
                r -> at(y, r, h) - at(y, r, h - 1)))
        for (name, (resp, oracle)) in lhs
            v = PLP._lhs_at(resp, ctx, h)
            @test all(same(v[r], oracle(r)) for r in eachindex(v))
        end
    end
end

@testitem "formula parsing" tags = [:core, :formula] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    b = fixture(:fx_balanced)
    build(f; shock = nothing, kw...) = PLP._build_design(f, b; id = :id, time = :t,
        frequency = nothing, subset = nothing, clusters = (), shock = shock, kw...)
    d = build(@formula(leads(log(y + 20)) ~ 0 + x + w))
    @test !d.intercept && d.coefnames == ["x", "w"] && d.responsename == "leads(log(y + 20))"
    d = build(@formula(leads(y) ~ x + fe(id) + fe(grp) & fe(t)))
    @test !d.intercept && d.fenames == ["id", "grp&t"] && d.groupnames == [:unit, :time, Symbol("grp&t")]
    @test_throws ArgumentError build(@formula(leads(y) ~ x + fe(id) & w))           # FE slopes
    @test_throws ArgumentError build(@formula(leads(y) ~ x + s & X))                # shock required
    @test build(@formula(leads(y) ~ x + s & X); shock = "s & X").shockcols == [3]
    @test_throws ArgumentError build(@formula(leads(y) ~ x); shock = "nope")
    @test_throws ArgumentError build(@formula(leads(y) ~ x + (w ~ X)))              # IV (v0.2)
    @test_throws ArgumentError build(@formula(leads(y) ~ leads(x)))                 # RHS response
    @test_throws ArgumentError build(@formula(leads(y) ~ x + lags(x, 0:-1)))
    # same-named helper from an unknown module
    m = Module(:Elsewhere)
    Core.eval(m, :(leads(x) = x))
    f = StatsModels.FormulaTerm(StatsModels.FunctionTerm(m.leads, [StatsModels.term(:y)], :(leads(y))),
        (StatsModels.term(:x),))
    @test_throws ArgumentError build(f)
    # categorical regressors are not supported in v0.1
    bt = merge(b, (; c = string.(b.grp)))
    @test_throws ArgumentError PLP._build_design(@formula(leads(y) ~ x + c), bt; id = :id,
        time = :t, frequency = nothing, subset = nothing, clusters = (), shock = nothing)
    # a plain response is horizon-0 only
    @test_throws ArgumentError panellp(@formula(y ~ x + fe(id)), b; id = :id, time = :t, horizon = 2)
    @test panellp(@formula(y ~ x + fe(id)), b; id = :id, time = :t, horizon = 0) isa PanelLocalProjection
end

@testitem "subset is applied after leads and lags" tags = [:core, :formula] setup = [PLPSupport, LSDV] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + lags(y, 1) + fe(id))
    # a window that ends before the data: the leads of its last rows read t = 31..33
    keep = b.t .<= 30
    r = panellp(f, b; id = :id, time = :t, horizon = 3, subset = keep)
    @test nobspath(r) == fill(count(keep .& (b.t .>= 2)), 4)
    p = PLP.PanelIndex(b.id, b.t)
    X = hcat(Float64.(b.x), PLP._lag(p, Float64.(b.y), 1))
    for h in 0:3
        yh = PLP._lead(p, Float64.(b.y), h)
        rows = findall(keep .& isfinite.(yh) .& isfinite.(X[:, 2]))
        @test rows == findall(esample(r, h))
        _, M, _ = LSDV.within(r, r.fits[h + 1], [b.id])
        @test coef(r, h) ≈ M(X[rows, :]) \ M(yh[rows]) rtol = 1e-10
    end
    # pre-filtering instead loses the leads past the window
    rpre = panellp(f, map(v -> v[keep], b); id = :id, time = :t, horizon = 3)
    @test nobs(rpre, 3) == nobs(r, 3) - 3 * 30
    # a subset with missing values: missing is false
    km = Vector{Union{Missing, Bool}}(keep)
    km[1:5] .= missing
    @test esample(panellp(f, b; id = :id, time = :t, horizon = 0, subset = km), 0) ==
          esample(panellp(f, b; id = :id, time = :t, horizon = 0,
        subset = coalesce.(km, false)), 0)
end

@testitem "missing inputs of function terms" tags = [:core, :formula] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)                   # w has missing values
    @test any(ismissing, u.w)
    build(f) = PLP._build_design(f, u; id = :id, time = :t, frequency = nothing,
        subset = nothing, clusters = (), shock = "x")
    miss = ismissing.(u.w)
    # a comparison of a missing value is undefined, not false
    d = build(@formula(leads(y) ~ x + (w > 0)))
    @test all(isnan, d.X[miss, 3]) && !any(isnan, d.X[.!miss, 3])
    @test d.X[.!miss, 3] == Float64.(u.w[.!miss] .> 0)
    # nested calls keep their natural types (Bool conditions)
    d = build(@formula(leads(y) ~ x + ifelse(w > 0, w, 0.0) + !(w > 0) + ((w > 0) | (x > 0))))
    @test all(isnan, d.X[miss, 3:5])
    ok = .!miss
    @test d.X[ok, 3] == ifelse.(u.w[ok] .> 0, u.w[ok], 0.0)
    @test d.X[ok, 4] == Float64.(.!(u.w[ok] .> 0))
    @test d.X[ok, 5] == Float64.((u.w[ok] .> 0) .| (u.x[ok] .> 0))
    # a lag of a function of a missing value is undefined too
    d = build(@formula(leads(y) ~ x + lags(w > 0, 1)))
    p = d.pidx
    lagrow = PLP._shift(p, -1)
    @test all(r -> lagrow[r] == 0 || !miss[lagrow[r]] || isnan(d.X[r, 3]), eachindex(u.w))
    # so does a response that is a function of a missing value
    uz = merge(u, (; z = u.w))
    r = panellp(@formula(leads(z > 0) ~ x + fe(id)), uz; id = :id, time = :t, horizon = 1)
    @test !any(esample(r, 0)[miss])
end

@testitem "formula edge cases" tags = [:core, :formula] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    b = fixture(:fx_balanced)
    # StatsModels' lag/lead shift by row position: rejected
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + lag(y, 1) + fe(id)), b;
        id = :id, time = :t, horizon = 0)
    # numeric columns of abstract element type
    ba = merge(b, (; wa = Any[b.w...], wn = Vector{Union{Nothing, Float64}}(b.w)))
    r = panellp(@formula(leads(y) ~ x + w + fe(id)), b; id = :id, time = :t, horizon = 0)
    for v in (:wa, :wn)
        f = StatsModels.FormulaTerm(leads(term(:y)), (term(:x), term(v), fe(term(:id))))
        @test coef(panellp(f, ba; id = :id, time = :t, horizon = 0), 0) == coef(r, 0)
    end
    # NaN in a fixed-effect or cluster column is a missing value
    gnan = Float64.(b.grp)
    gnan[b.t .== 3] .= NaN
    bn = merge(b, (; g = gnan))
    rn = panellp(@formula(leads(y) ~ x + fe(id) + fe(g)), bn; id = :id, time = :t,
        horizon = 0, clusters = (:g,))
    @test !any(esample(rn, 0)[b.t .== 3]) && rn.fits[1].dof.felevels[2] == 5
    # the only regressor is an interaction: it is the shock
    rs = panellp(@formula(leads(y) ~ s & X + fe(id)), b; id = :id, time = :t, horizon = 0)
    @test rs.shocks == ["s & X"]
end
