# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "summaries, accessors and display" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels, StatsFuns, Tables, DataFrames
    b = fixture(:fx_balanced)
    r = panellp(@formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id)), b;
        id = :id, time = :t, horizon = 4)
    @test horizons(r) == 0:4 && coefnames(r) == ["x", "L1.x", "L1.y", "w"]
    @test coefpath(r) == [coef(r, h)[1] for h in 0:4]
    @test coefpath(r; term = :w) == coefpath(r; term = "w")
    @test nobspath(r) == [nobs(r, h) for h in 0:4]
    @test count(esample(r, 2)) == nobs(r, 2)
    cv = vcov(Clustered(:unit), r)
    @test vcov(cv, 3) === cv.V[4]
    @test stderror(cv; term = "x") ≈ [sqrt(V[1, 1]) for V in cv.V]
    s = summarize(r, cv; level = 0.9, scale = 2)
    @test s isa IRFSummary && s.coef == 2 .* coefpath(r) && s.se ≈ 2 .* stderror(cv; term = "x")
    q = tdistinvcdf(29, 0.95)
    @test s.upper ≈ s.coef .+ q .* s.se
    sn = summarize(r, Clustered(:unit); dist = :normal)
    @test sn.upper ≈ sn.coef .+ norminvcdf(0.975) .* sn.se && all(isinf, sn.df)
    @test summarize(r, Clustered(:unit; ssc = :none)).df == fill(Inf, 5)
    df = DataFrame(s)
    @test names(df) == ["horizon", "coef", "se", "lower", "upper", "df", "nobs"] && nrow(df) == 5
    @test Tables.istable(s)
    for x in (r, cv, s)
        @test !isempty(sprint(show, MIME"text/plain"(), x))
        @test !isempty(sprint(show, x))
    end
    @test_throws ArgumentError coef(r, 9)
    @test_throws ArgumentError summarize(r, cv; term = "nope")
    @test_throws ArgumentError summarize(r, cv; level = 1.5)
    @test_throws ArgumentError vcov(Clustered(:nostored), r)
end

@testitem "a variance belongs to its result" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + lags(y, 1) + fe(id))
    r = panellp(f, b; id = :id, time = :t, horizon = 2)
    ra = panellp(f, b; id = :id, time = :t, horizon = 2, subset = b.t .<= 40)  # same shape
    rb = panellp(f, b; id = :id, time = :t, horizon = 4:6)
    @test summarize(r, vcov(Clustered(:unit), r)).se == summarize(r, Clustered(:unit)).se
    @test_throws ArgumentError summarize(r, vcov(Clustered(:unit), ra))
    @test_throws ArgumentError summarize(r, vcov(Clustered(:unit), rb))
end

@testitem "horizons and samples" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    u = fixture(:fx_unbalanced)
    f = @formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t))
    r06 = panellp(f, u; id = :id, time = :t, horizon = 6)
    r16 = panellp(f, u; id = :id, time = :t, horizon = 1:6)
    @test horizons(r16) == 1:6 && all(coef(r16, h) == coef(r06, h) for h in 1:6)
    @test horizons(panellp(f, u; id = :id, time = :t, horizon = [0, 1, 2])) == 0:2
    for bad in (0:2:6, [0, 2], 3.0, -1, Int[])
        @test_throws ArgumentError panellp(f, u; id = :id, time = :t, horizon = bad)
    end
    # the common sample: the same rows at every horizon, inside each horizon's own sample
    rc = panellp(f, u; id = :id, time = :t, horizon = 6, sample = :common)
    @test all(esample(rc, h) == esample(rc, 0) for h in 0:6)
    @test all(all(esample(rc, 0) .<= esample(r06, h)) for h in 0:6)
    # horizons beyond the data: empty (NaN), reported; none estimable: an error. At h = 58
    # each unit has one row (t = 2), a singleton of fe(id)
    b = fixture(:fx_balanced)
    fb = @formula(leads(y) ~ x + lags(y, 1) + fe(id))
    re = @test_logs (:warn, r"no observations at h = 58–60") panellp(fb, b; id = :id,
        time = :t, horizon = 57:60)
    @test [f.status for f in re.fits] == [:ok, :empty, :empty, :empty]
    @test re.fits[2].dof.singletons == 30
    @test all(isnan, summarize(re, Clustered(:unit)).se[2:4])
    @test_throws ArgumentError panellp(fb, b; id = :id, time = :t, horizon = 60:62)
    err = try
        panellp(fb, b; id = :id, time = :t, horizon = 60, sample = :common)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("horizon 59 has no observations", err.msg)
end

@testitem "Date periods through panellp" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using Dates, StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_dated)
    di = merge(d, (; q = [PLP._period_code(x, :quarter) for x in d.date]))
    f = @formula(leads(y) ~ x + lags(y, 1) + fe(id))
    rD = panellp(f, d; id = :id, time = :date, frequency = :quarter, horizon = 4)
    rI = panellp(f, di; id = :id, time = :q, horizon = 4)
    @test all(coef(rD, h) == coef(rI, h) for h in 0:4) && nobspath(rD) == nobspath(rI)
    @test vcov(DriscollKraaySE(2), rD).V == vcov(DriscollKraaySE(2), rI).V
    @test vcov(Clustered(:date), rD).V == vcov(Clustered(:time), rD).V
    @test coef(panellp(f, d; id = :id, time = :date, frequency = "quarter", horizon = 1), 1) ==
          coef(rD, 1)
    # dates on different days of the same quarter: time effects (and clusters) by quarter
    shifted = [isodd(i) ? lastdayofquarter(x) : x + Day(10) for (i, x) in zip(d.id, d.date)]
    ds = merge(d, (; date = shifted))
    ft = @formula(leads(y) ~ x + lags(y, 1) + fe(id) + fe(date))
    r1 = panellp(ft, ds; id = :id, time = :date, frequency = :quarter, horizon = 2,
        clusters = ((:id, :date),))
    r2 = panellp(@formula(leads(y) ~ x + lags(y, 1) + fe(id) + fe(q)), di; id = :id, time = :q,
        horizon = 2)
    @test all(coef(r1, h) ≈ coef(r2, h) for h in 0:2)
    @test r1.fits[1].dof.felevels == r2.fits[1].dof.felevels
end

@testitem "inputs" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels, Tables
    b = fixture(:fx_balanced)
    f = @formula(leads(y) ~ x + w + fe(id))
    r = panellp(f, b; id = :id, time = :t, horizon = 1)
    # row tables and generators of rows
    rows = Tables.rowtable(b)
    @test coef(panellp(f, rows; id = :id, time = :t, horizon = 1), 1) == coef(r, 1)
    @test coef(panellp(f, (x for x in rows); id = :id, time = :t, horizon = 1), 1) == coef(r, 1)
    # rows in any order
    perm = sortperm([hash((i, t)) for (i, t) in zip(b.id, b.t)])
    rp = panellp(f, map(v -> v[perm], b); id = :id, time = :t, horizon = 1)
    @test coef(rp, 1) ≈ coef(r, 1) rtol = 1e-12
    @test vcov(DriscollKraaySE(3), rp).V[2] ≈ vcov(DriscollKraaySE(3), r).V[2] rtol = 1e-12
    # columns of a different length than id and time (a column table does not check)
    short = merge(b, (; w = b.w[1:10]))
    @test_throws DimensionMismatch panellp(f, short; id = :id, time = :t, horizon = 0)
    @test_throws DimensionMismatch panellp(@formula(leads(y) ~ x + fe(g)), merge(b, (; g = [1, 2]));
        id = :id, time = :t, horizon = 0)
    # clusters as strings
    @test panellp(f, b; id = :id, time = :t, horizon = 0, clusters = "reg").groupnames ==
          [:unit, :time, :reg]
    @test_throws ArgumentError panellp(f, b; id = :id, time = :t, horizon = 0, clusters = 1)
    # :unit and :time name the panel dimensions
    bu = merge(b, (; unit = b.reg))
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + fe(unit)), bu; id = :id,
        time = :t, horizon = 0)
    @test_throws ArgumentError panellp(f, bu; id = :id, time = :t, horizon = 0, clusters = (:unit,))
    @test_throws ArgumentError panellp(f, merge(b, (; time = b.id)); id = :time, time = :t, horizon = 0)
    # formulas built programmatically
    fp = leads(term(:y)) ~ term(:x) + term(:w) + fe(term(:id))
    @test coef(panellp(fp, b; id = :id, time = :t, horizon = 1), 1) == coef(r, 1)
    p = 2
    fl = leads(term(:y)) ~ term(:x) + lags(term(:y), p) + fe(term(:id)) + fe(term(:grp)) & fe(term(:t))
    rl = panellp(fl, b; id = :id, time = :t, horizon = 0)
    @test coefnames(rl) == ["x", "L1.y", "L2.y"] && rl.fenames == ["id", "grp&t"]
    @test coef(rl, 0) == coef(panellp(@formula(leads(y) ~ x + lags(y, 2) + fe(id) + fe(grp) & fe(t)),
        b; id = :id, time = :t, horizon = 0), 0)
    @test sprint(show, fl) == "leads(y) ~ x + lags(y, 1:2) + fe(id) + fe(grp) & fe(t)"
end

@testitem "type stability of the hot paths" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels, Test
    const PLP = PanelLocalProjections
    b = fixture(:fx_unbalanced)
    r = panellp(@formula(leads(y) ~ x + lags(x, 1) + w + fe(id) + fe(t)), b; id = :id,
        time = :t, horizon = 1)
    rows = findall(r.fits[1].esample)
    fs = PLP._fe_sample([Int32.(b.id), Int32.(b.t)], rows)
    v = randn(length(rows))
    AP = AlternatingProjections()
    @test (@inferred PLP._absorb!(v, fs, AP, PLP._absorb_buffers(length(v), fs))) isa
          Tuple{Int, Bool}
    X = randn(50, 3)
    @test (@inferred PLP._ols(X, randn(50), X, false, false, [1])).status === :ok
    @test (@inferred PLP._score_stats(randn(50, 2), trues(2), [Int32.(rand(1:5, 50))], [5]))[2] ==
          [5]
    d = PLP._build_design(@formula(leads(y) ~ x + w + fe(id)), b; id = :id, time = :t,
        frequency = nothing, subset = nothing, clusters = (), shock = nothing)
    y = PLP._lhs_at(d.response, PLP._EvalContext(d.pidx, d.columns), 1)
    mask = d.rhsok .& isfinite.(y)
    @test (@inferred PLP._fit_horizon(WithinFE(), d, 1, y, mask, 0, AP)).status === :ok
    @test (@inferred PLP._fit_horizon(SplitPanelJackknife(), d, 1, y, mask, 0, AP)).status === :ok
    @test (@inferred vcov(Clustered(:unit, :time), r)) isa PanelLPVcov
    @test (@inferred vcov(DriscollKraaySE(), r)) isa PanelLPVcov
    @test (@inferred summarize(r, vcov(Robust(), r))) isa IRFSummary
end

@testitem "plot recipe" tags = [:core, :api] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels, RecipesBase
    b = fixture(:fx_balanced)
    r = panellp(@formula(leads(y) ~ x + fe(id)), b; id = :id, time = :t, horizon = 3)
    s = summarize(r, Clustered(:unit))
    series = RecipesBase.apply_recipe(Dict{Symbol, Any}(), s)
    @test length(series) == 2
    @test series[1].args[1] == 0:3 || series[1].args[1] == collect(0:3)
end
