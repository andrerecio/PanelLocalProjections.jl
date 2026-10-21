# Cross-checks of t-LAHR (Almuzara & Sancibrián 2026) on the fx_tlahr panels:
# - the authors' panel_LP.R, run on our panels (crosscheck/R/p4_as_panelLP.R →
#   test/fixtures/r/p4_as_panelLP.csv);
# - Stata (crosscheck/stata/p4_tlahr.do → test/fixtures/stata/p4_tlahr.csv): reghdfe for the
#   lag-augmented regression, and on the synthetic time-series regression, where the panel
#   regression is one, `regress …, vce(hc2, dfadjust)` (HC2 and Bell–McCaffrey df) and
#   `harreg …, estimator(ewc)`.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule CrossCheck4 begin
    using StatsModels, PanelLocalProjections
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    function readref(path)
        lines = readlines(path)
        hdr = Symbol.(split(lines[1], ","))
        map(lines[2:end]) do l
            NamedTuple{Tuple(hdr)}(Tuple(strip.(split(l, ","))))
        end
    end
    _num(s) = s == "." || isempty(s) ? NaN : parse(Float64, s)
    const DATA = Dict(p => read_csv(joinpath(FIX, "synthetic", "$p.csv"))
                      for p in ("fx_tlahr", "fx_tlahr_unbal"))
    const H = 6
    # panel_LP.R specs: right-hand side and shocks. Their regression always has the fixed
    # effects passed in FE, keeps singletons, and lags y and every s·X column min(h, p_max)
    # times: LagAugmented(CappedLags(p_max)).
    const AS = Dict(
        "s_twfe" => ("s & X + fe(id) + fe(t)", ["s & X"]),
        "ones_unit" => ("X + fe(id)", ["X"]),
        "sv_twfe" => ("sv & X + fe(id) + fe(t)", ["sv & X"]),
        "s_twfe_w" => ("s & X + sv + fe(id) + fe(t)", ["s & X"]),
        "s2_twfe" => ("s & X + sv & X + fe(id) + fe(t)", ["s & X", "sv & X"]))
    function asfit(panel, spec, pmax, cumulative)
        rhs, shock = AS[spec]
        resp = cumulative ? "cumul(y)" : "leads(y)"
        f = eval(Meta.parse("@formula($resp ~ $rhs)"))
        return panellp(f, DATA[panel]; id = :id, time = :t, horizon = H,
            estimator = LagAugmented(CappedLags(pmax)), shock, singletons = :keep)
    end
end

@testitem "R: Almuzara–Sancibrián panel_LP.R" tags = [:crosscheck, :r, :tlahr] setup = [CrossCheck4] begin
    using .CrossCheck4: readref, FIX, AS, H, asfit, _num
    const PLP = PanelLocalProjections
    ref = readref(joinpath(FIX, "r", "p4_as_panelLP.csv"))
    @test length(ref) == 672
    groups = unique((x.panel, x.spec, x.pmax, x.cumulative) for x in ref)
    @test length(groups) == 2 * length(AS) * 2 * 2
    for (panel, spec, pmax, cumul) in groups
        r = asfit(panel, spec, parse(Int, pmax), cumul == "1")
        cv = Dict((sm, w) => vcov(TLAHR(; refine = sm, weights = w), r)
                  for sm in (false, true), w in (:raw, :residual))
        # their code demeans with fixest (tolerance 1e-8): exact on the balanced panel
        tol = panel == "fx_tlahr" ? 1e-11 : 1e-7
        for x in filter(x -> (x.panel, x.spec, x.pmax, x.cumulative) == (panel, spec, pmax, cumul), ref)
            h, k, small = parse(Int, x.h), parse(Int, x.k), x.small == "1"
            j = findfirst(==(r.shocks[k]), coefnames(r))
            @test coef(r, h)[j] ≈ _num(x.estimate) rtol = tol
            c = cv[(small, :raw)]
            @test sqrt(vcov(c, h)[j, j]) ≈ _num(x.se) rtol = tol
            small || continue
            if spec == "s2_twfe" && k == 2
                # their degrees of freedom use the contrast e₁ for every coefficient
                # (panel_LP.R; bug B2 of docs/src/audit/as_code.md): ours use e₂
                @test !(c.df[h + 1][j] ≈ _num(x.df))
                fit = r.fits[h + 1]
                kidx = findall(fit.keep)
                sy = fit.synthetic[findfirst(s -> s.col == j, fit.synthetic)]
                @test PLP._iknu(PLP._ikdesign(sy.Xraw[:, kidx]), 1) ≈ _num(x.df) rtol = tol
            else
                @test c.df[h + 1][j] ≈ _num(x.df) rtol = tol
            end
            # the paper's weights (the demeaned exposure, periods weighted by ω_t) agree with
            # the code's (the raw exposure, unweighted) where ω_t is constant: on the balanced
            # panel, with a time-invariant exposure or s ≡ 1
            if panel == "fx_tlahr" && spec in ("ones_unit", "s_twfe", "s_twfe_w")
                cr = cv[(true, :residual)]
                @test sqrt(vcov(cr, h)[j, j]) ≈ _num(x.se) rtol = tol
                @test cr.df[h + 1][j] ≈ _num(x.df) rtol = tol
            end
        end
    end
end

@testitem "Stata: reghdfe, regress vce(hc2, dfadjust), harreg" tags = [:crosscheck, :stata, :tlahr] setup = [CrossCheck4] begin
    using .CrossCheck4: readref, FIX, DATA, H, _num
    using StatsModels
    ref = readref(joinpath(FIX, "stata", "p4_tlahr.csv"))
    @test length(ref) == 2 * (H + 1) + 4 * (H + 1)
    # A. the lag-augmented two-way regression, p(h) = min(h, 2), clustered by period (CR1,
    # reghdfe's factor; t(G − 1))
    for ds in ("fx_tlahr", "fx_tlahr_unbal")
        r = panellp(@formula(leads(y) ~ s & X + fe(id) + fe(t)), DATA[ds]; id = :id,
            time = :t, horizon = H, estimator = LagAugmented(CappedLags(2)))
        cv = vcov(Clustered(:time), r)
        rows = filter(x -> x.dataset == ds && x.spec == "la_cap2_twfe_cl_t", ref)
        @test length(rows) == H + 1
        for x in rows
            h = parse(Int, x.h)
            @test coef(r, h)[1] ≈ _num(x.b) rtol = 1e-10
            @test sqrt(vcov(cv, h)[1, 1]) ≈ _num(x.se) rtol = 1e-10
            @test nobs(r, h) == _num(x.N)
            @test cv.df[h + 1][1] == _num(x.df)
            @test cv.info[h + 1].groups == _num(x.G)
        end
    end
    # B. the synthetic regression: s ≡ 1 with unit effects, and ŝ_i with unit and time
    # effects; regressors the shock and two of its lags
    specs = Dict(
        "ones_unit" => (@formula(leads(y) ~ X + lags(X, 1:2) + fe(id)), "X"),
        "s_twfe" => (@formula(leads(y) ~ s & X + lags(s & X, 1:2) + fe(id) + fe(t)), "s & X"))
    for (spec, (f, shock)) in specs
        r = panellp(f, DATA["fx_tlahr"]; id = :id, time = :t, horizon = H, shock)
        j = findfirst(==(shock), coefnames(r))
        for (suffix, cv) in (("hc2", vcov(TLAHR(), r)), ("ewc", vcov(CosineHAR(), r)))
            rows = filter(x -> x.spec == "syn_$(spec)_$suffix", ref)
            @test length(rows) == H + 1
            for x in rows
                h = parse(Int, x.h)
                @test coef(r, h)[j] ≈ _num(x.b) rtol = 1e-10
                @test sqrt(vcov(cv, h)[j, j]) ≈ _num(x.se) rtol = 1e-10
                @test cv.df[h + 1][j] ≈ _num(x.df) rtol = 1e-10
                @test count(r.fits[h + 1].timepresent) == _num(x.N)
            end
        end
    end
end
