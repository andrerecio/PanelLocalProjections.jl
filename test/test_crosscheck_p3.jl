# Cross-checks of cumulative-difference OLS (Ugarte-Ruiz 2026) against Stata (regress,
# xtscc, locproj; xtreg and xtlp m(spj) as comparators) and fixest, on the xtarsim panels
# and the native fx_cdols panel: test/fixtures/stata/p3_cdols.csv
# (crosscheck/stata/p3_cdols.do) and test/fixtures/r/p3_fixest.csv (crosscheck/R/p3_fixest.R).
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule CrossCheck3 begin
    using StatsModels, PanelLocalProjections
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    _num(s) = (s = strip(s, ['"', ' ']); s == "." || isempty(s) ? NaN : parse(Float64, s))
    function readref(path)
        lines = readlines(path)
        hdr = Symbol.(split(lines[1], ","))
        map(lines[2:end]) do l
            fields = [strip(m.match, ['"', ' ']) for m in eachmatch(r"\"[^\"]*\"|[^,]+", l)]
            NamedTuple{Tuple(hdr)}(Tuple(fields))
        end
    end
    const H = 9
    const DATASETS = ("xtarsim_s8y0", "xtarsim_s0y8", "xtarsim_s8y8", "fx_cdols", "fx_unbalanced")
    # fx_unbalanced: gaps, singleton units, a period absent for every unit; its unit-level
    # shock x is the CD-OLS shock
    function dataset(ds)
        d = read_csv(joinpath(FIX, startswith(ds, "fx_") ? "synthetic" : "stata", "$ds.csv"))
        return ds == "fx_unbalanced" ? merge(d, (; s = d.x)) : d
    end
    const DATA = Dict(ds => dataset(ds) for ds in DATASETS)
    # specs not run on fx_unbalanced (xtscc refuses a period absent for every unit)
    const BALANCED_ONLY = ("cdols_dk3", "xtlp_spj_level_cl_id")
    const TERMS = Dict("s" => "s", "L.s" => "L1.s", "L2.s" => "L2.s", "LD.y" => "L1.D.y",
        "l(s, 1)" => "L1.s")
    const CD = @formula(ldiff(y) ~ s + lags(s, 1))
    const DKX = DriscollKraaySE(3; ssc = :xtscc, timelags = :position)
    # spec => (formula, estimator, singletons, variance, number of terms)
    const STATA = Dict(
        "cdols_cl_id" => (CD, CDOLS(), :drop, Clustered(:unit), 2),
        "cdols_robust" => (CD, CDOLS(), :drop, Robust(), 2),
        "cdols_dy1_cl_id" => (@formula(ldiff(y) ~ s + lags(s, 1) + lags(firstdiff(y), 1)),
            CDOLS(), :drop, Clustered(:unit), 3),
        "cdols_s2_cl_id" => (@formula(ldiff(y) ~ s + lags(s, 1:2)), CDOLS(), :drop,
            Clustered(:unit), 3),
        # regress … i.t keeps singleton periods
        "cdols_tfe_cl_id" => (@formula(ldiff(y) ~ s + lags(s, 1) + fe(t)), CDOLS(), :keep,
            Clustered(:unit), 2),
        "cdols_tfe_robust" => (@formula(ldiff(y) ~ s + lags(s, 1) + fe(t)), CDOLS(), :keep,
            Robust(), 2),
        "cdols_dk3" => (CD, CDOLS(), :drop, DKX, 2),
        "locproj_cmlt_cl_id" => (CD, CDOLS(), :drop, Clustered(:unit), 1),
        "locproj_cmlt_yl1_cl_id" => (@formula(ldiff(y) ~ s + lags(s, 1) + lags(firstdiff(y), 1)),
            CDOLS(), :drop, Clustered(:unit), 1),
        "locproj_cmlt_robust" => (CD, CDOLS(), :drop, Robust(), 1),
        # comparators: fixed effects on the long difference (CMFE1, xtreg keeps singletons)
        # and the split-panel jackknife on the level response
        "cmfe1_cl_id" => (@formula(ldiff(y) ~ s + lags(s, 1) + fe(id)), WithinFE(), :keep,
            Clustered(:unit), 2),
        "xtlp_spj_level_cl_id" => (@formula(leads(y) ~ s + fe(id)), SplitPanelJackknife(),
            :drop, Clustered(:unit), 1))
end

@testitem "Stata: regress, xtscc, locproj (CD-OLS)" tags = [:crosscheck, :stata, :cdols] setup = [CrossCheck3] begin
    using .CrossCheck3: readref, FIX, DATA, DATASETS, STATA, TERMS, H, _num, BALANCED_ONLY
    ref = readref(joinpath(FIX, "stata", "p3_cdols.csv"))
    pairs = unique((x.dataset, x.spec) for x in ref)
    @test Set(pairs) == Set((ds, sp) for ds in DATASETS for sp in keys(STATA)
                            if !(ds == "fx_unbalanced" && sp in BALANCED_ONLY))
    for (ds, sp) in pairs
        f, est, sing, vc, nterms = STATA[sp]
        rows = filter(x -> x.dataset == ds && x.spec == sp, ref)
        @test length(rows) == nterms * (H + 1)
        r = panellp(f, DATA[ds]; id = :id, time = :t, horizon = H, estimator = est,
            singletons = sing)
        cv = vcov(vc, r)
        for x in rows
            k = parse(Int, x.h) + 1
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            @test sqrt(cv.V[k][j, j]) ≈ _num(x.se) rtol = 1e-10
            isnan(_num(x.N)) || @test nobs(r, k - 1) == Int(_num(x.N))
            isnan(_num(x.df_r)) || @test cv.df[k][j] == _num(x.df_r)
            isnan(_num(x.G)) || @test cv.info[k].groups == Int(_num(x.G))
        end
    end
end

@testitem "R: fixest (CD-OLS)" tags = [:crosscheck, :r, :cdols] setup = [CrossCheck3] begin
    using .CrossCheck3: readref, FIX, DATA, DATASETS, TERMS, H, CD, _num
    ref = readref(joinpath(FIX, "r", "p3_fixest.csv"))
    specs = Dict("fixest_cdols_cl_id" => Clustered(:unit), "fixest_cdols_hetero" => Robust(),
        "fixest_cdols_cl_id_none" => Clustered(:unit; ssc = :none))
    pairs = unique((x.dataset, x.spec) for x in ref)
    @test Set(pairs) == Set((ds, sp) for ds in DATASETS for sp in keys(specs))
    for (ds, sp) in pairs
        rows = filter(x -> x.dataset == ds && x.spec == sp, ref)
        @test length(rows) == 2 * (H + 1)
        r = panellp(CD, DATA[ds]; id = :id, time = :t, horizon = H, estimator = CDOLS())
        cv = vcov(specs[sp], r)
        for x in rows
            k = parse(Int, x.h) + 1
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            @test sqrt(cv.V[k][j, j]) ≈ _num(x.se) rtol = 1e-10
            @test nobs(r, k - 1) == Int(_num(x.nobs))
        end
    end
end
