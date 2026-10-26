# Cross-checks of panel LP-IV (Phase 6a) on the fx_lpiv panels, against Stata
# (crosscheck/stata/p6_lpiv.do → test/fixtures/stata/p6_lpiv.csv): ivreghdfe and xtivreg2;
# and R (crosscheck/R/p6_fixest.R → test/fixtures/r/p6_fixest.csv): fixest's feols with an
# IV part, every variance without small-sample correction.
#
# Map (crosscheck/README.md):
# - ivreghdfe … robust (with or without `small`, which it ignores): Robust();
# - ivreghdfe … cluster(c): Clustered(c) when no fixed effect is nested in c. When one is
#   (unit effects with cluster(id), time effects with cluster(t)), ivreghdfe's factor
#   (N−1)/(N − K − df′)·G/(G−1) counts df′ = df_a minus the levels of the nested effect,
#   two fewer than reghdfe's rule (`ssc = :stata`); the test applies it to CR0;
# - ivreghdfe … cluster(id t): Clustered(:unit, :time; psd = false);
# - ivreghdfe … dkraay(b): DriscollKraaySE(b − 1) (b is the bandwidth, lags + 1);
# - xtivreg2 …, fe cluster(id): Clustered(:unit; ssc = :none); fe robust small: Robust().
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "Stata: ivreghdfe and xtivreg2" tags = [:crosscheck, :stata, :lpiv] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    num(s) = (s = strip(s); s == "." || isempty(s) ? NaN : parse(Float64, s))
    lines = readlines(joinpath(FIX, "stata", "p6_lpiv.csv"))
    hdr = Tuple(Symbol.(split(lines[1], ",")))
    ref = [NamedTuple{hdr}(Tuple(strip.(split(l, ",")))) for l in lines[2:end]]
    @test length(ref) == 160
    kw = (; id = :id, time = :t, horizon = 0:4)
    # ivreghdfe's clustered factor when the fixed effect `f` is nested in the cluster
    function ivfac(fit, f, G)
        n = fit.dof.nobs
        return (n - 1) / (n - fit.dof.k - (fit.dof.df_a - fit.dof.felevels[f])) * G / (G - 1)
    end
    for ds in ("fx_lpiv", "fx_lpiv_unbal")
        d = read_csv(joinpath(FIX, "synthetic", "$ds.csv"))
        r2 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id) + fe(t)), d; kw...)
        r1 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id)), d; kw...)
        rg = panellp(@formula(leads(y) ~ (g1 & ffr + g2 & ffr + g3 & ffr ~ g1 & mp + g2 & mp + g3 & mp) +
                w + fe(id)), d; kw...)
        cv = Dict(
            "robust" => vcov(Robust(), r2),
            "cl_id0" => vcov(Clustered(:unit; ssc = :none), r2),
            "cl_t0" => vcov(Clustered(:time; ssc = :none), r2),
            "cl_id_t" => vcov(Clustered(:unit, :time; psd = false), r2),
            "dk3" => vcov(DriscollKraaySE(2), r2),
            "fe_cl_t" => vcov(Clustered(:time), r1),
            "xt_cl_id" => vcov(Clustered(:unit; ssc = :none), r1),
            "xt_robust" => vcov(Robust(), r1),
            "grp" => vcov(Clustered(:time), rg))
        for x in filter(x -> x.dataset == ds, ref)
            h = parse(Int, x.h)
            spec = replace(x.spec, "_small" => "")
            r = spec == "grp_fe_cl_t" ? rg : startswith(spec, "twfe") ? r2 : r1
            j = spec == "grp_fe_cl_t" ? parse(Int, x.term[2:2]) : 1
            @test coef(r, h)[j] ≈ num(x.b) rtol = 1e-12
            fit = r.fits[h + 1]
            v = if spec == "twfe_robust"
                cv["robust"].V[h + 1][j, j]
            elseif spec == "twfe_cl_id"
                ivfac(fit, 1, parse(Int, x.G)) * cv["cl_id0"].V[h + 1][j, j]
            elseif spec == "twfe_cl_t"
                ivfac(fit, 2, parse(Int, x.G)) * cv["cl_t0"].V[h + 1][j, j]
            elseif spec == "twfe_cl_id_t"
                cv["cl_id_t"].V[h + 1][j, j]
            elseif spec == "twfe_dk3"
                cv["dk3"].V[h + 1][j, j]
            elseif spec == "fe_cl_t"
                cv["fe_cl_t"].V[h + 1][j, j]
            elseif spec == "xtivreg2_fe_cl_id"
                cv["xt_cl_id"].V[h + 1][j, j]
            elseif spec == "xtivreg2_fe_robust"
                cv["xt_robust"].V[h + 1][j, j]
            else
                @test spec == "grp_fe_cl_t"
                cv["grp"].V[h + 1][j, j]
            end
            @test sqrt(v) ≈ num(x.se) rtol = 1e-9
            @test fit.dof.nobs == parse(Int, x.N)
        end
    end
end

@testitem "R: fixest feols IV" tags = [:crosscheck, :r, :lpiv] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    lines = readlines(joinpath(FIX, "r", "p6_fixest.csv"))
    hdr = Tuple(Symbol.(split(lines[1], ",")))
    ref = [NamedTuple{hdr}(Tuple(strip.(split(l, ",")))) for l in lines[2:end]]
    @test length(ref) == 150
    spec = Dict("twfe_hetero0" => Robust(; ssc = :none),
        "twfe_cl_id0" => Clustered(:unit; ssc = :none),
        "twfe_cl_t0" => Clustered(:time; ssc = :none),
        "twfe_cl_id_t0" => Clustered(:unit, :time; ssc = :none, psd = false),
        # fixest lags over the periods present (period 25 is absent in fx_lpiv_unbal)
        "twfe_dk2_0" => DriscollKraaySE(2; ssc = :none, timelags = :position))
    term = Dict("fit_sffr" => "s & ffr", "w" => "w", "Ly" => "L1.y")
    for ds in ("fx_lpiv", "fx_lpiv_unbal")
        d = read_csv(joinpath(FIX, "synthetic", "$ds.csv"))
        # fixef.rm = "none" keeps singletons
        r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id) + fe(t)), d;
            id = :id, time = :t, horizon = 0:4, singletons = :keep)
        cv = Dict(k => vcov(v, r) for (k, v) in spec)
        for x in filter(x -> x.dataset == ds, ref)
            h = parse(Int, x.h)
            j = findfirst(==(term[x.term]), coefnames(r))
            # fixest demeans iteratively (exact on the balanced panel)
            tol = ds == "fx_lpiv" ? 1e-11 : 5e-9
            @test coef(r, h)[j] ≈ parse(Float64, x.b) rtol = tol
            @test sqrt(cv[x.spec].V[h + 1][j, j]) ≈ parse(Float64, x.se) rtol = tol
            @test nobs(r, h) == parse(Int, x.nobs)
        end
    end
end

# Lag-augmented panel LP-IV (Phase 6b): p = min(h, 2) lags of y, of s & ffr (endogenous)
# and of s & mp (instruments): LagAugmented(CappedLags(2)). Stata:
# crosscheck/stata/p6b_lpiv_la.do; R: crosscheck/R/p6b_fixest_la.R.
@testitem "Stata and R: lag-augmented LP-IV" tags = [:crosscheck, :stata, :r, :lpiv, :tlahr] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    function readref(path)
        lines = readlines(path)
        hdr = Tuple(Symbol.(split(lines[1], ",")))
        return [NamedTuple{hdr}(Tuple(strip.(split(l, ",")))) for l in lines[2:end]]
    end
    num(s) = s == "." || isempty(s) ? NaN : parse(Float64, s)
    st = readref(joinpath(FIX, "stata", "p6b_lpiv_la.csv"))
    fx = readref(joinpath(FIX, "r", "p6b_fixest_la.csv"))
    @test length(st) == 40 && length(fx) == 96
    est = LagAugmented(CappedLags(2))
    kw = (; id = :id, time = :t, horizon = 0:4, estimator = est)
    ivfac(fit, f, G) = (n = fit.dof.nobs;
        (n - 1) / (n - fit.dof.k - (fit.dof.df_a - fit.dof.felevels[f])) * G / (G - 1))
    term = Dict("fit_sffr" => "s & ffr", "fit_Lsffr1" => "L1.(s & ffr)",
        "fit_Lsffr2" => "L2.(s & ffr)", "w" => "w", "Ly1" => "L1.y", "Ly2" => "L2.y")
    for ds in ("fx_lpiv", "fx_lpiv_unbal")
        d = read_csv(joinpath(FIX, "synthetic", "$ds.csv"))
        r2 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...)
        r1 = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id)), d; kw...)
        @test lagorder(r2) == [0, 1, 2, 2, 2]
        cv0 = vcov(Clustered(:time; ssc = :none), r2)
        cv = Dict("la_twfe_robust" => vcov(Robust(), r2), "la_twfe_dk3" => vcov(DriscollKraaySE(2), r2),
            "la_fe_cl_t" => vcov(Clustered(:time), r1))
        for x in filter(x -> x.dataset == ds, st)
            h = parse(Int, x.h)
            r = x.spec == "la_fe_cl_t" ? r1 : r2
            fit = r.fits[h + 1]
            @test coef(r, h)[1] ≈ num(x.b) rtol = 1e-12
            @test fit.dof.nobs == parse(Int, x.N)
            v = x.spec == "la_twfe_cl_t" ? ivfac(fit, 2, parse(Int, x.G)) * cv0.V[h + 1][1, 1] :
                cv[x.spec].V[h + 1][1, 1]
            @test sqrt(v) ≈ num(x.se) rtol = 1e-9
        end
        # fixest keeps singletons; clustered by period without correction it is t-LAHR
        rk = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) + fe(t)), d; kw...,
            singletons = :keep)
        cvk = Dict("la_twfe_cl_t0" => vcov(TLAHR(refine = false), rk),
            "la_twfe_hetero0" => vcov(Robust(; ssc = :none), rk))
        tol = ds == "fx_lpiv" ? 1e-11 : 5e-9
        for x in filter(x -> x.dataset == ds, fx)
            h = parse(Int, x.h)
            j = findfirst(==(term[x.term]), coefnames(rk))
            @test coef(rk, h)[j] ≈ num(x.b) rtol = tol
            @test sqrt(cvk[x.spec].V[h + 1][j, j]) ≈ num(x.se) rtol = tol
            @test nobs(rk, h) == parse(Int, x.nobs)
        end
    end
end
