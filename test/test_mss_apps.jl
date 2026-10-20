# The four empirical applications of Mei, Sheng & Shi (2026): Romer–Romer (GDP and
# unemployment), Baron–Verner–Xiong, Mian–Sufi–Verner, Cerra–Saxena. The data are the
# authors' replication CSVs, kept locally under docs/references/ (not redistributed), so
# these items are tagged :local_mss and skipped when the data are absent. Each horizon is a
# pre-built outcome column, estimated with horizon = 0 as in the authors' notebook.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule MSSApps begin
    using PanelLocalProjections
    const DIR = joinpath(@__DIR__, "..", "docs", "references", "code", "panel-lp-replication",
        "applications")
    const LOCAL = joinpath(@__DIR__, "fixtures", "local", "mss_apps_R.csv")

    function readapp(file)
        lines = readlines(joinpath(DIR, file))
        hdr = Symbol.(strip.(split(lines[1], ","), '"'))
        cols = [Any[] for _ in hdr]
        for l in lines[2:end], (j, f) in enumerate(split(l, ","))
            f = strip(f, '"')
            push!(cols[j], f == "NA" || isempty(f) ? missing : parse(Float64, f))
        end
        return NamedTuple{Tuple(hdr)}(Tuple(any(ismissing, c) ? Vector{Union{Missing, Float64}}(c) :
                                            Vector{Float64}(c) for c in cols))
    end

    # file, id, time, shocks, controls, outcome columns, horizons, time FE, variance (as in
    # the notebook), outcome multiplier
    const APPS = Dict(
        "RR_gdp" => ("empirical_RR_f4_lngdp_1980.csv", :COUNTDUMS, :halfyear, ["CRISIS"],
            vcat(["l$(k)LNGDP" for k in 1:4], ["l$(k)CRISIS" for k in 1:4]),
            ["f$(h)LNGDP" for h in 0:10], 0:10, true, Clustered(:unit; ssc = :none), 1.0),
        "RR_unemp" => ("empirical_RR_f4_lnunemp_1980.csv", :COUNTDUMS, :halfyear, ["CRISIS"],
            vcat(["l$(k)UNEMP" for k in 1:4], ["l$(k)CRISIS" for k in 1:4]),
            ["f$(h)UNEMP" for h in 0:10], 0:10, true, Clustered(:unit; ssc = :none), 1.0),
        "BVX" => ("empirical_BVX_t1_y.csv", :ccode, :year, ["R_B", "R_N"],
            vcat(["L$(k)R_B" for k in 1:3], ["L$(k)R_N" for k in 1:3], "D1y",
                ["L$(k)D1y" for k in 1:3], "D1d_y", ["L$(k)D1d_y" for k in 1:3]),
            ["Fd$(h)y" for h in 1:6], 1:6, false,
            Clustered(:unit, :time; ssc = :mss, psd = :raw), 100.0),
        "MSV" => ("empirical_MSV_f2.csv", :CountryCode, :year, ["L0HHD_L1GDP", "L0NFD_L1GDP"],
            vcat(["L$(k)y" for k in 0:4], ["L$(k)HHD_L1GDP" for k in 1:4],
                ["L$(k)NFD_L1GDP" for k in 1:4]),
            ["F$(h)y" for h in 1:10], 1:10, false,
            Clustered(:unit, :time; ssc = :none, psd = false), 1.0),
        "CS" => ("empirical_CS_f3.csv", :cntyNUM, :obs, ["CRISIS"],
            vcat(["l$(k)CRISIS" for k in 1:4], ["l$(k)GRRT_WB" for k in 1:4]),
            ["cf$(h)GRRT_WB" for h in 1:10], 1:10, false, Clustered(:unit; ssc = :none), 1.0))

    "Coefficient, standard error and observations of every shock, horizon and estimator."
    function estimate_all()
        out = Dict{Tuple{String, String, Int, String}, NTuple{3, Float64}}()
        for (name, (file, id, tm, X, C, Ys, hs, te, vc, ymult)) in APPS
            d = readapp(file)
            f = eval(Meta.parse("@formula(Yh ~ " * join(X, " + ") * " + " * join(C, " + ") *
                                " + fe($id)" * (te ? " + fe($tm)" : "") * ")"))
            for (Y, h) in zip(Ys, hs), (est, lab) in ((WithinFE(), "FE"), (SplitPanelJackknife(), "SPJ"))
                dd = merge(d, (; Yh = ymult .* d[Symbol(Y)]))
                r = panellp(f, dd; id, time = tm, horizon = 0, estimator = est, singletons = :keep,
                    absorb = te ? OneSweep() : AlternatingProjections(), shock = X)
                cv = vcov(vc, r)
                for (j, x) in enumerate(X)
                    out[(name, lab, h, x)] = (coef(r, 0)[j], sqrt(cv.V[1][j, j]), nobs(r, 0))
                end
            end
        end
        return out
    end
end

@testitem "MSS applications: Table 1 and sample sizes" tags = [:local_mss, :spj] setup = [MSSApps] begin
    res = MSSApps.estimate_all()
    # Mei, Sheng & Shi (2026), Table 1: SPJ and FE at the horizon where SPJ peaks
    # (Romer–Romer scaled by 7, a moderate crisis)
    table1 = [("RR_gdp", 7, "CRISIS", -6.285, -5.432, 7), ("RR_unemp", 7, "CRISIS", 2.632, 2.128, 7),
        ("BVX", 4, "R_B", -4.213, -3.422, 1), ("BVX", 6, "R_N", -5.876, -4.561, 1),
        ("MSV", 7, "L0HHD_L1GDP", -0.558, -0.387, 1), ("MSV", 5, "L0NFD_L1GDP", -0.136, -0.124, 1),
        ("CS", 10, "CRISIS", -7.893, -4.507, 1)]
    for (app, h, x, spj, fe, scale) in table1
        @test round(scale * res[(app, "SPJ", h, x)][1]; digits = 3) == spj
        @test round(scale * res[(app, "FE", h, x)][1]; digits = 3) == fe
    end
    # complete-case observations per horizon (the authors keep singletons)
    nobs_ref = Dict(
        "RR_gdp" => [fill(1584, 7); 1560; 1536; 1512; 1488],
        "RR_unemp" => [fill(1331, 7); 1307; 1283; 1259; 1235],
        "BVX" => [2548, 2548, 2548, 2466, 2384, 2302],
        "MSV" => [752, 722, 692, 662, 632, 602, 572, 543, 515, 487],
        "CS" => [4121, 3961, 3808, 3650, 3502, 3357, 3215, 3078, 2941, 2807])
    for (app, ns) in nobs_ref
        hs = MSSApps.APPS[app][7]
        x = MSSApps.APPS[app][4][1]
        @test [Int(res[(app, "FE", h, x)][3]) for h in hs] == ns
    end
end

@testitem "MSS applications: the authors' R notebook calls" tags = [:local_mss, :spj, :r] setup = [MSSApps] begin
    # the authors' outputs on their data: with the data present, a missing output is a
    # failure (run `Rscript crosscheck/R/mss_apps_local.R`), not a silent pass
    @test isfile(MSSApps.LOCAL)
    if isfile(MSSApps.LOCAL)
        res = MSSApps.estimate_all()
        lines = readlines(MSSApps.LOCAL)
        @test length(lines) == 129
        for l in lines[2:end]
            app, method, h, term, b, se = split(l, ",")
            ours = res[(app, method, parse(Int, h), term)]
            @test ours[1] ≈ parse(Float64, b) rtol = 1e-11
            @test ours[2] ≈ parse(Float64, se) rtol = 1e-11
        end
    end
end

@testitem "MSS applications: exact two-way FE against xtlp tfe" tags = [:local_mss, :spj, :stata] setup = [MSSApps] begin
    # the Romer–Romer applications with exact two-way fixed effects (the default absorption)
    # against xtlp, whose tfe demeaning iterates: an independent check of the "exact"
    # columns of the one-sweep table (docs/src/replication/mss2026.md). The xtlp output is
    # local: run crosscheck/stata/mss_apps_local.do
    path = joinpath(dirname(MSSApps.LOCAL), "mss_apps_xtlp.csv")
    @test isfile(path)
    if isfile(path)
        lines = readlines(path)
        @test length(lines) == 45
        for app in ("RR_unemp", "RR_gdp")
            (file, id, tm, X, C, Ys, hs, te, vc, ymult) = MSSApps.APPS[app]
            d = MSSApps.readapp(file)
            f = eval(Meta.parse("@formula(Yh ~ " * join(X, " + ") * " + " * join(C, " + ") *
                                " + fe($id) + fe($tm))"))
            for l in lines[2:end]
                a, method, h, term, b, se, n = strip.(split(l, ","))
                a == app || continue
                dd = merge(d, (; Yh = d[Symbol(Ys[parse(Int, h) + 1])]))
                est = method == "FE" ? WithinFE() : SplitPanelJackknife()
                r = panellp(f, dd; id, time = tm, horizon = 0, estimator = est,
                    singletons = :keep, shock = X)
                @test coef(r, 0)[1] ≈ parse(Float64, b) rtol = 1e-10
                # xtlp's two-way demeaning stops on an absolute step of 1e-9
                @test sqrt(vcov(Clustered(:unit; ssc = :none), r).V[1][1, 1]) ≈
                      parse(Float64, se) rtol = 1e-8
                @test nobs(r, 0) == parse(Int, n)
            end
        end
    end
end
