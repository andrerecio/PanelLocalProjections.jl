# Cross-checks of the fixed-effects panel LP against Stata (xtlp, reghdfe, xtreg, regress,
# xtscc) and R (fixest; the Mei–Sheng–Shi LP_panel). The reference outputs are committed
# fixtures computed on the synthetic panels: test/fixtures/stata/p1_fe.csv
# (crosscheck/stata/p1_fe.do), test/fixtures/r/p1_fixest.csv (crosscheck/R/p1_fixest.R) and
# test/fixtures/r/p1p2_mss_lppanel.csv (crosscheck/R/p1p2_mss_lppanel.R).
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule CrossCheck begin
    using StatsModels, PanelLocalProjections
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const PLP = PanelLocalProjections
    const FIX = joinpath(@__DIR__, "fixtures")

    _num(s) = s == "." || isempty(s) ? NaN : parse(Float64, s)
    # rows of a reference CSV; quoted fields may contain commas (fixest term names)
    function readref(path)
        lines = readlines(path)
        hdr = Symbol.(split(lines[1], ","))
        map(lines[2:end]) do l
            fields = [strip(m.match, ['"', ' ']) for m in eachmatch(r"\"[^\"]*\"|[^,]+", l)]
            NamedTuple{Tuple(hdr)}(Tuple(fields))
        end
    end

    const H = 6
    const RHS = "x + lags(x, 1) + lags(y, 1) + w"
    const NOLAG = "x + w"
    formula(lhs, fes; rhs = RHS) = eval(Meta.parse("@formula($lhs ~ $rhs" *
        join(" + fe($f)" for f in fes) * ")"))
    const TERMS = Dict("x" => "x", "L.x" => "L1.x", "L.y" => "L1.y", "w" => "w",
        "l(x, 1)" => "L1.x", "l(y, 1)" => "L1.y")

    # A Stata spec: fixed effects, singleton policy, variance, tolerance on SEs, whether the
    # critical df is compared, right-hand side.
    spec(fes, sing, vc, tol = 1e-12; df = true, rhs = RHS) = (; fes, sing, vc, tol, df, rhs)
    const ID, IDT = ["id"], ["id", "t"]
    DKP(m; kw...) = DriscollKraaySE(m; timelags = :position, kw...)
    const STATA = Dict(
        "xtlp_fe_cl_id" => spec(ID, :drop, Clustered(:unit)),
        "xtlp_fe_cl_id_ase" => spec(ID, :drop, Clustered(:unit; ssc = :none); df = false),
        "xtlp_fe_robust" => spec(ID, :drop, Robust()),
        "xtlp_fe_cl_t" => spec(ID, :drop, Clustered(:time)),
        "xtlp_fe_cl_idt" => spec(ID, :drop, Clustered(:unit, :time)),
        "xtlp_fe_cl_idt_ase" => spec(ID, :drop, Clustered(:unit, :time; ssc = :none); df = false),
        "xtlp_fe_dk" => spec(ID, :drop, DKP(:nw94)),
        "xtlp_fe_dk3" => spec(ID, :drop, DKP(3)),
        "xtlp_fe_dk3_nodfadj" => spec(ID, :drop, DKP(3; ssc = :xtscc)),
        "xtlp_fe_dk3_ase" => spec(ID, :drop, DKP(3; ssc = :none); df = false),
        "xtlp_fe_cl_id_keepsing" => spec(ID, :keep, Clustered(:unit)),
        # xtlp's two-way demeaning stops on an absolute step of 1e-9
        "xtlp_tfe_cl_id" => spec(IDT, :drop, Clustered(:unit), 1e-9),
        "xtlp_tfe_cl_id_ase" => spec(IDT, :drop, Clustered(:unit; ssc = :none), 1e-9; df = false),
        "xtlp_tfe_cl_t" => spec(IDT, :drop, Clustered(:time), 1e-9),
        "xtlp_tfe_cl_idt" => spec(IDT, :drop, Clustered(:unit, :time), 1e-9),
        "xtlp_tfe_cl_idt_ase" => spec(IDT, :drop, Clustered(:unit, :time; ssc = :none), 1e-9; df = false),
        "xtlp_tfe_robust" => spec(IDT, :drop, Robust(), 1e-9),
        "xtlp_tfe_dk" => spec(IDT, :drop, DKP(:nw94), 1e-9),
        "xtlp_tfe_dk3" => spec(IDT, :drop, DKP(3), 1e-9),
        "xtlp_tfe_dk3_nodfadj" => spec(IDT, :drop, DKP(3; ssc = :xtscc), 1e-9),
        "xtlp_tfe_cl_id_keepsing" => spec(IDT, :keep, Clustered(:unit), 1e-9),
        "reghdfe_idt_cl_idt" => spec(IDT, :drop, Clustered(:unit, :time), 1e-9),
        "reghdfe_idt_cl_id" => spec(IDT, :drop, Clustered(:unit), 1e-9),
        "reghdfe_id_cl_reg" => spec(ID, :drop, Clustered(:reg)),
        "reghdfe_3fe_cl_id" => spec(["id", "t", "grp) & fe(t"], :drop, Clustered(:unit), 1e-9),
        "reghdfe_id_dk4" => spec(ID, :drop, DKP(3)),
        "reghdfe_idt_dk4" => spec(IDT, :drop, DKP(3), 1e-9),
        "reghdfe_idt_robust" => spec(IDT, :drop, Robust(), 1e-9),
        "reghdfe_id_cl_regt" => spec(ID, :drop, Clustered(:reg, :time)),
        "reghdfe_noabs_cl_regt" => spec(String[], :drop, Clustered(:reg, :time)),
        "xtreg_fe_cl_id" => spec(ID, :keep, Clustered(:unit)),
        "xtreg_tfe_cl_id" => spec(IDT, :keep, Clustered(:unit), 1e-10),
        "regress_cl_id" => spec(String[], :drop, Clustered(:unit)),
        "xtscc_fe3" => spec(ID, :keep, DKP(3; ssc = :xtscc)),
        "xtlp_fe_nolag_cl_id" => spec(ID, :drop, Clustered(:unit); rhs = NOLAG),
        "xtlp_fe_nolag_cl_id_keepsing" => spec(ID, :keep, Clustered(:unit); rhs = NOLAG),
        "xtlp_tfe_nolag_cl_id" => spec(IDT, :drop, Clustered(:unit), 1e-9; rhs = NOLAG),
        "reghdfe_idt_nolag_cl_id" => spec(IDT, :drop, Clustered(:unit), 1e-9; rhs = NOLAG))
    # specs chosen because their two-way variance is not PSD on the fixtures
    const PSD_TRIGGERS = ("reghdfe_id_cl_regt", "reghdfe_noabs_cl_regt")
    const TRANSFORMS = Dict("level" => "leads(y)", "diff" => "leads(firstdiff(y))",
        "cmltdiff" => "ldiff(y)", "cmltsum" => "cumul(y)")

    # The synthetic panels and the panels derived from them in crosscheck/stata/p1_fe.do
    # and crosscheck/R/p1_fixest.R, by the same filters.
    const DATA = let d = Dict{String, Any}(name => read_csv(joinpath(FIX, "synthetic", "$name.csv"))
                                           for name in ("fx_balanced", "fx_unbalanced"))
        sub(nt, keep) = map(v -> v[keep], nt)
        b = d["fx_balanced"]
        for name in ("fx_balanced", "fx_unbalanced")
            d[name * "_t29"] = sub(d[name], d[name].t .<= 29)
        end
        g = sub(b, .!((b.id .== 30) .& .!in.(b.t, Ref((10, 11)))))
        d["fx_balanced_gap"] = merge(g, (; w = [t == 20 ? missing : w for (t, w) in zip(g.t, g.w)]))
        d["fx_balanced_2blocks"] = sub(b, .!(((b.id .<= 15) .& (b.t .> 30)) .|
                                               ((b.id .> 15) .& (b.t .<= 30))))
        d["fx_balanced_id4"] = sub(b, b.id .<= 4)
        d
    end
    estimate(ds, lhs, fes; singletons = :drop, absorb = AlternatingProjections(),
        estimator = WithinFE(), rhs = RHS) =
        panellp(formula(lhs, fes; rhs), DATA[ds]; id = :id, time = :t, horizon = H,
            singletons, absorb, estimator, clusters = (:reg,))
end

@testitem "Stata: xtlp, reghdfe, xtreg, regress, xtscc" tags = [:crosscheck, :stata, :vcov] setup = [CrossCheck] begin
    using .CrossCheck: readref, estimate, STATA, TERMS, TRANSFORMS, FIX, H, RHS, PSD_TRIGGERS, _num
    const PLP = PanelLocalProjections
    ref = readref(joinpath(FIX, "stata", "p1_fe.csv"))
    pairs = unique((x.dataset, x.spec) for x in ref)
    # every reference spec is mapped, and every mapping has reference rows
    @test Set(last.(pairs)) ==
          union(Set(keys(STATA)), Set("xtlp_fe_tr_$tr" for tr in keys(TRANSFORMS)))
    for (ds, sp) in pairs
        startswith(sp, "xtlp_fe_tr_") && continue
        s = STATA[sp]
        rows = filter(x -> x.dataset == ds && x.spec == sp, ref)
        @test length(rows) == (s.rhs == RHS ? 4 : 2) * (H + 1)
        r = estimate(ds, "leads(y)", s.fes; singletons = s.sing, rhs = s.rhs)
        cv = vcov(s.vc, r)
        clipped = false     # a horizon where the PSD repair fired and the SEs agree
        for x in rows
            k = parse(Int, x.h) + 1
            fit = r.fits[k]
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            if sp == "reghdfe_noabs_cl_regt" && _num(x.se) == 0
                # reghdfe's repaired V (slope block replaced) is refused by Stata when
                # posted ("variance matrix is nonsymmetric or highly singular") and zeroed
                @test cv.info[k].psdfixed
                continue
            end
            @test sqrt(cv.V[k][j, j]) ≈ _num(x.se) rtol = s.tol
            clipped |= cv.info[k].psdfixed
            @test nobs(r, k - 1) == Int(_num(x.N))
            s.df && @test cv.df[k][j] == _num(x.df_r)
            G = _num(x.G)    # clusters, or periods for reghdfe's Driscoll–Kraay
            isnan(G) || @test G == (s.vc isa DriscollKraaySE ? cv.info[k].periods : cv.info[k].groups)
            if startswith(sp, "reghdfe")
                # reghdfe counts the constant of a pooled regression as an absorbed level
                gs = s.vc isa Robust ? Int[] :
                     s.vc isa DriscollKraaySE ? [2] : [PLP._dimindex(r, d) for d in s.vc.dims]
                dfa = isempty(gs) ? fit.dof.df_a : PLP._dfa_nested(fit, gs)[1]
                @test dfa + r.intercept == Int(_num(x.df_a))
                @test sum(fit.dof.felevels) + r.intercept == Int(_num(x.df_a_initial))
            end
        end
        sp in PSD_TRIGGERS && @test clipped
    end
    # response transforms built by xtlp itself (hor(0 6) ytransf(…))
    for (tr, lhs) in TRANSFORMS, ds in ("fx_balanced", "fx_unbalanced")
        r = estimate(ds, lhs, ["id"])
        cv = vcov(Clustered(:unit), r)
        rows = filter(x -> x.dataset == ds && x.spec == "xtlp_fe_tr_$tr", ref)
        @test length(rows) == 4 * (H + 1)
        for x in rows
            k = parse(Int, x.h) + 1
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            @test sqrt(cv.V[k][j, j]) ≈ _num(x.se) rtol = 1e-12
            @test cv.df[k][j] == _num(x.df_r)
        end
    end
    # the derived panels exercise what the plain fixtures cannot
    dk = estimate("fx_balanced_t29", "leads(y)", ["id"])
    @test [i.lags for i in vcov(STATA["xtlp_fe_dk"].vc, dk).info] == [3; fill(2, H)]
    gap = estimate("fx_balanced_gap", "leads(y)", ["id"]; singletons = :keep)
    @test gap.fits[1].dof.nobs == estimate("fx_balanced_gap", "leads(y)", ["id"]).fits[1].dof.nobs + 1
    @test !(vcov(DriscollKraaySE(3; timelags = :position), gap).V[1] ≈
            vcov(DriscollKraaySE(3), gap).V[1])
    blocks = estimate("fx_balanced_2blocks", "leads(y)", ["id", "t"])
    @test blocks.fits[1].dof.fecomponents[1, 2] == 2
end

@testitem "R: fixest" tags = [:crosscheck, :r, :vcov] setup = [CrossCheck] begin
    using .CrossCheck: readref, estimate, TERMS, FIX, H, _num
    ref = readref(joinpath(FIX, "r", "p1_fixest.csv"))
    # fixest's default vcov_fix = TRUE clips negative eigenvalues of V itself: psd = :raw
    specs = Dict(
        "fixest_fe_cl_id_none" => (["id"], Clustered(:unit; ssc = :none)),
        "fixest_fe_hetero_none" => (["id"], Robust(; ssc = :none)),
        "fixest_tfe_cl_id_none" => (["id", "t"], Clustered(:unit; ssc = :none)),
        "fixest_tfe_twoway_none" => (["id", "t"], Clustered(:unit, :time; ssc = :none, psd = :raw)),
        # fixest 0.14's default small-sample correction coincides with Stata's here
        "fixest_tfe_cl_id_default" => (["id", "t"], Clustered(:unit)),
        # four units: the two-way variance is not PSD (with and without fixest's repair)
        "fixest_fe_twoway_none_id4" => (["id"], Clustered(:unit, :time; ssc = :none, psd = :raw)),
        "fixest_fe_twoway_none_nofix_id4" => (["id"], Clustered(:unit, :time; ssc = :none, psd = false)))
    pairs = unique((x.dataset, x.spec) for x in ref)
    @test Set(replace(last(p), "_keepsing" => "") for p in pairs) == Set(keys(specs))
    for (ds, sp) in pairs
        keep = endswith(sp, "_keepsing")
        fes, vc = specs[replace(sp, "_keepsing" => "")]
        rows = filter(x -> x.dataset == ds && x.spec == sp, ref)
        @test length(rows) == 4 * (H + 1)
        r = estimate(ds, "leads(y)", fes; singletons = keep ? :keep : :drop)
        cv = vcov(vc, r)
        for x in rows
            k = parse(Int, x.h) + 1
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            se = _num(x.se)
            if isnan(se)    # fixest without repair: the square root of a negative variance
                @test cv.V[k][j, j] < 0
            else
                @test sqrt(cv.V[k][j, j]) ≈ se rtol = 1e-10
            end
            @test nobs(r, k - 1) == Int(_num(x.nobs))
        end
        if occursin("_id4", sp)
            @test any(i -> i.psdfixed, vcov(Clustered(:unit, :time; ssc = :none, psd = :raw), r).info)
        end
    end
end

@testitem "R: Mei–Sheng–Shi LP_panel (FE)" tags = [:crosscheck, :r, :vcov] setup = [CrossCheck] begin
    using .CrossCheck: readref, estimate, FIX, _num
    ref = readref(joinpath(FIX, "r", "p1p2_mss_lppanel.csv"))
    AP, OS = AlternatingProjections(), OneSweep()
    specs = Dict(
        "mss_fe_unit" => (["id"], AP, Clustered(:unit; ssc = :none)),
        "mss_fe_unit_robust" => (["id"], AP, Clustered(:unit; ssc = :mss)),
        "mss_te_unit" => (["id", "t"], OS, Clustered(:unit; ssc = :none)),
        "mss_te_unit_robust" => (["id", "t"], OS, Clustered(:unit; ssc = :mss)),
        "mss_fe_tw" => (["id"], AP, Clustered(:unit, :time; ssc = :none, psd = false)),
        "mss_fe_tw_robust_eig" => (["id"], AP, Clustered(:unit, :time; ssc = :mss, psd = :raw)),
        "mss_te_tw" => (["id", "t"], OS, Clustered(:unit, :time; ssc = :none, psd = false)))
    for (spec, (fes, ab, vc)) in specs, ds in ("fx_balanced", "fx_unbalanced")
        r = estimate(ds, "leads(y)", fes; singletons = :keep, absorb = ab)
        cv = vcov(vc, r)
        rows = filter(x -> x.dataset == ds && x.spec == spec, ref)
        @test length(rows) == 4
        for x in rows
            k = parse(Int, x.h) + 1
            @test coef(r, k - 1)[1] ≈ _num(x.b) rtol = 1e-12
            @test sqrt(cv.V[k][1, 1]) ≈ _num(x.se) rtol = 1e-12
        end
    end
end
