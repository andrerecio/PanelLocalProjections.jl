# Cross-checks of the split-panel jackknife against Stata (xtlp m(spj); xtspj) and the
# Mei–Sheng–Shi R function LP_panel(method = "SPJ"), on the synthetic panels:
# test/fixtures/stata/p2_spj.csv (crosscheck/stata/p2_spj.do) and
# test/fixtures/r/p1p2_mss_lppanel.csv (crosscheck/R/p1p2_mss_lppanel.R).
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "Stata: xtlp m(spj) and xtspj" tags = [:crosscheck, :stata, :spj] setup = [CrossCheck] begin
    using .CrossCheck: readref, estimate, TERMS, FIX, RHS, NOLAG, _num
    ref = readref(joinpath(FIX, "stata", "p2_spj.csv"))
    SPJ = SplitPanelJackknife()
    DKP(m; kw...) = DriscollKraaySE(m; timelags = :position, kw...)
    ID, IDT = ["id"], ["id", "t"]
    # spec => (fixed effects, singletons, variance, tolerance on SEs, compare df, RHS)
    sp(fes, sing, vc, tol = 1e-12; df = true, rhs = RHS) = (; fes, sing, vc, tol, df, rhs)
    specs = Dict(
        "xtlp_spj_fe_cl_id_ase" => sp(ID, :drop, Clustered(:unit; ssc = :none); df = false),
        "xtlp_spj_fe_cl_id" => sp(ID, :drop, Clustered(:unit)),
        "xtlp_spj_fe_cl_t" => sp(ID, :drop, Clustered(:time)),
        "xtlp_spj_fe_cl_reg" => sp(ID, :drop, Clustered(:reg)),
        "xtlp_spj_fe_cl_idt" => sp(ID, :drop, Clustered(:unit, :time)),
        "xtlp_spj_fe_cl_idt_ase" => sp(ID, :drop, Clustered(:unit, :time; ssc = :none); df = false),
        "xtlp_spj_fe_dk" => sp(ID, :drop, DKP(:nw94)),
        "xtlp_spj_fe_dk2_ase" => sp(ID, :drop, DKP(2; ssc = :none); df = false),
        "xtlp_spj_fe_dk2_nodfadj" => sp(ID, :drop, DKP(2; ssc = :xtscc)),
        "xtlp_spj_fe_robust" => sp(ID, :drop, Robust()),
        "xtlp_spj_fe_cl_id_ase_keepsing" => sp(ID, :keep, Clustered(:unit; ssc = :none); df = false),
        "xtlp_spj_fe_nolag_cl_id" => sp(ID, :drop, Clustered(:unit); rhs = NOLAG),
        "xtlp_spj_fe_nolag_cl_id_keepsing" => sp(ID, :keep, Clustered(:unit); rhs = NOLAG),
        # xtlp's two-way demeaning stops on an absolute step of 1e-9
        "xtlp_spj_tfe_cl_id_ase" => sp(IDT, :drop, Clustered(:unit; ssc = :none), 1e-9; df = false),
        "xtlp_spj_tfe_cl_id" => sp(IDT, :drop, Clustered(:unit), 1e-9),
        "xtlp_spj_tfe_cl_t" => sp(IDT, :drop, Clustered(:time), 1e-9),
        "xtlp_spj_tfe_cl_idt_ase" => sp(IDT, :drop, Clustered(:unit, :time; ssc = :none), 1e-9; df = false),
        "xtlp_spj_tfe_robust" => sp(IDT, :drop, Robust(), 1e-9),
        "xtlp_spj_tfe_dk2_ase" => sp(IDT, :drop, DKP(2; ssc = :none), 1e-9; df = false),
        "xtlp_spj_tfe_cl_id_keepsing" => sp(IDT, :keep, Clustered(:unit), 1e-9),
        "xtlp_spj_tfe_nolag_cl_id" => sp(IDT, :drop, Clustered(:unit), 1e-9; rhs = NOLAG),
        "xtlp_spj_tfe_nolag_cl_id_keepsing" => sp(IDT, :keep, Clustered(:unit), 1e-9; rhs = NOLAG))
    # every xtlp spec of the fixture is mapped, with all its rows (terms × h = 0..3)
    @test Set(x.spec for x in ref if startswith(x.spec, "xtlp")) == Set(keys(specs))
    for (spec, s) in specs, ds in ("fx_balanced", "fx_unbalanced")
        r = estimate(ds, "leads(y)", s.fes; singletons = s.sing, estimator = SPJ, rhs = s.rhs)
        cv = vcov(s.vc, r)
        rows = filter(x -> x.dataset == ds && x.spec == spec, ref)
        @test length(rows) == (s.rhs == RHS ? 16 : 8)
        for x in rows
            k = parse(Int, x.h) + 1
            j = findfirst(==(TERMS[x.term]), coefnames(r))
            @test coef(r, k - 1)[j] ≈ _num(x.b) rtol = 1e-10
            @test sqrt(cv.V[k][j, j]) ≈ _num(x.se) rtol = s.tol
            @test nobs(r, k - 1) == Int(_num(x.N))
            s.df && @test cv.df[k][j] == _num(x.df_r)
            G = _num(x.G)
            isnan(G) || s.vc isa DriscollKraaySE || @test cv.info[k].groups == G
        end
    end
    # xtspj (Sun & Dhaene 2019), model(regress) method(parm): point estimates on the
    # balanced panel where every unit has the same even number of periods (h = 1, 3)
    r = estimate("fx_balanced", "leads(y)", ["id"]; estimator = SPJ)
    xs = filter(x -> x.spec == "xtspj_parm", ref)
    @test length(xs) == 8
    for x in xs
        h = parse(Int, x.h)
        @test coef(r, h)[findfirst(==(TERMS[x.term]), coefnames(r))] ≈ _num(x.b) rtol = 1e-10
    end
end

@testitem "R: Mei–Sheng–Shi LP_panel (SPJ)" tags = [:crosscheck, :r, :spj] setup = [CrossCheck] begin
    using .CrossCheck: readref, estimate, FIX, _num
    ref = readref(joinpath(FIX, "r", "p1p2_mss_lppanel.csv"))
    SPJ = SplitPanelJackknife()
    AP, OS = AlternatingProjections(), OneSweep()
    specs = Dict(
        "mss_spj_unit" => (["id"], AP, Clustered(:unit; ssc = :none)),
        "mss_spj_unit_robust" => (["id"], AP, Clustered(:unit; ssc = :mss)),
        "mss_spj_te_unit" => (["id", "t"], OS, Clustered(:unit; ssc = :none)),
        "mss_spj_te_unit_robust" => (["id", "t"], OS, Clustered(:unit; ssc = :mss)),
        "mss_spj_tw" => (["id"], AP, Clustered(:unit, :time; ssc = :none, psd = false)),
        "mss_spj_te_tw" => (["id", "t"], OS, Clustered(:unit, :time; ssc = :none, psd = false)),
        "mss_spj_tw_robust_eig" => (["id"], AP, Clustered(:unit, :time; ssc = :mss, psd = :raw)),
        "mss_spj_te_tw_robust_eig" => (["id", "t"], OS,
            Clustered(:unit, :time; ssc = :mss, psd = :raw)))
    @test Set(x.spec for x in ref if startswith(x.spec, "mss_spj")) == Set(keys(specs))
    # fx_balanced_id4: units 1–4, where the two-way variances are not PSD
    for (spec, (fes, ab, vc)) in specs, ds in ("fx_balanced", "fx_unbalanced", "fx_balanced_id4")
        r = estimate(ds, "leads(y)", fes; singletons = :keep, absorb = ab, estimator = SPJ)
        cv = vcov(vc, r)
        rows = filter(x -> x.dataset == ds && x.spec == spec, ref)
        @test length(rows) == 4
        for x in rows
            k = parse(Int, x.h) + 1
            @test coef(r, k - 1)[1] ≈ _num(x.b) rtol = 1e-12
            if isnan(_num(x.se))    # R's sqrt of a negative variance
                @test cv.V[k][1, 1] < 0
            else
                @test sqrt(cv.V[k][1, 1]) ≈ _num(x.se) rtol = 1e-12
            end
        end
    end
    # the eigenvalue repair (eigens = T) is exercised: it fires on the four-unit panel
    r4 = estimate("fx_balanced_id4", "leads(y)", ["id"]; singletons = :keep, estimator = SPJ)
    @test any(i -> i.psdfixed, vcov(Clustered(:unit, :time; ssc = :mss, psd = :raw), r4).info)
end
