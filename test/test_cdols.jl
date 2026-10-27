# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testmodule URDGP begin
    using StableRNGs
    # One draw of the Ugarte-Ruiz (2026) design, eqs. (1)–(2):
    #     y_it = ρy y_i,t−1 + β s_it + u_i + e_it,   s_it = μ_i + ρs (s_i,t−1 − μ_i) + v_it,
    # u_i = α(1 − ρy)(1 + s̄_i) (the paper's correlated effects (1 + s̄_i − s̄) without the
    # grand mean s̄ ≈ 0, which only shifts every level), e ~ N(0, 1), v ~ N(0, σv²), a
    # burn-in from the stationary means. The paper sets σv by a signal-to-noise ratio of 20
    # (`σξ`); σv = 1 is not its design. μ_i ~ N(0, μsd²) are unit-specific shock means
    # (zero in the paper's design).
    function draw(rng; N, T, ρy, ρs, β = 1.0, α = 10.0, burn = 200, μsd = 0.0, σv = 1.0)
        id, tt, ys, ss = Int[], Int[], Float64[], Float64[]
        for i in 1:N
            μ = μsd * randn(rng)
            S = zeros(burn + T)
            dev = 0.0
            for k in 1:(burn + T)
                dev = ρs * dev + σv * randn(rng)
                S[k] = μ + dev
            end
            u = α * (1 - ρy) * (1 + sum(view(S, (burn + 1):(burn + T))) / T)
            yp = (u + β * μ) / (1 - ρy)
            for k in 1:(burn + T)
                yk = ρy * yp + β * S[k] + u + randn(rng)
                if k > burn
                    push!(id, i); push!(tt, k - burn); push!(ys, yk); push!(ss, S[k])
                end
                yp = yk
            end
        end
        return (id = id, t = tt, y = ys, s = ss)
    end
    # true response, and the probability limits of the coefficient on s_t−1 and of the
    # response when s_t−1 is omitted
    βh(h, ρy, ρs; β = 1.0) = β * sum(ρy^(h - j) * ρs^j for j in 0:h)
    γh(h, ρy, ρs; β = 1.0) = (ρy^(h + 1) - 1) * β / (1 - ρy * ρs)
    βnolag(h, ρy, ρs; β = 1.0) = βh(h, ρy, ρs; β) - (1 - ρy^(h + 1)) * β * ρs / (1 - ρy * ρs)
    # the level LP without s_t−1 (the paper's FE and SPJ comparators, eq. 24), as T → ∞
    βlevel(h, ρy, ρs; β = 1.0) = βh(h, ρy, ρs; β) + β * ρy^(h + 1) * ρs / (1 - ρy * ρs)
    # xtarsim's σ_ξ for a signal-to-noise ratio (xtarsim.ado, lines 41–49 and 57), with
    # γ = ρy and ρ = ρs: the shock innovation sd of the paper's designs (SNR 20)
    function σξ(γ, ρ, snr; β = 1.0)
        Vφ = 1 / (1 + (γ + ρ)^2 * (γ * ρ - 1) / (1 + γ * ρ) - γ^2 * ρ^2)
        return sqrt((snr - γ^2 / (1 - γ^2)) / (β^2 * Vφ))
    end
end

@testitem "CDOLS specification rules" tags = [:core, :cdols] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using StatsModels
    const PLP = PanelLocalProjections
    d = fixture(:fx_cdols)
    est(f; kw...) = panellp(f, d; id = :id, time = :t, horizon = 3, estimator = CDOLS(), kw...)
    f = @formula(ldiff(y) ~ s + lags(s, 1))
    r = est(f)
    # the estimator is pooled OLS with an intercept on the long difference
    @test coefnames(r) == ["(Intercept)", "s", "L1.s"]
    @test coef(r, 2) == coef(panellp(f, d; id = :id, time = :t, horizon = 3), 2)
    @test est(@formula(ldiff(y) ~ s + lags(s, 1:2) + lags(firstdiff(y), 1:2) + fe(t))) isa
          PanelLocalProjection
    # the rules
    bad = (@formula(leads(y) ~ s + lags(s, 1)),                 # not a long difference
           @formula(ldiff(y) ~ s),                              # no s_t−1
           @formula(ldiff(y) ~ s + lags(s, 2:3)),               # no s_t−1
           @formula(ldiff(y) ~ s + lags(s, 1) + lags(y, 1)),    # a level of y
           @formula(ldiff(y) ~ s + lags(s, 1) + firstdiff(y)),  # Δy_t is not lagged
           @formula(ldiff(y) ~ s + lags(s, 1) + fe(id)),        # unit fixed effects
           # the same long difference, spelled as a cumulated difference
           @formula(cumul(firstdiff(y)) ~ s + lags(s, 1) + lags(y, 1)),
           # both columns of lags(s, 0:1) are shocks by default: L1.s is dated t − 1
           @formula(ldiff(y) ~ lags(s, 0:1)))
    for g in bad
        @test_throws ArgumentError est(g)
    end
    # a lagged response, y_t+h−1 − y_t−2, contains Δy_t−1
    @test_throws ArgumentError est(
        @formula(ldiff(lags(y, 1)) ~ lags(s, 1:2) + lags(firstdiff(y), 1)); shock = "L1.s")
    # a shock dated t − 1 with the base t − 1 would estimate β_{h+1} − β_0
    @test_throws ArgumentError est(@formula(ldiff(y) ~ lags(s, 1:2)); shock = "L1.s")
    # with a lagged response the shock is dated t − 1, and its lag must be there
    @test_throws ArgumentError est(@formula(ldiff(lags(y, 1)) ~ lags(s, 1) + lags(s, 3:3));
        shock = "L1.s")
    # the shock's first lag is found however it is written
    rs = est(@formula(ldiff(y) ~ s + lags(s, 1)))
    @test coef(est(@formula(ldiff(y) ~ lags(s, 0:1)); shock = "L0.s"), 2) ≈ coef(rs, 2)
    @test est(@formula(ldiff(y) ~ firstdiff(s) + firstdiff(lags(s, 1))); shock = "D.s") isa
          PanelLocalProjection
    @test est(@formula(ldiff(lags(y, 1)) ~ lags(s, 1:2) + lags(firstdiff(y), 2:2));
        shock = "L1.s") isa PanelLocalProjection
    # cumul(firstdiff(y)) is ldiff(y)
    rd = est(@formula(ldiff(y) ~ s + lags(s, 1) + lags(firstdiff(y), 1)))
    @test all(coef(est(@formula(cumul(firstdiff(y)) ~ s + lags(s, 1) +
                                lags(firstdiff(y), 1))), h) ≈ coef(rd, h) for h in 0:3)
    # the shock list
    @test_throws ArgumentError est(f; shock = String[])
    @test_throws ArgumentError est(f; shock = ["s", "s"])
    # fixed effects: nested in the units in a horizon's sample (here through `subset`),
    # and groups with a single unit
    dp = merge(d, (; pair = cld.(d.id, 2), almost = [i <= 2 ? 1 : i for i in d.id]))
    estp(g; kw...) = panellp(g, dp; id = :id, time = :t, horizon = 2, estimator = CDOLS(), kw...)
    @test estp(@formula(ldiff(y) ~ s + lags(s, 1) + fe(pair))) isa PanelLocalProjection
    @test_throws ArgumentError estp(@formula(ldiff(y) ~ s + lags(s, 1) + fe(pair));
        subset = isodd.(dp.id))
    @test_logs (:warn, r"single unit") match_mode = :any estp(
        @formula(ldiff(y) ~ s + lags(s, 1) + fe(almost)))
    # unit fixed effects in disguise: an interaction nested in the units
    dg = merge(d, (; half = Int.(d.t .> 25)))
    @test_throws ArgumentError panellp(@formula(ldiff(y) ~ s + lags(s, 1) + fe(id) & fe(half)),
        dg; id = :id, time = :t, horizon = 1, estimator = CDOLS())
    # cumul(g) of a growth rate g is the long difference of its level
    p = PLP.PanelIndex(d.id, d.t)
    dd = merge(d, (; g = Float64.(d.y) .- PLP._lag(p, Float64.(d.y), 1)))
    rc = panellp(@formula(cumul(g) ~ s + lags(s, 1) + lags(g, 1)), dd; id = :id, time = :t,
        horizon = 3, estimator = CDOLS())
    rl = est(@formula(ldiff(y) ~ s + lags(s, 1) + lags(firstdiff(y), 1)))
    @test all(coef(rc, h) ≈ coef(rl, h) for h in 0:3)
    @test_throws ArgumentError panellp(@formula(cumul(g) ~ s + lags(s, 1) + g), dd; id = :id,
        time = :t, horizon = 1, estimator = CDOLS())
    # a pre-built outcome column needs validate = false
    db = merge(dd, (; ld0 = dd.g))
    @test_throws ArgumentError panellp(@formula(ld0 ~ s + lags(s, 1)), db; id = :id,
        time = :t, horizon = 0, estimator = CDOLS())
    @test coef(panellp(@formula(ld0 ~ s + lags(s, 1)), db; id = :id, time = :t, horizon = 0,
        estimator = CDOLS(; validate = false)), 0) ≈ coef(r, 0)
end

@testitem "xtarsim panels have the documented moments" tags = [:core, :cdols, :fixtures] setup = [PLPSupport, URDGP] begin
    using Statistics
    const PLP = PanelLocalProjections
    using .PLPSupport: PLPDGPs
    using .URDGP: σξ
    for (d, ρs, ρy) in (("s8y0", 0.8, 0.0), ("s0y8", 0.0, 0.8), ("s8y8", 0.8, 0.8))
        x = PLPDGPs.read_csv(joinpath(dirname(@__DIR__), "test", "fixtures", "stata",
            "xtarsim_$d.csv"))
        @test length(x.id) == 5000
        p = PLP.PanelIndex(x.id, x.t)
        ls, ly = PLP._lag(p, Float64.(x.s), 1), PLP._lag(p, Float64.(x.y), 1)
        ok = isfinite.(ls)
        ρ̂ = ls[ok] \ x.s[ok]
        ξ = x.s[ok] .- ρ̂ .* ls[ok]
        @test abs(ρ̂ - ρs) < 4 * sqrt((1 - ρs^2) / count(ok)) + 0.01
        @test std(ξ) ≈ σξ(ρy, ρs, 20) rtol = 0.03
        # the outcome equation with σ_ε = 1 and the paper's correlated effects
        ε = (x.y .- ρy .* ly .- x.s .- x.eta)[ok]
        @test abs(mean(ε)) < 0.05
        @test std(ε) ≈ 1 rtol = 0.05
        s̄ = [mean(x.s[x.id .== i]) for i in 1:100]
        @test x.eta ≈ [10 * (1 - ρy) * (1 + s̄[i] - mean(s̄)) for i in x.id] rtol = 1e-12
    end
end

@testitem "the native design reproduces the xtarsim panels" tags = [:cdols, :mc] setup = [URDGP, PLPSupport] begin
    using .PLPSupport: PLPDGPs
    using StableRNGs, Statistics, StatsModels
    const PLP = PanelLocalProjections
    # URDGP with the paper's shock scale (SNR 20) against each committed xtarsim panel:
    # the sd and AR(1) coefficient of s, the sd of Δy and CD-OLS at h = 0 and 4, within
    # 4 sd of 40 native draws
    function moments(d)
        p = PLP.PanelIndex(d.id, d.t)
        ls = PLP._lag(p, Float64.(d.s), 1)
        dy = Float64.(d.y) .- PLP._lag(p, Float64.(d.y), 1)
        ok = isfinite.(ls)
        r = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), d; id = :id, time = :t, horizon = 4,
            estimator = CDOLS())
        return [std(d.s), ls[ok] \ d.s[ok], std(dy[isfinite.(dy)]), coef(r, 0)[2], coef(r, 4)[2]]
    end
    R = 40
    for (name, ρs, ρy) in (("s8y0", 0.8, 0.0), ("s0y8", 0.0, 0.8), ("s8y8", 0.8, 0.8))
        x = PLPDGPs.read_csv(joinpath(dirname(@__DIR__), "test", "fixtures", "stata",
            "xtarsim_$name.csv"))
        rng = StableRNG(14)
        sims = reduce(hcat, [moments(URDGP.draw(rng; N = 100, T = 50, ρy, ρs,
                                σv = URDGP.σξ(ρy, ρs, 20))) for _ in 1:R])
        μ, σ = vec(mean(sims; dims = 2)), vec(std(sims; dims = 2))
        @test all(abs.(moments(x) .- μ) .<= 4 .* σ .* sqrt(1 + 1 / R))
    end
end

@testitem "CD-OLS probability limits" tags = [:cdols, :mc] setup = [URDGP] begin
    using StableRNGs, StatsModels
    for (ρy, ρs) in ((0.0, 0.8), (0.8, 0.0), (0.8, 0.8))
        d = URDGP.draw(StableRNG(11); N = 4000, T = 40, ρy, ρs)
        r = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), d; id = :id, time = :t,
            horizon = 9, estimator = CDOLS())
        r0 = panellp(@formula(ldiff(y) ~ s), d; id = :id, time = :t, horizon = 9,
            estimator = CDOLS(; validate = false))
        V, V0 = vcov(Clustered(:unit), r), vcov(Clustered(:unit), r0)
        for h in (0, 4, 9)
            k = h + 1
            @test abs(coef(r, h)[2] - URDGP.βh(h, ρy, ρs)) < 4 * sqrt(V.V[k][2, 2])
            @test abs(coef(r, h)[3] - URDGP.γh(h, ρy, ρs)) < 4 * sqrt(V.V[k][3, 3])
            @test abs(coef(r0, h)[2] - URDGP.βnolag(h, ρy, ρs)) < 4 * sqrt(V0.V[k][2, 2])
        end
    end
    # The estimands. The fixed-effects LP in levels without s_t−1 (the paper's FE and SPJ
    # comparators) converges to β_h + βρ_y^{h+1}ρ_s/(1 − ρ_yρ_s), +178% of β at h = 0 under
    # joint persistence (the paper's SPJ "bias" of 177.45% in its Table 9). With s_t−1 it
    # targets β_h, as CD-OLS. A long panel makes the Nickell bias small.
    d = URDGP.draw(StableRNG(13); N = 100, T = 2000, ρy = 0.8, ρs = 0.8)
    kw = (; id = :id, time = :t, horizon = 4)
    lv = panellp(@formula(leads(y) ~ s + fe(id)), d; kw...)
    ll = panellp(@formula(leads(y) ~ s + lags(s, 1) + fe(id)), d; kw..., shock = "s")
    for h in (0, 4)
        @test coef(lv, h)[1] ≈ URDGP.βlevel(h, 0.8, 0.8) rtol = 0.03
        @test coef(ll, h)[findfirst(==("s"), coefnames(ll))] ≈ URDGP.βh(h, 0.8, 0.8) rtol = 0.03
    end
end

@testitem "CD-OLS is inconsistent with unit-specific shock means" tags = [:cdols, :mc] setup = [URDGP] begin
    using StableRNGs, StatsModels
    # the documented failure: shock means that differ across units stay in the pooled
    # residual of s_t on s_t−1
    d = URDGP.draw(StableRNG(12); N = 4000, T = 40, ρy = 0.0, ρs = 0.8, μsd = 2.0)
    r = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), d; id = :id, time = :t, horizon = 9,
        estimator = CDOLS())
    se = stderror(vcov(Clustered(:unit), r); term = "s")
    @test (coef(r, 9)[2] - URDGP.βh(9, 0.0, 0.8)) / se[10] > 10
end

@testitem "CD-OLS standard errors: robust vs clustered" tags = [:cdols, :mc] setup = [URDGP] begin
    using StableRNGs, StatsModels, Statistics
    # The paper's designs (N = 100, T = 50, signal-to-noise ratio 20): the long-difference
    # error is serially correlated, so heteroskedasticity-robust standard errors are not
    # valid, while unit-clustered ones match the sampling sd. The paired ratio
    # mean(se_robust / se_clustered) is ten times more precise than a ratio to the sd.
    # (2000 replications, seed 22: with ρ_s = 0, ρ_y = 0.8, robust/sd 1.19, 1.51, 1.33 at
    # h = 2, 4, 9, paired 1.20, 1.46, 1.33; with ρ_s = ρ_y = 0.8, robust/sd 0.91, 0.88, 0.93
    # at h = 1, 2, 4, paired 0.93, 0.89, 0.91; clustered/sd within 0.98–1.04 in both.)
    R, hs = 200, (0, 1, 2, 4, 9)
    function mc(ρs, ρy)
        rng = StableRNG(22)
        σv = URDGP.σξ(ρy, ρs, 20)
        est, seR, seC = zeros(R, 5), zeros(R, 5), zeros(R, 5)
        for rep in 1:R
            d = URDGP.draw(rng; N = 100, T = 50, ρy, ρs, burn = 100, σv)
            r = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), d; id = :id, time = :t,
                horizon = 9, estimator = CDOLS())
            sr = stderror(vcov(Robust(), r); term = "s")
            sc = stderror(vcov(Clustered(:unit), r); term = "s")
            for (k, h) in enumerate(hs)
                est[rep, k], seR[rep, k], seC[rep, k] = coef(r, h)[2], sr[h + 1], sc[h + 1]
            end
        end
        sd = vec(std(est; dims = 1))
        return (; paired = vec(mean(seR ./ seC; dims = 1)), clust = vec(mean(seC; dims = 1)) ./ sd)
    end
    # a ratio of standard deviations has a Monte Carlo error of about 1/√(2R)
    tol = 4 / sqrt(2 * (R - 1))
    # ρ_s = 0, ρ_y = 0.8 (Table 5): robust standard errors overstate at h ≥ 2
    a = mc(0.0, 0.8)
    @test abs(a.paired[1] - 1) < 0.05
    @test a.paired[4] > 1.25 && a.paired[5] > 1.2
    @test all(abs.(a.clust .- 1) .< tol)
    # ρ_s = ρ_y = 0.8 (the paper's main design): they understate at h = 1–4
    b = mc(0.8, 0.8)
    @test all(b.paired[2:4] .< 0.97)
    @test all(abs.(b.clust .- 1) .< tol)
end

@testitem "Ugarte-Ruiz (2026), Fig. 9: CD-OLS on the Cerra–Saxena data" tags = [:local_mss, :cdols] setup = [PLPSupport] begin
    using .PLPSupport: PLPDGPs
    using StatsModels
    const PLP = PanelLocalProjections
    cs = PLPDGPs.read_csv(joinpath(dirname(@__DIR__), "docs", "references", "code",
        "panel-lp-replication", "applications", "empirical_CS_f3.csv"))
    p = PLP.PanelIndex(cs.cntyNUM, cs.obs)
    # the h = 0 outcome is the growth rate g_t, the lead of l1GRRT_WB
    d = merge(cs, (; cf0GRRT_WB = PLP._lead(p, PLP._nanfloat(cs.l1GRRT_WB), 1)))
    rhs = "CRISIS + l1CRISIS + l2CRISIS + l3CRISIS + l4CRISIS + " *
          "l1GRRT_WB + l2GRRT_WB + l3GRRT_WB + l4GRRT_WB"
    # read off the paper's Fig. 9 (CMOLS), to about ±0.2
    fig9 = [-1.1, -2.5, -3.4, -4.1, -4.5, -4.65, -5.2, -4.9, -5.1, -5.3, -5.7]
    for h in 0:10
        f = eval(Meta.parse("@formula(cf$(h)GRRT_WB ~ $rhs)"))
        r = panellp(f, d; id = :cntyNUM, time = :obs, horizon = 0,
            estimator = CDOLS(; validate = false))
        @test abs(coef(r, 0)[2] - fig9[h + 1]) < 0.1
    end
end
