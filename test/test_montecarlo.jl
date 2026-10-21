# Monte Carlo smoke reproductions of published simulation results, at reduced size.
# Tolerances combine our Monte Carlo error with the reference's (4 standard errors), and
# the seeds are fixed once: never retune a seed to make a gate pass.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "MSS toy design: FE Nickell bias and SPJ (Figures 1–3)" tags = [:mc, :spj] setup = [PLPSupport] begin
    using .PLPSupport: PLPDGPs
    using StableRNGs, Statistics, StatsFuns
    # Mei, Sheng & Shi (2026) §2.4 / simulations/toy_-0.6.xlsx: ρ = 0.8, (N, T) = (30, 60),
    # β₀ = −0.6, 1000 replications; reference mean IRFs and 95% coverage of FE / SPJ.
    ref = Dict(1 => (-0.4603, -0.4825, 0.796, 0.942),
               5 => (-0.1277, -0.1950, 0.374, 0.912),
               10 => (0.0271, -0.0371, 0.211, 0.856))
    R, Rref = 300, 1000
    rng = StableRNG(20261030)
    f = @formula(leads(y) ~ x + fe(id))
    est = Dict(:FE => zeros(R, 11), :SPJ => zeros(R, 11))
    ses = Dict(:FE => zeros(R, 11), :SPJ => zeros(R, 11))
    cover = Dict(:FE => falses(R, 11), :SPJ => falses(R, 11))
    z = norminvcdf(0.975)
    truth = [-0.6 * 0.8^h for h in 0:10]
    for rep in 1:R
        d = PLPDGPs.mss_toy(rng)
        for (lab, e) in ((:FE, WithinFE()), (:SPJ, SplitPanelJackknife()))
            r = panellp(f, d; id = :id, time = :t, horizon = 10, estimator = e)
            # their simulation code always applies N/(N−1)·(n−1)/(n−K) (LP_panel_fe.R:259)
            se = stderror(vcov(Clustered(:unit; ssc = :mss), r); term = "x")
            b = coefpath(r)
            est[lab][rep, :] = b
            ses[lab][rep, :] = se
            cover[lab][rep, :] = abs.(b .- truth) .<= z .* se
        end
    end
    for (h, (mfe, mspj, cfe, cspj)) in ref
        for (lab, m, c) in ((:FE, mfe, cfe), (:SPJ, mspj, cspj))
            b = est[lab][:, h + 1]
            @test abs(mean(b) - m) <= 4 * std(b) * sqrt(1 / R + 1 / Rref)
            p = mean(cover[lab][:, h + 1])
            @test abs(p - c) <= 4 * sqrt(c * (1 - c) * (1 / R + 1 / Rref))
        end
    end
    # the FE bias at h = 5 is large relative to its Monte Carlo error; SPJ removes most of it
    @test mean(est[:FE][:, 6]) - truth[6] > 0.05
    @test abs(mean(est[:SPJ][:, 6]) - truth[6]) < 0.01
    # the SPJ variance itself: at h = 1 and 5, where the SPJ bias is negligible, the mean
    # standard error matches the sampling sd (a sd ratio has relative error ≈ 1/√(2R)).
    # Coverage alone would not detect a variance 15–20% too small, such as one built from
    # FE-type scores X̃·ê_spj (mean se / sd ≈ 0.86 at h = 1, 0.78 at h = 5)
    for h in (1, 5)
        @test abs(mean(ses[:SPJ][:, h + 1]) / std(est[:SPJ][:, h + 1]) - 1) <= 3 / sqrt(2R)
    end
end

@testmodule Table1 begin
    # Almuzara & Sancibrián (2026) Table 1: coverage of 90% intervals at h = 0, 1, 2, 5,
    # T = 30, N = 1500, ρ = 0.7, 5000 replications.
    const REF = Dict(
        0.9 => (oneway = (0.494, 0.372, 0.344, 0.376), twoway = (0.819, 0.726, 0.725, 0.811),
            dk = (0.716, 0.725, 0.695, 0.782), plain = (0.851, 0.842, 0.811, 0.806),
            recommended = (0.903, 0.899, 0.872, 0.875)),
        0.1 => (oneway = (0.888, 0.878, 0.866, 0.880), twoway = (0.841, 0.829, 0.814, 0.820),
            dk = (0.809, 0.791, 0.789, 0.810), plain = (0.848, 0.849, 0.836, 0.835),
            recommended = (0.903, 0.902, 0.891, 0.895)))
    const RREF = 5000
    coverage(out, k) = vec(sum(reduce(hcat, getfield.(out, k)); dims = 2)) ./ length(out)
end

@testitem "A&S stylized design: t-LAHR coverage (Table 1)" tags = [:mc, :tlahr] setup = [TLAHRSupport, Table1] begin
    using .TLAHRSupport: stylized, table1
    using .Table1: REF, RREF, coverage
    using Logging, StableRNGs, Statistics
    # Reduced to N = 200 and R = 300: the calibration σ_u² = N(1 − R̄²)/R̄² (their fn 10)
    # keeps the macro share of the variation as N changes. R̄² = 0.9.
    R = 300
    rng = StableRNG(20261031)
    out = with_logger(NullLogger()) do
        [table1(stylized(rng)) for _ in 1:R]
    end
    ref = REF[0.9]
    tol(c) = 4 * sqrt(c * (1 - c) * (1 / R + 1 / RREF))
    for k in (:plain, :recommended)
        for (p, c) in zip(coverage(out, k), getfield(ref, k))
            @test abs(p - c) <= tol(c)
        end
    end
    # unit clustering fails with a strong macro signal, the recommendation does not
    @test all(coverage(out, :oneway) .<= 0.6)
    @test all(coverage(out, :recommended) .- coverage(out, :oneway) .>= 0.3)
    @test all(coverage(out, :recommended) .>= coverage(out, :twoway) .- 0.03)
    # the refinement on the same draws: a paired comparison, whose Monte Carlo error is far
    # below that of a cell (the published gain is .052–.069; without the Imbens–Kolesár df
    # it falls to about .03)
    gain = coverage(out, :recommended) .- coverage(out, :plain)
    @test all(gain .>= 0) && mean(gain) >= 0.045
    # the AIC rule mostly picks one lag (the design is a VAR(1))
    @test mean(getfield.(out, :phat) .== 1) >= 0.6
end

@testitem "A&S stylized design: Table 1 at the paper's size" tags = [:mc_full, :tlahr] setup = [TLAHRSupport, Table1] begin
    using .TLAHRSupport: stylized, table1
    using .Table1: REF, RREF, coverage
    using Logging, StableRNGs, Statistics
    # N = 1500, R = 2000; both calibrations of σ_u² at R̄² = 0.9 (which one reproduces the
    # unit-cluster row tells which the paper used). Reported with @info.
    R = 2000
    for (R2, cal) in ((0.9, :fn10), (0.9, :exact), (0.1, :fn10))
        rng = StableRNG(20261101)
        out = with_logger(NullLogger()) do
            [table1(stylized(rng; N = 1500, R2, cal)) for _ in 1:R]
        end
        ref = REF[R2]
        rows = Dict(k => coverage(out, k) for k in (:oneway, :twoway, :dk, :plain, :recommended))
        fmt(v) = join((string(round(x; digits = 3)) for x in v), ", ")
        @info "Table 1, R̄² = $R2, calibration $cal" oneway = fmt(rows[:oneway]) twoway = fmt(rows[:twoway]) dk = fmt(rows[:dk]) plain = fmt(rows[:plain]) recommended = fmt(rows[:recommended]) share_p1 = mean(getfield.(out, :phat) .== 1)
        cal === :fn10 || continue
        tol(c) = max(0.03, 3.5 * sqrt(c * (1 - c) * (1 / R + 1 / RREF)))
        for k in (:plain, :recommended), (p, c) in zip(rows[k], getfield(ref, k))
            @test abs(p - c) <= tol(c)
        end
        if R2 == 0.9
            @test all(rows[:recommended] .>= rows[:twoway] .- 0.03)
            @test all(rows[:recommended] .>= rows[:oneway] .- 0.03)
        end
    end
end
