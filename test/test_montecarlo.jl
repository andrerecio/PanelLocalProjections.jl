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
