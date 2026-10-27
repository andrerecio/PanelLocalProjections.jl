# The figure of the README and of docs/src/examples/us_states.md: the panel LP-IV and the
# reduced form of the US state example, with t-LAHR and state-clustered 90% bands.
#
#     julia --project examples/us_states/build.jl          # once: downloads and builds the panel
#     julia --project=docs examples/us_states/figure.jl    # docs env: CairoMakie
#
# Writes docs/src/assets/us_states_irf.png.

using PanelLocalProjections, CairoMakie   # @formula is re-exported
include(joinpath(@__DIR__, "csv.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
d = readtable(joinpath(REPO, "docs", "data", "us_states", "processed", "state_panel.csv"))
kw = (; id = :fips, time = :t, horizon = 12, estimator = LagAugmented())

iv = panellp(@formula(ldiff(lpi) ~ (s & firstdiff(gs1) ~ s & mp) + fe(fips) + fe(t)), d; kw...)
rf = panellp(@formula(ldiff(lpi) ~ s & mp + fe(fips) + fe(t)), d; kw..., shock = "s & mp")

fig = Figure(; size = (1000, 400), fontsize = 14)
panels = ((1, iv, TLAHR(refine = false), "Panel LP-IV: s × Δ1-year rate, instrumented",
        "percent, per pp of the rate and sd of exposure", "t-LAHR (A&S A.4), 90%"),
    (2, rf, TLAHR(), "Reduced form: s × monetary policy surprise",
        "percent, per pp surprise and sd of exposure", "t-LAHR (HC2 + IK df), 90%"))
for (c, r, spec, title, ylabel, lab) in panels
    ax = Axis(fig[1, c]; title, xlabel = "quarters after the shock", ylabel)
    irfplot!(ax, summarize(r, spec; level = 0.9); color = :steelblue, bandalpha = 0.3)
    u = summarize(r, Clustered(:unit); level = 0.9)
    for b in (u.lower, u.upper)
        lines!(ax, u.horizon, b; color = :black, linestyle = :dash, linewidth = 1)
    end
    axislegend(ax, [PolyElement(; color = (:steelblue, 0.3)),
            LineElement(; color = :black, linestyle = :dash)],
        [lab, "clustered by state, 90%"]; position = :lb, framevisible = false)
end
out = joinpath(REPO, "docs", "src", "assets", "us_states_irf.png")
mkpath(dirname(out))
save(out, fig; px_per_unit = 2)
println("wrote ", relpath(out, REPO))
