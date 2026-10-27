# Panel LP-IV on the US state panel (docs/src/examples/us_states.md), in the setting of
# Almuzara & Sancibrián (2026, §3.4): the response of state personal income to the 1-year
# Treasury rate, by exposure to manufacturing, instrumented by monetary policy surprises.
#
#     julia --project examples/us_states/build.jl      # once: downloads and builds the panel
#     julia --project examples/us_states/estimate.jl
#
# Writes docs/data/us_states/processed/results.csv (gitignored, like the data).

using PanelLocalProjections, StatsModels, Printf
include(joinpath(@__DIR__, "csv.jl"))

const DATA = normpath(joinpath(@__DIR__, "..", "..", "docs", "data", "us_states", "processed"))
d = readtable(joinpath(DATA, "state_panel.csv"))
H = 12
kw = (; id = :fips, time = :t, horizon = 0:H, estimator = LagAugmented())

# Panel LP-IV: y_{i,t+h} − y_{i,t−1} on s_i ΔGS1_t, instrumented by s_i MP_t, with state and
# quarter effects; lag augmentation adds lags of Δy and of s_i ΔGS1_t (instrumented by lags
# of s_i MP_t), with p̂ chosen by AIC on the reduced form.
iv = panellp(@formula(ldiff(lpi) ~ (s & firstdiff(gs1) ~ s & mp) + fe(fips) + fe(t)), d; kw...)

# The reduced form: the LP on the observed shock s_i MP_t (A&S's main setting).
rf = panellp(@formula(ldiff(lpi) ~ s & mp + fe(fips) + fe(t)), d; kw..., shock = "s & mp")

display(iv)
println()
# First-stage strength: by period (t-LAHR, A&S supplement A.4) and robust (heteroskedastic
# across units and periods, which treats the N·T observations as independent)
fs_t = firststage(iv)
fs_r = firststage(iv; vcov = Robust())
display(fs_t)
println()
display(fs_r)
println()

specs = (TLAHR(refine = false), Clustered(:unit), Clustered(:unit, :time), DriscollKraaySE())
specs_rf = (TLAHR(), Clustered(:unit), Clustered(:unit, :time), DriscollKraaySE())
open(joinpath(DATA, "results.csv"), "w") do io
    println(io, "model,spec,h,coef,se,lower,upper,df,nobs")
    for (model, r, sp) in (("iv", iv, specs), ("reduced_form", rf, specs_rf)), spec in sp
        s = summarize(r, spec; level = 0.9)
        for k in eachindex(s.horizon)
            println(io, join((model, "\"$(s.spec)\"", s.horizon[k], s.coef[k], s.se[k],
                    s.lower[k], s.upper[k], s.df[k], s.nobs[k]), ","))
        end
        display(s)
        println()
    end
end
@printf("lags: p̂ = %d (IV), %d (reduced form)\n", iv.lags.p, rf.lags.p)

# The other numbers quoted on the docs page.
fs_u = firststage(iv; vcov = Clustered(:unit))
display(fs_u)
println()
# the aggregate first stage: ΔGS1_t on MP_t (one state's rows carry the aggregates)
one = findall(==(first(d.fips)), d.fips)
dg, mpq = diff(d.gs1[one]), d.mp[one][2:end]
X = hcat(ones(length(dg)), mpq)
e = dg - X * (X \ dg)
@printf("aggregate first stage: R² of ΔGS1 on MP = %.3f\n",
    1 - sum(abs2, e) / sum(abs2, dg .- sum(dg) / length(dg)))
ratio = [sqrt(vcov(TLAHR(refine = false), iv).V[k][1, 1] / vcov(Clustered(:unit), iv).V[k][1, 1])
         for k in 1:(H + 1)]
@printf("t-LAHR / state-clustered SEs of the IV: %.2f to %.2f (max at h = %d)\n",
    minimum(ratio), maximum(ratio), argmax(ratio) - 1)
mp = d.mp[one]
big = sortperm(abs.(mp); rev = true)[1:6]
@printf("six largest surprise quarters: %s, %.0f%% of the sum of squares\n",
    join(("$(d.year[one][k])Q$(d.quarter[one][k])" for k in sort(big)), ", "),
    100 * sum(abs2, mp[big]) / sum(abs2, mp))
# Imbens–Kolesár df of the reduced form's shock, and the share of its squared influence
# weights a_t² (A&S supplement B.2) in those quarters (package internals: the synthetic
# regression of `TLAHR()`)
cv = vcov(TLAHR(), rf)
for h in (0, 4, 8, 12)
    fit = rf.fits[h + 1]
    kidx = findall(fit.keep)
    sy = fit.synthetic[1]
    des = PanelLocalProjections._ikdesign(sqrt.(sy.weight) .* sy.X[:, kidx])
    j = findfirst(==("s & mp"), rf.coefnames)
    c = findfirst(==(j), kidx)
    a = des.U * (des.V[c, :] ./ des.Sr)
    periods = findall(fit.timepresent)
    inbig = [p in big for p in periods]
    @printf("h = %2d: IK df %.2f, %d periods, share of a_t² in the six quarters %.0f%%\n",
        h, cv.df[h + 1][j], length(periods), 100 * sum(abs2, a[inbig]) / sum(abs2, a))
end
