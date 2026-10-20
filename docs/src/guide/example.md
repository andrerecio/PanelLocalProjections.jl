# Getting started

A fixed-effects panel local projection on `fx_balanced`, one of the synthetic panels in
`test/fixtures/synthetic/` (30 units, 60 periods).

```@example start
using PanelLocalProjections, DelimitedFiles
path = joinpath(pkgdir(PanelLocalProjections), "test", "fixtures", "synthetic", "fx_balanced.csv")
raw, hdr = readdlm(path, ','; header = true)
data = (; (Symbol(h) => (h in ("id", "t", "grp", "reg") ? Int.(raw[:, j]) : Float64.(raw[:, j]))
           for (j, h) in enumerate(vec(hdr)))...)

r = panellp(@formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + w + fe(id) + fe(t)), data;
            id = :id, time = :t, horizon = 8)
```

Standard errors clustered by unit, with Stata's small-sample correction and t(G − 1)
critical values:

```@example start
s = summarize(r, Clustered(:unit))
```

The same coefficients with Driscoll–Kraay standard errors:

```@example start
summarize(r, DriscollKraaySE(); level = 0.9)
```

Plot with Makie (load a backend such as CairoMakie), or with Plots.jl (`plot(s)`):

```@example start
using CairoMakie
irfplot(s)
```
