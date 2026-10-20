# Formula

```julia
panellp(@formula(leads(y) ~ x + lags(x, 1:4) + lags(y, 1:4) + fe(id) + fe(t)), data;
        id = :id, time = :t, horizon = 10)
```

## Response (left-hand side)

| term | response at horizon ``h`` | Stata `xtlp` / `locproj` |
|---|---|---|
| `leads(y)` | ``y_{i,t+h}`` | `ytransf(level)` |
| `ldiff(y)` | ``y_{i,t+h} - y_{i,t-1}`` | `ytransf(cmltdiff)`, `tr(cmlt)` |
| `anchor(y, z)` | ``y_{i,t+h} - z_{i,t}``; `anchor(y, y)` is ``y_{i,t+h} - y_{i,t}`` | — |
| `cumul(y)` | ``\sum_{k=0}^{h} y_{i,t+k}`` | `ytransf(cmltsum)` |
| `leads(firstdiff(y))` | ``y_{i,t+h} - y_{i,t+h-1}`` | `ytransf(diff)` |
| `y` | a pre-built outcome; only with `horizon = 0` | `hor(0)` |

Functions of columns work inside the response, e.g. `leads(100 * log(y))`.

## Regressors (right-hand side)

- Columns and functions of columns: `x`, `log(x)`, `x^2`, `(x > 0)`,
  `ifelse(x > 0, x, 0)`. A function of a missing value is missing, so `(x > 0)` is
  undefined where `x` is missing, not `false`. StatsModels' `lag(x, k)` shifts by row
  position and is rejected; use `lags(x, k)`.
- Interactions: `s & X`, typically a unit-level exposure times an aggregate shock.
- Lags: `lags(x, p)` gives ``x_{t-1}, …, x_{t-p}``, named `L1.x`, …; `lags(x, a:b)` with
  ``a ≥ 0`` gives ``x_{t-a}, …, x_{t-b}``. They also apply to interactions,
  `lags(s & X, 1)`, and to first differences, `lags(firstdiff(y), 1:4)` (named `L1.D.y`).
- First differences: `firstdiff(x)`, named `D.x`.
- Fixed effects: `fe(id)`, `fe(t)`, and interactions such as `fe(sector) & fe(quarter)`.
  They are absorbed, not estimated.
- An intercept is included when there are no fixed effects, unless the formula has `0 +`.

Categorical regressors are not supported in v0.1: absorb them with `fe()` or build dummy
columns. Panel LP-IV (`(endo ~ inst)`) is planned for v0.2.

## The shock

`shock` names the coefficient(s) whose impulse response you want. By default it is the
first regressor. StatsModels lists main effects before interactions, so a formula with an
interaction and other regressors requires `shock` explicitly:

```julia
panellp(@formula(ldiff(k) ~ lev & mp + lev + fe(firm) + fe(sector) & fe(quarter)), data;
        id = :firm, time = :quarter, horizon = 12, shock = "lev & mp")
```

A shock that the fixed effects absorb at every horizon is an error, not a zero
coefficient. An example is an aggregate shock with `fe(time)`. A shock absorbed or
collinear only in some horizons' samples makes those horizons `:singular`, with `NaN`
coefficients. Other regressors absorbed by the fixed effects, or collinear with the rest,
are omitted (their coefficient is `NaN`). Each call reports all of these in a warning,
and `show(r)` lists them.

## Formulas built in code

The helpers return StatsModels terms, so `~`, `+` and `&` compose them, for example to
loop over lag lengths:

```julia
using StatsModels: term
for p in 1:4
    f = ldiff(term(:y)) ~ term(:x) + lags(term(:y), p) + fe(term(:id)) + fe(term(:t))
    r = panellp(f, data; id = :id, time = :t, horizon = 8)
end
```

## Helpers from other packages

The helper names (`leads`, `lags`, `fe`, …) are shared with LocalProjections.jl,
Regress.jl and FixedEffectModels.jl, and their helpers are accepted inside `panellp`
formulas. A same-named function from any other module is an error. If two packages export
the same helper, qualify it: `PanelLocalProjections.leads(y)`.
