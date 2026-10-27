# Panel data

`panellp` accepts any [Tables.jl](https://github.com/JuliaData/Tables.jl) table in long
format: one row per unit and period. You name the unit and period columns with `id` and
`time`.

## Periods

- **Integer periods** (years, Stata `%tq`/`%tm` codes, any consecutive numbering) are used
  as they are: periods `t` and `t + 1` are adjacent.
- **`Date`/`DateTime` periods** need `frequency`, one of `:year`, `:halfyear`, `:quarter`
  or `:month`. Each date is mapped to its period, so quarterly data dated at the start,
  middle or end of the quarter all work, even when units date the same quarter
  differently. Time effects `fe(date)` and clusters on the date column group by period as
  well. Two rows of the same unit in one period are an error.
- **Float periods** must be integral. Strings (such as `"2001Q1"`) must be converted to
  dates first.

A calendar code such as `YYYYMM` is *not* a consecutive numbering: most of its range is
never observed, so every lead or lag across a year would be undefined. `panellp` warns
when most of the span is unobserved. Convert such codes to dates and pass `frequency`.

## Gaps and missing values

Leads, lags and differences are taken on the (unit, period) grid. If period `t + h` of a
unit is absent, the lead at `t` is undefined and that row leaves the horizon-`h` sample.
Nothing is ever shifted by row position, and no rows are dropped before shifting.
`missing` (and `NaN`) values in the formula variables are treated the same way. A
function of a missing value is missing too: `(w > 0)` is undefined where `w` is
missing, not `false`.

Each horizon's sample is the set of rows where the response and every regressor are
defined, and where every fixed-effect and declared cluster variable is non-missing
(`missing`, `nothing` or `NaN`). With `sample = :common` the intersection over all
horizons is used instead; it is an error if no row survives in all of them.

A horizon without observations has `NaN` coefficients and is reported in a warning. The
same holds for a horizon whose sample leaves the shock absorbed by the fixed effects or
collinear (status `:singular`). It is an error when no horizon can be estimated.

## Restricting the sample

`subset` (a Bool column name or vector) restricts the estimation sample *after* leads and
lags are computed. A lag at the first period of the window therefore still reads the
period before it. Dropping those rows from the data before calling `panellp` would instead
lose that observation:

```julia
r = panellp(f, data; id = :firm, time = :quarter, horizon = 8, subset = :post1990)
```

## Singletons

With fixed effects, a row that is alone in its group of some fixed-effect dimension is
dropped. The removal is repeated until none is left, as `reghdfe` and `xtlp` do.
`singletons = :keep` keeps them, as `xtreg`, `xtscc`, fixest with `fixef.rm = "none"`, and
the Mei–Sheng–Shi replication code do. Under exact absorption (the default) `WithinFE`
point estimates are unaffected. `OneSweep()` estimates are affected, because one sweep is
not the exact within transformation once the sample is unbalanced. So are
`SplitPanelJackknife` estimates with fixed effects other than the unit: a singleton in
another dimension still counts in its unit's rows and so moves the halves. The number of
observations and clusters, and therefore the small-sample factors, change in every case.
