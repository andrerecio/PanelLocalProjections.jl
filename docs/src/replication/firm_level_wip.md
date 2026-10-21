# Firm-level applications (work in progress)

Almuzara & Sancibrián (2026, §5) re-estimate three firm-level studies of monetary
policy with t-LAHR: Ottonello & Winberry (2020, OW20), Cloyne, Ferreira, Froemel &
Surico (2023, CFFS23) and Jeenas & Lagos (2024, JL24). All three use quarterly Compustat;
OW20 and JL24 also use CRSP, CFFS23 WorldScope. These data are proprietary, so this page
maps the specifications to `panellp` calls for users who have them. **None of them has
been run on the data**: the calls below were run only on synthetic panels of the same
shape, to check that they estimate. The specifications are summarised from the papers
and from Almuzara & Sancibrián's Table 2; check them against the papers before use.

Almuzara & Sancibrián's common equation (eq. 16) is

```math
\Delta k_{i,t+h} = a_{i,h} + d_{h,g(i)t} + \beta(h)' s_{i,t-1} X_t + \eta(h)' W_{it} + v_{h,it},
```

with firm effects, group × quarter effects, a lagged firm characteristic ``s_{i,t-1}``
interacted with the monetary shock ``X_t``, and controls. Their re-estimation adds
``\hat p(h)`` lags of the outcome and of ``s_{i,t-1}X_t`` (`LagAugmented()`) and reports
t-LAHR intervals (`TLAHR()`); quarters can be integers or `Date`s with
`frequency = :quarter`.

## Ottonello & Winberry (2020)

- Outcome: ``\log K_{i,t+h} - \log K_{i,t-1}``, the long difference of the capital
  stock (`ldiff(logk)`).
- Exposure: leverage or distance to default at ``t - 1``, standardised.
- Shock: high-frequency monetary-policy surprises aggregated to quarters
  (expansionary positive).
- Effects: firm, and sector (8 SIC groups) × quarter.
- Controls: lagged firm controls (sales growth, size, current assets) and the exposure
  interacted with lagged GDP growth.
- Their standard errors: clustered by firm and quarter, 1990Q1–2007Q4.

```julia
f = @formula(ldiff(logk) ~ lags(dd, 1) & mp + lags(dd, 1) + lags(salesgr, 1) +
             lags(size, 1) + lags(curassets, 1) + lags(dd, 1) & lags(gdpgr, 1) +
             fe(firm) + fe(sector) & fe(quarter))
r = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 12,
            shock = "L1.dd & mp", estimator = LagAugmented())
summarize(r, TLAHR(); level = 0.90)          # their two-way clustering, for comparison:
w = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 12,
            shock = "L1.dd & mp")
summarize(w, Clustered(:unit, :time); level = 0.90)
```

The exposure varies over time, so the synthetic weights of `TLAHR()` are the lagged
characteristic residualised within sector × quarter cells.

## Jeenas & Lagos (2024), reduced form

- Outcome: the investment rate in levels, with one lag of it as a control.
- Exposure: lagged stock turnover, interacted with an indicator of the bottom half of
  the liquidity distribution (``I_{L,t-1}``). Every term of the regression is interacted
  with it, the fixed effects included.
- Effects: firm, and 2-digit SIC × quarter.
- Shock: Jarociński–Karadi surprises.
- Their standard errors: clustered by firm and 3-digit SIC × quarter, 1990Q1–2016Q4.

With two exposures of the same shock (turnover × high and × low liquidity), `TLAHR()`
weights each by its exposure residualised on the other (see
[t-LAHR](../methods/tlahr.md)). The interactions of the fixed effects with the liquidity
indicator are fixed effects of the pairs, `fe(firm) & fe(lowliq)` and
`fe(sic2) & fe(quarter) & fe(lowliq)`; the firm × liquidity effects are nested in the
firms, so they satisfy `LagAugmented`'s requirement of unit effects. The interacted
slopes are products such as `lags(x, 1) & lowliq`. The lag of the outcome that
`LagAugmented` adds, `L1.x`, is the sum of the two interacted lags, so it is left out of
the regression (silently: the formula's columns span it).

```julia
f = @formula(leads(x) ~ lags(turn, 1) & highliq & mp + lags(turn, 1) & lowliq & mp +
             lags(x, 1) & highliq + lags(x, 1) & lowliq + lags(turn, 1) & highliq +
             lags(turn, 1) & lowliq + fe(firm) & fe(lowliq) +
             fe(sic2) & fe(quarter) & fe(lowliq))
r = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 12,
            shock = ["L1.turn & highliq & mp", "L1.turn & lowliq & mp"],
            estimator = LagAugmented())
```

## Blocked: instrumental variables

Two specifications need panel LP-IV, planned for version 0.2:

- **CFFS23.** The change in the one-year rate, interacted with four age × dividend groups,
  instrumented by the Gertler–Karadi proxy shock; firm effects and seasonal dummies;
  Driscoll–Kraay standard errors with 12 lags.
- **JL24, IV.** Tobin's q instrumented by the turnover × shock interactions.

Almuzara & Sancibrián give the t-LAHR variance for panel LP-IV in their supplement (§B),
but not its small-sample refinement.
