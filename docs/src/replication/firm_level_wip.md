# Firm-level applications (work in progress)

Almuzara & Sancibrián (2026, §5) re-estimate three firm-level studies of monetary
policy with t-LAHR: Ottonello & Winberry (2020, OW20), Cloyne, Ferreira, Froemel &
Surico (2023, CFFS23) and Jeenas & Lagos (2024, JL24). All three use quarterly Compustat;
OW20 and JL24 also use CRSP, CFFS23 WorldScope, and OW20 Datastream incorporation dates.
These data are licensed (Compustat and CRSP through WRDS; Datastream and WorldScope from
LSEG/Refinitiv) and cannot be redistributed. The authors' replication packages contain
their code and results but not the raw data. This page therefore maps the specifications
to `panellp` calls for users who have those licences (**requires WRDS**, and LSEG for
Datastream/WorldScope). **None of them has been run on
the data**: the calls below were run only on synthetic panels of the same shape, to check
that they estimate. For a panel LP-IV on public data, see the
[US state example](../examples/us_states.md). The specifications are summarised from the papers
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

Source of the specification: the replication package of Ottonello & Winberry (2020),
"Financial Heterogeneity and the Investment Channel of Monetary Policy", *Econometrica*
88(6), 2473–2502 (the Econometrica replication package, June 2020). Its `construct_panel.do`
builds the sample from Compustat quarterly fundamentals, the CRSP–Compustat merge, CRSP
market values (for distance to default) and Datastream incorporation dates. Its
`results/` folder holds their estimates, the targets of a replication. The package
excludes the raw data "to respect potential copyright restrictions".

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

## Cloyne, Ferreira, Froemel & Surico (2023): panel LP-IV

- Outcome: ``X_{i,t+h} - X_{i,t-1}``, with X capital expenditure over beginning-of-period
  net PP&E (`ldiff(capex)`).
- Endogenous: the change in the one-year rate, ``\Delta R_t`` (end of quarter, +25 bp
  normalisation), interacted with four age × dividend groups at ``t - 1`` (`g1`–`g4`, dummies
  of the lagged group).
- Instrument: the Gertler–Karadi (2015) proxy-SVAR shock × the same groups.
- Effects: firm and quarter-of-year (seasonal). There are no time effects, so the group
  responses include general-equilibrium effects. The group dummies enter as controls (three
  of them, as the four sum to one).
- Their standard errors: Driscoll–Kraay with 12 lags, 1986Q1–2016Q4.

```julia
f = @formula(ldiff(capex) ~ (g1 & dR + g2 & dR + g3 & dR + g4 & dR ~
                             g1 & gk + g2 & gk + g3 & gk + g4 & gk) +
             g2 + g3 + g4 + fe(firm) + fe(qofyear))
r = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 20)
summarize(r, DriscollKraaySE(12); term = "g1 & dR", level = 0.90)   # their bands
firststage(r; vcov = DriscollKraaySE(12))                          # KP F, 4 endogenous
# A&S's re-estimation: lag augmentation and t-LAHR (A.4, CR0 by period)
ra = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 20,
             estimator = LagAugmented(CappedLags(4)))
summarize(ra, TLAHR(refine = false); term = "g1 & dR", level = 0.90)
```

With four endogenous shocks, `ICLags` weights the reduced form with the exposure of the
first one (`g1`), so `LagAugmented()` works too. `TLAHR(refine = true)` is not defined for
IV fits (see [panel LP-IV](../methods/lpiv.md)).

## Jeenas & Lagos (2024): instrumental variables

- Endogenous: ``\log q_{it}`` (Tobin's q), a firm-level variable.
- Instruments: ``\mathcal T_{i,t-1}\varepsilon^m_t`` and
  ``I_{L,t-1}\mathcal T_{i,t-1}\varepsilon^m_t`` (the turnover × shock interactions of the
  reduced form).
- Their first-stage F: 16.0 at h = 0 for the investment rate (fn. 50).

```julia
f = @formula(leads(x) ~ (logq ~ lags(turn, 1) & highliq & mp + lags(turn, 1) & lowliq & mp) +
             lags(x, 1) & highliq + lags(x, 1) & lowliq + lags(turn, 1) & highliq +
             lags(turn, 1) & lowliq + fe(firm) & fe(lowliq) +
             fe(sic2) & fe(quarter) & fe(lowliq))
r = panellp(f, data; id = :firm, time = :quarter, frequency = :quarter, horizon = 12,
            clusters = ((:sic3, :quarter),))
summarize(r, Clustered(:unit, (:sic3, :quarter)))        # their clustering
summarize(r, TLAHR(refine = false))                      # by period
firststage(r; vcov = TLAHR(refine = false))
```

The endogenous regressor is not an aggregate series times an exposure, so `ICLags` (whose
weights need one) does not apply; `LagAugmented(ConstantLags(p))` lags ``\log q`` as an
endogenous regressor, instrumented by the lags of the instruments.
