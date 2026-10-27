# Ugarte-Ruiz (2026)

*Panel local projections without fixed-effects: a cumulative-difference OLS estimator*,
BBVA Research. The paper reports its applications as figures only. The values below are
read off them, to about ±0.2.

## Cerra–Saxena (2008): currency crises (Fig. 9)

The data are the Mei–Sheng–Shi replication file `empirical_CS_f3.csv` (local only; see
the [Mei–Sheng–Shi page](mss2026.md)). Its outcomes are pre-built long differences
`cf{h}GRRT_WB` = ``y_{t+h} - y_{t-1}`` of log real GDP (× 100). The h = 0 outcome is
the growth rate ``g_t``, which is the lead of `l1GRRT_WB`. The specification is the
crisis dummy, four of its lags and four lags of growth, pooled with an intercept and
without country effects:

```julia
f = @formula(cf4GRRT_WB ~ CRISIS + l1CRISIS + l2CRISIS + l3CRISIS + l4CRISIS +
                          l1GRRT_WB + l2GRRT_WB + l3GRRT_WB + l4GRRT_WB)
panellp(f, data; id = :cntyNUM, time = :obs, horizon = 0,
        estimator = CDOLS(validate = false))      # a pre-built outcome column
```

| ``h`` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `panellp` | −1.07 | −2.53 | −3.34 | −4.10 | −4.49 | −4.65 | −5.18 | −4.91 | −5.11 | −5.32 | −5.69 |
| Fig. 9, CMOLS | −1.1 | −2.5 | −3.4 | −4.1 | −4.5 | −4.65 | −5.2 | −4.9 | −5.1 | −5.3 | −5.7 |

Every horizon is within 0.06 of the figure (a local test item checks 0.1).

## Mian–Sufi–Verner (2017): household and firm debt (Fig. 10)

Not reproduced. The paper does not state the CD-OLS outcome, its base period or how the
outcome lags enter for this application. With the MSS file `empirical_MSV_f2.csv`, the
closest variant tried (outcome ``y_{t+h} - y_t``, both debt shocks and four lags of each,
four lags of ``\Delta y``) matches Fig. 10 up to ``h \approx 3`` and drifts away after.
The difference is not only the CD-OLS choices: even Fig. 10's fixed-effects and jackknife
lines differ from FE and SPJ on the MSS file (household debt at ``h = 7``: −0.47 and
−0.65, against −0.387 and −0.558), whereas Fig. 9's match the Cerra–Saxena file. So the
paper's Mian–Sufi–Verner data or specification differ from the MSS file.

## Furceri et al. (2018): tariffs (Fig. 8)

Not attempted: the data are not part of any replication kit available here.
