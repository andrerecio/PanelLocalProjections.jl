# Fixed-effect absorption

Fixed effects are removed by the within transformation on each horizon's own sample. The
samples differ across horizons because leads run out at the end of each unit's series.

- **One dimension.** Each column is demeaned by group once, which is exact.
- **Several dimensions** ([`AlternatingProjections`](@ref)). Gauss–Seidel sweeps demean
  the dimensions in formula order. The limit is the residual of the projection on all the
  fixed-effect dummies (von Neumann–Halperin). After every two sweeps an Irons–Tuck step
  extrapolates along the last two updates. That direction lies in the span of the fixed
  effects, so the limit is unchanged. The step is skipped when the second difference is
  negligible relative to the first. Iteration stops when a sweep moves no element by more
  than `tol · max|v₁|` (default `tol = 1e-12`). Here `v₁` is the column after its first
  sweep, which has already removed the column's level, so the rule depends on neither
  the units nor the level of the data. Rescaling a variable, or adding a constant to it,
  changes neither the estimates nor the number of sweeps. On slowly mixing panels, such
  as short rotating spells, the remaining error can exceed the last step by a factor
  ``1/(1-\rho)``, where ``\rho`` is the contraction rate of the sweeps. Non-convergence
  within `maxiter` sweeps is reported, with the horizons affected.
- **[`OneSweep`](@ref).** One pass per dimension. With `fe(id) + fe(t)` this is the
  unit-then-time demeaning of the Mei–Sheng–Shi replication code, which is exact only on
  balanced samples. It exists to reproduce their published numbers.

## Absorbed degrees of freedom

``df_a`` is the rank of the fixed-effect dummies, counted as reghdfe does:

- the levels of the first dimension;
- for each later dimension, its levels minus ``\max_{i<j} m_{ij}``. Here ``m_{ij}`` is the
  number of connected components of the graph linking the groups of dimensions ``i`` and
  ``j``.

This is exact for one or two dimensions. With three or more it is an upper bound, and so
gives conservative standard errors; reghdfe reports the same bound. The test suite checks
``df_a`` against the rank of the dummy matrix (one and two dimensions, including a
sample with two disconnected blocks) and against reghdfe's `e(df_a)` (one, two and three
dimensions, nested clusters).

## Absorbed and collinear regressors

A regressor is *absorbed* and omitted when its within-transformed squared norm is at most
1e-8 times its pre-absorption centred squared norm (a norm ratio of 1e-4): its
coefficient is `NaN` and it is not counted in ``K``. The test runs before any rescaling,
which could otherwise inflate round-off into a spurious unit-norm column.

The remaining columns are scaled to unit norm and factorised by a pivoted QR. Columns
beyond the numerical rank (``|R_{jj}| ≤ 10^{-9}|R_{11}|``) are omitted as collinear.
Omitted regressors are reported once per call, with the horizons where they were omitted.

A *shock* that is absorbed or collinear in a horizon's sample makes that horizon
`:singular`: its coefficients are `NaN`, and the horizons are reported. This happens, for
example, with a staggered treatment dummy whose later horizons only see untreated
periods. When the shock is not identified at any horizon `panellp` throws an error, as
for an aggregate shock with time effects.

The design follows FixedEffectModels.jl (M. Gomez) and Regress.jl (G. Ragusa); the
acceleration follows fixest (L. Bergé) and reghdfe (S. Correia).
