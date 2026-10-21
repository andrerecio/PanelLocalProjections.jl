# PanelLocalProjections.jl

Panel local projections (Jordà 2005): for each horizon ``h``, a panel regression of
``y_{i,t+h}`` on a shock. The package implements:

- the fixed-effects (within/LSDV) panel LP, with unit-, time- and two-way-clustered,
  heteroskedasticity-robust and Driscoll–Kraay standard errors;
- the **split-panel jackknife** correction of the Nickell bias of the fixed-effects
  estimator (Mei, Sheng & Shi, 2026, *Journal of International Economics*);
- the **cumulative-difference OLS** estimator, which avoids the within transformation
  (Ugarte-Ruiz, 2026);
- **t-LAHR** inference for panel LPs with aggregate shocks: time-clustered,
  lag-augmented, heteroskedasticity-robust, with an HC2 + Imbens–Kolesár small-sample
  refinement (Almuzara & Sancibrián, 2026).

!!! note "Status"
    Version 0.1.0. Panel LP-IV is planned for version 0.2; the roadmap is in `TODO.md` at
    the repository root.

Fixed-effect absorption, the OLS solver and all variance estimators are implemented in the
package. Their design takes inspiration from
[Regress.jl](https://github.com/gragusa/Regress.jl) and
[CovarianceMatrices.jl](https://github.com/gragusa/CovarianceMatrices.jl)
(Giuseppe Ragusa), [FixedEffectModels.jl](https://github.com/FixedEffects/FixedEffectModels.jl)
(Matthieu Gomez), and the Stata commands `reghdfe` (Sergio Correia) and `xtlp`
(Mei, Sheng & Shi). None of them is a dependency.

## Package

```@docs
PanelLocalProjections
```
