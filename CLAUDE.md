# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

PanelLocalProjections.jl estimates panel local projections: for each horizon h it runs a panel regression of y_{i,t+h} on a shock. Methods:

- **Mei–Sheng–Shi (2026, *JIE*):** the Nickell bias of the FE estimator and the split-panel jackknife (SPJ) correction. Their Stata command is `xtlp`.
- **Ugarte-Ruiz (2026):** cumulative-difference OLS (CD-OLS). It regresses y_{t+h} − y_{t−1} on s_t and s_{t−1} by pooled OLS, without unit FE.
- **Almuzara–Sancibrián (2026):** **t-LAHR** inference. It is time-clustered, lag-augmented and heteroskedasticity-robust, with an HC2 + Imbens–Kolesár df refinement. `starting_idea.md` calls it a "HAR estimator", but it is **not** HAR: Driscoll–Kraay, EWC and unit/two-way clustering are its *comparators*.

The phased roadmap, with exit criteria per phase, is in `TODO.md`. Panel LP-IV is planned for v0.2.

## Hard rules

- **Autonomy.**
  - The package must not depend on Regress.jl, CovarianceMatrices.jl, FixedEffects.jl/FixedEffectModels.jl, Vcov.jl or ShiftedArrays. FE absorption, OLS and every variance estimator are implemented in-house.
  - Allowed deps: StatsModels, Tables, StatsAPI, StatsFuns, RecipesBase and stdlibs; Makie only as a weakdep.
  - Code adapted from Regress.jl/CovarianceMatrices.jl/FixedEffectModels.jl, from `~/Workspace/LocalProjections.jl` (the user's fork; its `src/terms.jl` is the source of the LHS formula DSL), or from MIT-licensed reference code (MSS R, `xtlp`) carries a source comment naming the project and file.
- **Licences.**
  - The Almuzara–Sancibrián MATLAB/R code has **no licence** and `pLP` is CC BY-NC-SA. Implement those methods from the papers' formulas only; never translate the code.
  - Never commit third-party code or data.
  - Outputs of any tool computed on *our own* synthetic panels may be committed, including A&S `panel_LP.R` outputs. Each such file is listed in `test/fixtures/README.md`.
- **Local-only reference material.** It is gitignored, read-only, and never formatted or committed:
  - `docs/paper/`: papers as markdown, with figure JPEGs.
  - `docs/references/code/`:
    - `panel-lp-replication/` is the current MSS replication kit (`applications/LP_panel_all.r` + CSVs + `data_preparation/*.dta`);
    - `nickell_bias/*.zip` is an older snapshot;
    - `PanelLocalProjections/` is the A&S code;
    - `ramey_govinfrastructure_replication/` is Ramey's code.
  - `docs/notes/`: planning research reports, the design draft, 52 adversarial review findings, and `design/approved_plan.md` (**the authoritative design**). Claims tagged [I]/[INF]/[W] there are unverified.
  - The markdown conversions have garbled spots: MSS main.md:132 and eq. (22); A&S Table 3 row labels. Check the notes before trusting a garbled formula.
- **Git.** Work on `dev`. Commit only when a phase gate passes (see `TODO.md`). Never push.

## Commands

```bash
# Package environment (Julia 1.12; manifests are per version: Manifest-v1.12.toml / Manifest-v1.10.toml, gitignored)
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using Pkg; Pkg.test()'                      # full suite: TestItemRunner items + Aqua
julia +lts --project -e 'using Pkg; Pkg.test()'                 # compat floor 1.10 (first time: Pkg.resolve() before Pkg.test())

# Run a tagged subset. Test deps are in test/Project.toml; TestEnv lives in the dev env, stacked onto LOAD_PATH.
julia --project -e 'push!(LOAD_PATH, "dev"); using TestEnv; TestEnv.activate(); using TestItemRunner; @run_package_tests filter = ti -> (:fixtures in ti.tags)'

julia --project=dev dev/format.jl [--check]                     # SciML formatting; --check exits 1 if anything would change
julia --project=dev dev/make_fixtures.jl [--check]              # (re)generate / verify test/fixtures/synthetic/*.csv
julia --project=docs -e 'using Pkg; Pkg.instantiate()' && julia --project=docs docs/make.jl

# Cross-checks (see crosscheck/README.md; the .venv and crosscheck/ado are local-only)
uv venv .venv --python 3.13 && uv pip install --python .venv/bin/python -r crosscheck/requirements.txt
.venv/bin/python crosscheck/run_do.py crosscheck/stata/setup_stata.do       # installs pinned ados into crosscheck/ado
.venv/bin/python crosscheck/ado_manifest.py --check                          # installed ados unchanged (SSC drift detector)
.venv/bin/python crosscheck/run_do.py crosscheck/stata/<file>.do --log crosscheck/out/<file>.log
Rscript crosscheck/R/<file>.R
```

- **Test tags:**
  - `:core`, `:aqua`, `:fixtures`, `:panel`, `:formula`, `:absorb`, `:ols`, `:vcov`, `:api`, `:spj`, `:mc` (Monte Carlo, ~20 s), `:crosscheck` (with `:stata`/`:r`); `:cdols` and `:tlahr` arrive with later phases;
  - `:mc_full` runs only with `PLP_MC_FULL=1`;
  - `:local_mss` items are skipped automatically when their local data is absent (`test/runtests.jl`).
- **Formatting:** files that use `@testitem … setup = [...]` start with `#! format: off`, because the SciML style mangles them.
- **Kaimon.** When the `mcp__kaimon__*` tools are available, use the live REPL (session for this repo via `start_session(project_path = …)`, then `ex`) instead of spawning `julia` for exploration. Revise picks up `src/` edits automatically, but a struct change in `src/types.jl` needs `manage_repl(command = "restart")`. Load test-only packages and the fixtures with `pushfirst!(LOAD_PATH, joinpath(pkgdir(PanelLocalProjections), "dev"))` and `include("test/support/PLPDGPs.jl")` (`read_csv`). `run_tests` runs the whole suite, because TestItemRunner ignores its pattern argument.
- **Scale benchmark:** `julia --project=. -e 'push!(LOAD_PATH, "dev"); include("dev/benchmark_scale.jl")'` (about 15 s; Phase 1 reference: 12.8 s, 14 sweeps per horizon and 1.6 GB max RSS for 1.1M rows and 21 horizons).

## Architecture (Phase 1 implemented; later phases per `docs/notes/design/approved_plan.md`)

One module and no submodules. Every struct lives in `src/types.jl`. Include order: `utils` → `types` → `panelindex` → `formula` → `design` → `absorb` → `ols` → `scores` → `lagselect` → `estimators` → `panellp` → `vcov` → `summary` → `plotting`; Makie lives in `ext/`.

**Pipeline.** `panellp(f, data; id, time, horizon, estimator, subset, sample, singletons, absorb, clusters)` runs these stages:
1. `PanelIndex`: unit codes, integer calendar periods (Dates need `frequency`), and a (unit, period) → row lookup.
2. Formula rewrite into our own term types (`leads`, `ldiff`, `anchor`, `cumul`, `lags`, `firstdiff`, `fe`).
3. `PanelDesign`: a NaN float table, with the RHS built once by gathers.
4. Per-horizon mask, then the singleton policy.
5. FE absorption (`absorb.jl`), then an absorbed-column screen and pivoted QR (`ols.jl`).
6. The estimator's fit: `WithinFE` / `SplitPanelJackknife` / `CDOLS` / `LagAugmented`.
7. `HorizonFit`, which stores per-group score sums, the bread and the synthetic series. It does **not** store the design matrix.

`vcov(spec, r)` is computed post-estimation (`Robust`, `Clustered`, `DriscollKraaySE`, `TLAHR`, `CosineHAR`); `summarize` and plotting work from its result.

**Invariants:**
- Never shift by row position. Never drop missing values before shifting. Use NaN, not `missing`.
- Re-absorb every FE on each horizon's sample and inside each SPJ half. Screen each half for absorbed columns as `_ols` does (`_absorbed_cols`): rescaling a round-off column turns it into an arbitrary coefficient. A failing half makes the horizon `:spj_half_singular`, and with no `:ok` horizon `panellp` throws.
- A shock absorbed by the FE (e.g. an aggregate shock with `fe(time)`) is never a zero coefficient. It is an error when it happens at every horizon. At some horizons only, those fits are `:singular` with NaN coefficients.
- Problems at some horizons (empty, singular, omitted regressors, non-convergence) are reported once per `panellp` call (`_fitnotes`) and listed by `show(r)`. Never use `maxlog`, which silences later calls for the rest of the session.
- Undefined variances (a single cluster or period, no residual df) are NaN with a warning, never an exception. Negative variances give NaN standard errors.
- DK, EWC and t-LAHR aggregate scores by calendar period.

**Compat switches** exist only where a published number depends on them:
- `OneSweep()`: the one-pass two-way demeaning of the MSS R code;
- `ssc = :mss`;
- `SplitPanelJackknife(split = :floor)`;
- `CappedLags`: the A&S code's lag rule;
- `timelags = :position`: xtlp/xtscc DK.

The map from each reference command to the matching Julia options is in `crosscheck/README.md`.

## Cross-check gotchas

- **Stata.**
  - Every do-file starts with `do "$PLP_ROOT/crosscheck/stata/common.do"`. It sets `set type double` and redirects `sysdir PLUS/PERSONAL` to `crosscheck/ado`, so the user's own ados are never used.
  - Use `import delimited …, asdouble`.
  - reghdfe 6 has no `compile` or `version` subcommands (they replay and fail with r(301)), and its DK `bw` = lags + 1.
  - Read ado sources as files under `crosscheck/ado/plus/<letter>/`, not with `viewsource`.
- **MSS application data.** The CSVs are float32-origin, so expect about 1e-6 agreement. MSS `LP_panel` is reliable only in its `H = 0` mode with pre-built leads and lags, and it needs a base `data.frame`.
- **fixest ≥ 0.13.** Use `ssc(K.adj = , G.adj = )`, set `fixef.rm` explicitly and pass `vcov =` explicitly. Its default `vcov_fix = TRUE` clips negative eigenvalues of V for every vcov type: that is `psd = :raw` (two-way).
- **reghdfe two-way clustering with `noabsorb`.** reghdfe clips the slope block on its own. Stata can then refuse the repaired V when posting it ("nonsymmetric or highly singular") and report zero standard errors.
- **A&S `panel_LP.R`.** It zero-fills rows when the FE columns are character, and its `Date` branch is dead code. Pass numeric ids and periods.
- **Fixture DGPs** (`test/support/PLPDGPs.jl`) must use explicit left-to-right loops for reductions. `sum` reassociates, and `Pkg.test` runs with `--check-bounds=yes`, which changed results in the last ulp.
