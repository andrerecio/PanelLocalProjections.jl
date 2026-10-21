# Roadmap: PanelLocalProjections.jl v0.1

Work happens on `dev`, with one commit per phase once its gate passes; nothing is pushed.

**Every phase gate requires all of the following:**
- `Pkg.test()` is green on Julia 1.12 and on `julia +lts` (1.10), with Aqua's checks all enabled.
- `julia --project=dev dev/format.jl --check` is clean.
- `julia --project=dev dev/make_fixtures.jl --check` passes.
- The phase's cross-check scripts reproduce the committed fixtures numerically.
- The docs build.
- CLAUDE.md, this file and CHANGELOG.md are updated.

The full design and the rationale for every convention are in the approved plan. Local copy: `docs/notes/design/approved_plan.md` (gitignored).

## Phase 0: Bootstrap
- [x] `dev` branch, `.gitignore` (docs/paper, docs/references, docs/notes, starting_idea.md are local-only)
- [x] Project.toml (`julia = "1.10"`), module skeleton, versioned manifests (`Manifest-v1.12.toml`, `Manifest-v1.10.toml`)
- [x] `test/` (TestItemRunner + standalone `test/Project.toml`, Aqua), `dev/` (format `--check`, `make_fixtures --check`)
- [x] Synthetic fixtures, `test/support/PLPDGPs.jl`
- [x] docs skeleton, GitHub workflows (CI on lts and 1; check-only Format; Documentation; CompatHelper; TagBot; dependabot)
- [x] README, CITATION.bib, CHANGELOG, CLAUDE.md
- [x] `crosscheck/`: PyStata runner, uv venv, project-local Stata ado dir, `setup_stata.do`, smoke checks, `VERSIONS.md`, ado manifest
- [x] Read the installed sources of xtlp, reghdfe, xtscc, harreg, locproj, xtarsim and xtspj; claims verified (`docs/notes/stata_sources/`, local)

**Exit:** all of the above; the Stata and R smoke checks pass.

## Phase 1: Core FE panel LP and inference
- [x] `PanelIndex`: calendar periods, Dates with `frequency`, gap statistics, gathers
- [x] Formula terms: `leads`, `ldiff`, `anchor`, `cumul`, `lags`, `firstdiff`, `fe`; `@formula` re-exported
- [x] Rewrite helpers by own function, accepting known foreign helpers (error on any other same-named function)
- [x] Design construction; masks (`subset`, `sample`, singletons)
- [x] FE absorption: one-way exact, alternating projections with guarded Irons–Tuck, `OneSweep`; df_a by reghdfe's pairwise rule; nesting
- [x] OLS: absorbed-column screen before equilibration; pivoted QR; absorbed/collinear-shock errors
- [x] `WithinFE`; `HorizonFit` with group score sums (global codes)
- [x] `Robust`, `Clustered` (1/2-way, `ssc` `:none`/`:stata`/`:mss`, `psd` `:stata`/`:raw`), `DriscollKraaySE` (`:nw94`, `:as2026`, `:mss2026`, Int; `:calendar`/`:position`)
- [x] `summarize`/`IRFSummary` (Tables source), `show`, RecipesBase recipe, Makie extension (`irfplot`, `irfplot!`)
- [x] Cross-checks: xtlp, reghdfe, xtscc, xtreg, regress (`p1_fe.do`, h = 0..6, 2380 rows incl. derived edge-case panels); fixest (incl. a non-PSD two-way case) and MSS `LP_panel` FE (R): b within 5e-13, se within 2e-10
- [x] Benchmark (`dev/benchmark_scale.jl`): WithinFE, firm + sector×quarter FE, 1.10M rows, H = 20: 12.8 s (14 sweeps per horizon), 1.6 GB max RSS
- [x] Adversarial review (4 lenses, 46 findings, verified): fixes applied (see CHANGELOG)
- [x] Phase 1 gate; commit
- Not covered by the test suite: the Makie extension (exercised by the docs build only)

## Phase 2: Split-panel jackknife (Mei–Sheng–Shi)
- [x] `SplitPanelJackknife(split = :ceil | :floor)`; halves re-absorbed, never dropping rows, screened for absorbed columns; D-score variance
- [x] Cross-checks: xtlp `m(spj)` (22 specs, b 1e-13, se 2e-11 tfe), xtspj (balanced, even n_h; 6e-14), MSS `LP_panel(method = "SPJ")` incl. `te = T` and the eigens repair (1e-14); hand-split identity and explicit-dummy variance on fx_unbalanced (`:ceil`, `:floor`)
- [x] MSS applications (local): exact n per h, Table 1 peaks to 3 decimals, the authors' notebook calls at about 1e-12
  (the EPS-extracted paths are not used: they come from the same R code, which is reproduced exactly).
  Deviations from the plan: BVX uses `psd = :raw` (the plan said `:stata`), which is what `eigens = T` does;
  the xtlp help examples on the `.dta` files are deferred (no `.dta` outputs exist yet)
- [x] xtlp `tfe` on the Romer–Romer data (local) confirms the exact two-way FE/SPJ columns of the replication page
- [x] Monte Carlo smoke run of the toy design (ρ = 0.8, N = 30, T = 60; 300 replications; FE and SPJ means, coverage, SPJ SE calibration at h = 1, 5)
- [x] Docs: `methods/spj.md`, `replication/mss2026.md` (one-sweep, split rule, eigens repair), `audit/mss_r.md`
- [x] Adversarial review (4 lenses, verified findings) addressed
- [x] Phase 2 gate; commit

## Phase 3: Cumulative-difference OLS (Ugarte-Ruiz)
- [x] `CDOLS(validate = true)` and its validation rules
- [x] xtarsim panels (committed; moments tested; the native DGP at SNR 20 reproduces them); regress / xtscc / locproj / xtreg / xtlp cross-checks, also on the unbalanced `fx_unbalanced` (b 1e-12, se 5e-14); fixest (4e-14)
- [x] Analytic plims (β_h, γ_h, omitted s_{t−1}) and a negative test with heterogeneous shock means
- [x] Monte Carlo of robust vs clustered standard errors (the paper's over-coverage)
- [x] Local: Cerra–Saxena CD-OLS reproduces Fig. 9 (within 0.06); Mian–Sufi–Verner not reproduced (specification not stated; Fig. 10's FE/SPJ lines also differ from the replication file)
- [x] Docs: `methods/cdols.md`, `replication/ugarte2026.md`, `audit/ugarte.md`
- [x] Adversarial review (4 lenses, verified findings) addressed: validator rules made structural, docs and Monte Carlo corrected, unbalanced cross-check added
- [x] Phase 3 gate; commit

## Phase 4: t-LAHR (Almuzara–Sancibrián)
- [x] `LagAugmented(ICLags | CappedLags | ConstantLags)`; synthetic AIC/BIC (B.3/B.4)
- [x] `TLAHR(refine, weights)`: CR0 by time; HC2 on synthetic regressors (thin SVD); IK ν with contrast e_k; exposure residualisation
- [x] `CosineHAR`, `synthetic`, `lagorder`
- [x] Cross-checks: A&S `panel_LP.R` grid (672 rows, 2e-12 / 2e-9); Stata `reghdfe` (the lag-augmented regression, CR1 by time; 2e-14, unbalanced se 5e-12), `regress, vce(hc2, dfadjust)` (`e(adj_df)`) and harreg on the synthetic regression (3e-14).
  The plan's xtlp `vce(cluster t, ase)` check is covered by reghdfe and by the identity `TLAHR(refine = false)` = `Clustered(:time; ssc = :none)`
- [x] Monte Carlo smoke run of the stylized Table 1 design (`:mc`); paper-size `:mc_full` run locally
- [x] Scale benchmark: `LagAugmented(ICLags())` 19 s and `TLAHR()` 0.7 s on 1.1M rows and 21 horizons, max RSS 2.6 GB; `CappedLags()` (up to 5 lags each) 25 s, 3.2 GB. Both above the 2 GB target of the fixed-effects run (1.6 GB): transient copies of the design in the criterion's six h = 0 fits and in every fit
- [x] Docs: `methods/tlahr.md`, `audit/as_code.md` (bugs B1–B4 verified by running the code). `audit/tidymacro.md` moves to Phase 5 (its findings need a source check)
- [x] Adversarial review (4 lenses, 35 findings, verified) addressed: rank-revealing per-period residualisation and exact zero weights (a crash in plain fits), weights residualised on every same-aggregate regressor and structural ownership, synthetic regression weighted by ω_t (CR2 by period) and ω-weighted criterion, numeric lag deduplication and lag priority in the collinearity screen, closed-form ν linear in T, messages, tests (CR2 oracle, paired Monte Carlo assertions), docs (sensitivity section)
- [x] Phase 4 gate: `Pkg.test()` on 1.12 and 1.10 (35709 tests), format, fixtures, docs build, p4 R/Stata fixtures byte-identical, ado manifest
- [x] Commit
- [ ] Later (research): exact CR2 by period on the panel regression, an alternative to the synthetic regression when unit-level controls or lags of time-varying exposures make it conservative (`methods/tlahr.md`)

## Phase 5: Documentation, audit, release
- [ ] Guide, methods and cross-check pages; audit pages (xtlp, Ramey, tidyMacro); `deviations.md`
- [ ] `replication/firm_level_wip.md` (OW20 and the JL24 reduced form; CFFS23 and the JL24 IV blocked on LP-IV)
- [ ] Memory: fewer transient copies of the design per fit (`_within_fit`: Xraw, A, Xt, Xk, the scores and `_eqr`'s copy), to bring `LagAugmented(ICLags())` (2.6 GB) and `CappedLags()` (3.2 GB) at scale under 2 GB
- [ ] Full cross-check re-run from an empty `crosscheck/ado`; strict docs build; version 0.1.0

## After v0.1
- Panel LP-IV, with t-LAHR for IV (A&S supp. §B) and ivreghdfe/xtivreg2 cross-checks
- Joint (sup-t) bands across horizons, built on the stored group score sums
- Negative horizons (pre-trends), lagged-FE interactions, a `panelshift` helper
- Firm-level replications (OW20, CFFS23, JL24) with user-supplied Compustat data
