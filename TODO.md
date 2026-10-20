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
- [ ] `CDOLS(validate = true)` and its validation rules
- [ ] xtarsim panels (committed); regress / xtscc / locproj / xtreg / xtlp cross-checks; fixest
- [ ] Analytic plims and a negative test with heterogeneous shock means

## Phase 4: t-LAHR (Almuzara–Sancibrián)
- [ ] `LagAugmented(ICLags | CappedLags | ConstantLags)`; synthetic AIC/BIC (B.3/B.4)
- [ ] `TLAHR(refine)`: CR0 by time; HC2 on synthetic regressors (thin SVD); IK ν with contrast e_k; exposure residualisation
- [ ] `CosineHAR`, `synthetic`, `lagorder`
- [ ] Cross-checks: xtlp `vce(cluster t, ase)`; Stata `regress, vce(hc2, dfadjust)` (`e(adj_df)`); harreg; A&S `panel_LP.R` grid
- [ ] Monte Carlo smoke run of the stylized Table 1 design (`:mc_full` locally)

## Phase 5: Documentation, audit, release
- [ ] Guide, methods and cross-check pages; audit pages (MSS R, A&S code, xtlp, Ugarte-Ruiz, Ramey, tidyMacro); `deviations.md`
- [ ] `replication/firm_level_wip.md` (OW20 and the JL24 reduced form; CFFS23 and the JL24 IV blocked on LP-IV)
- [ ] Full cross-check re-run from an empty `crosscheck/ado`; strict docs build; version 0.1.0

## After v0.1
- Panel LP-IV, with t-LAHR for IV (A&S supp. §B) and ivreghdfe/xtivreg2 cross-checks
- Joint (sup-t) bands across horizons, built on the stored group score sums
- Negative horizons (pre-trends), lagged-FE interactions, a `panelshift` helper
- Firm-level replications (OW20, CFFS23, JL24) with user-supplied Compustat data
