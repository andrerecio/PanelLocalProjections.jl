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
- [x] Guide, methods and cross-check pages (`crosscheck.md` generated from `crosscheck/README.md`); audit pages (xtlp, Ramey, tidyMacro); `deviations.md`
- [x] `replication/firm_level_wip.md` (OW20 and the JL24 reduced form; CFFS23 and the JL24 IV blocked on LP-IV)
- [x] Memory: the design absorbed in place in one buffer (within transform, kept columns and scores are views; no raw copy), responses built per horizon. Max RSS at scale: fixed effects 1.3 GB, `LagAugmented(ICLags())` 1.8 GB (was 2.6), `CappedLags()` 2.2 GB (was 3.2; live data with up to 10 added lag columns, not uncollected garbage)
- [x] Full cross-check re-run from an empty `crosscheck/ado` (manifest matches; every committed output byte-identical); strict docs build; version 0.1.0
- [x] Phase 5 gate: `Pkg.test()` on 1.12 and 1.10 (35709 tests), format, fixtures, strict docs build, cross-checks from an empty ado tree
- [x] Commit

## Phase 6: Panel LP-IV (v0.2)

Design: `docs/notes/design/lpiv_design.md` (gitignored). Method: A&S main §3.4 and supplement A.4.
Per horizon, 2SLS on the FE-absorbed sample. The fit stores the bread (X̂'X̂)⁻¹ and the scores X̂ ξ̂, where ξ̂ is the
structural residual. Every existing variance then applies unchanged, and `TLAHR(refine = false)` is A.4 exactly.

### 6a: IV core
- [x] Synthetic fixtures `fx_lpiv` (balanced) and `fx_lpiv_unbal`, using the eq. (14) DGP:
  ffr_t = 0.5 ffr_{t−1} + X_t + 0.6 Z_t and mp_t = X_t + 0.5 ν_t, with group dummies `g1`–`g3`
- [x] Stata: `ivreghdfe` pinned (`bfb5577`); `ivreg2`, `ranktest`, `xtivreg2`, `weakivtest` (SSC) in
  `setup_stata.do`; ado manifest re-baselined (new files only); `VERSIONS.md`
- [x] Formula: one `(endo ~ inst)` block of RHS terms (products, `lags`, `firstdiff`). Errors for:
  - `fe()` or a constant inside the block;
  - two blocks, or a block nested in a term;
  - an exogenous regressor listed as an instrument;
  - fewer instruments than endogenous columns.

  An instrument numerically collinear with a regressor is dropped from the first stage with a note.
- [x] Design: endogenous columns first and `PanelDesign.Z`; `shock` defaults to the endogenous columns. IV is
  rejected with `SplitPanelJackknife`, `CDOLS` and (until 6b) `LagAugmented`
- [x] Fit: `[y X Z]` absorbed in one buffer.
  - Screen: an absorbed endogenous shock makes the horizon `:singular`; absorbed instruments are dropped, and
    too few left gives `:underidentified`.
  - First stage by FWL: pivoted QR of the exogenous columns, then the instruments in order by Gram–Schmidt
    (the later of collinear instruments is dropped, with a note).
  - Second stage: bread (X̂'X̂)⁻¹, residual ỹ − X̃β, scores X̂ξ̂.
  - Identification check (`_identified`): endogenous regressors the instruments do not predict, or with
    collinear fitted values, make the horizon `:underidentified`.
- [x] `Robust`, `Clustered` (1- and 2-way), `DriscollKraaySE`, `CosineHAR` and `TLAHR(refine = false)` on IV
  fits; `TLAHR(refine = true)` errors until 6b; `show` lists the instrumented and the instruments
- [x] Identities and oracles (`test/test_lpiv.jl`):
  - coefficients and every variance equal to a dummy-variable 2SLS oracle;
  - (x ~ x) is OLS (coefficients, every variance and df, to round-off);
  - just-identified equals (Z'X)⁻¹Z'y;
  - several endogenous regressors.
- [x] Cross-checks (`test/test_crosscheck_p6.jl`):
  - `ivreghdfe` (robust, cluster id / t / two-way, `dkraay`), `xtivreg2, fe`, and several endogenous
    regressors (160 rows);
  - fixest `feols(y ~ w | fe | endo ~ inst)` (150 rows);
  - b 5e-15, se 2e-10.

  ivreghdfe's clustered factor with a fixed effect nested in the cluster differs from reghdfe's, as
  documented.
- [x] Adversarial review (2 lenses: code; econometrics and tests with 14 mutations). Fixed:
  - round-off fitted values reached `_ols` and became arbitrary coefficients; now `:underidentified`;
  - `synthetic()` mixed X̂ and the structural residual on IV fits; now IV fits store none.

  Tests added: over-identified 2SLS, collinear exogenous regressors, the instrument tolerance, per-horizon
  underidentification. The surviving mutations M8, M9, M12 and M14 are now killed. Non-IV fits are
  bit-identical to HEAD (180 checks).
- [ ] Phase 6a gate; commit

### 6b: lag augmentation and t-LAHR for IV
- [ ] `LagAugmented` with IV: the lags of the endogenous columns are endogenous, instrumented by the same lags of
  the instruments; the lags of y are exogenous
- [ ] `ICLags` on IV fits: the criterion on the reduced-form synthetic regression (documented extension)
- [ ] Synthetic series from X̂; `TLAHR(refine = false)` = A&S A.4; `TLAHR(refine = true)` on IV fits stays an
  error pointing to `refine = false` (user decision 2026-10-25, "fine for the moment"; the X̂ extension or A.4 as
  the IV default can be revisited)
- [ ] Monte Carlo of the A&S LP-IV design (Figures 1–2, panels d–f): refined vs A.4 vs DK vs two-way

### 6c: diagnostics, docs, release
- [ ] `firststage(r; vcov)`: first-stage coefficients, the Kleibergen–Paap rk Wald F with the chosen variance,
  and the Montiel Olea–Pflueger effective F (one endogenous column); cross-checks against `ivreg2`'s `widstat`
  and `weakivtest`
- [ ] Docs: `methods/lpiv.md`; CFFS23 / JL24 IV mappings in `replication/firm_level_wip.md`;
  Jordà–Schularick–Taylor trilemma (if the data are public)
- [ ] v0.2.0

## After v0.1
- [ ] Adversarial review of the Phase 4 review fixes and of the Phase 5 memory changes (deferred; started twice, stopped for time). Scope:
  - econometrics: the ω_t-weighted synthetic regression (CR2 by period), the ω-weighted criterion with its k(p) penalty and its rule for repeated regressions, residualising on every same-aggregate regressor (`_sameagg`), and the ownership rules (`_owners`);
  - code: `_shockweights`/`_projresid!` and the zero-weight rule, `_dropcollinear!`, the `_samecol` lag dedup, and the in-place `_within_fit` buffer (aliasing, the split-panel jackknife path);
  - tests and docs: mutation tests of the fixes, and every number changed in `methods/tlahr.md` and `audit/as_code.md`.

  Re-run the first review's reproductions and fuzzers on the current code. Fix the confirmed findings in a follow-up commit.
- Joint (sup-t) bands across horizons, built on the stored group score sums
- Negative horizons (pre-trends), lagged-FE interactions, a `panelshift` helper
- Firm-level replications (OW20, CFFS23, JL24) with user-supplied Compustat data
