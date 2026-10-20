# crosscheck/

Scripts that compute the reference numbers the Julia tests compare against: Stata through
PyStata, and R. The scripts are committed. The local toolchain and the run outputs
(`.venv/`, `crosscheck/ado/`, `crosscheck/out/`) are not.

## Setup (once per machine)

```bash
uv venv .venv --python 3.13
uv pip install --python .venv/bin/python -r crosscheck/requirements.txt
.venv/bin/python crosscheck/run_do.py crosscheck/stata/setup_stata.do --log crosscheck/out/setup_stata.log
.venv/bin/python crosscheck/ado_manifest.py --check     # installed files match the committed checksums
.venv/bin/python crosscheck/run_do.py crosscheck/stata/smoke.do
Rscript crosscheck/R/smoke.R
```

`setup_stata.do` installs every package into the project-local tree `crosscheck/ado/`.
`common.do` points `sysdir PLUS/PERSONAL` there, so the user's own Stata setup is never
touched. Versions and pins are in `VERSIONS.md`.

## Running a cross-check

```bash
.venv/bin/python crosscheck/run_do.py crosscheck/stata/<script>.do --log crosscheck/out/<script>.log
Rscript crosscheck/R/<script>.R
```

Each do-file starts with `do "$PLP_ROOT/crosscheck/stata/common.do"`. `run_do.py` sets
`PLP_ROOT` to the repository root and runs the do-file from its own directory, in a fresh
Stata session each time (~1–2 s). The exit code is the Stata return code.

## Conventions every script follows

- `set type double` (in `common.do`) and `import delimited …, asdouble`. Stata's float
  storage otherwise costs about 1e-7 relative error.
- Inputs are the committed synthetic panels in `test/fixtures/synthetic/`, or local-only
  data under `docs/references/` whose outputs go to `test/fixtures/local/` (gitignored).
- Outputs are written at full precision: Stata `%21x` or `%21.0g`; R `format(x, digits = 17)`
  or `data.table::fwrite`. Committed outputs go to `test/fixtures/{stata,r}/` and are
  listed with their provenance in `test/fixtures/README.md`.
- Re-runs are compared numerically (≤ 1e-12), ignoring header lines.
- fixest ≥ 0.13: `ssc(K.adj = , G.adj = )`, explicit `fixef.rm`, explicit `vcov =`.
- MSS `LP_panel` needs a base `data.frame` (`as.data.frame()` first). It is reliable only in
  its `H = 0` mode with pre-built leads and lags.

## Command ↔ Julia call map

Filled in phase by phase. Each row gives the exact `panellp`/`vcov` options that reproduce
a reference command: `singletons`, `timelags`, `ssc`, `psd`, `absorb`, `dist`.

| Reference command | Julia call | Phase |
|---|---|---|
| (Phase 1+) | | |
