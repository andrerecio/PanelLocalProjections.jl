* Install the Stata packages used by the cross-checks into crosscheck/ado (gitignored).
* GitHub packages are pinned by commit SHA; SSC packages cannot be pinned, so their files
* are checksummed by crosscheck/ado_manifest.py (see crosscheck/VERSIONS.md).
do "$PLP_ROOT/crosscheck/stata/common.do"
capture mkdir "$PLP_ROOT/crosscheck/ado"
capture mkdir "$PLP_ROOT/crosscheck/ado/plus"
capture mkdir "$PLP_ROOT/crosscheck/ado/personal"

local gh "https://raw.githubusercontent.com"
net install require, from("`gh'/sergiocorreia/stata-require/7ee750b12f9c51ec8a879445cbf39357c308bdee/src/") replace
net install ftools,  from("`gh'/sergiocorreia/ftools/7b3663e49ea5c5b81638c55be29edf416e68e8b7/src/") replace
net install reghdfe, from("`gh'/sergiocorreia/reghdfe/29b2b203f534413369d675d0c0f674d61291f536/src/") replace
net install xtlp,    from("`gh'/shenshuuu/panel-local-projection-stata/cd3fa2468ada916d9c05323840b0a2abd48e0173/") replace
* panel LP-IV (Phase 6): ivreghdfe needs ivreg2 and ranktest (SSC, below)
net install ivreghdfe, from("`gh'/sergiocorreia/ivreghdfe/bfb5577a6dbdfb029ab4ab6a7e93f7257a827b42/src/") replace

* weakivtest needs avar and distinct
foreach p in xtscc locproj xtspj xtarsim harreg ivreg2 ranktest xtivreg2 weakivtest avar distinct {
    ssc install `p', replace
}

ftools, compile
* reghdfe 6 compiles its Mata code on first use (it has no `compile` subcommand).
display as result "setup_stata: done"
