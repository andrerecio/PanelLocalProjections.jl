* Smoke check of the cross-check toolchain (Phase 0 gate).
do "$PLP_ROOT/crosscheck/stata/common.do"
display c(stata_version) " " c(edition_real)
foreach c in require ftools reghdfe xtlp xtscc locproj xtspj xtarsim harreg {
    which `c'
}
sysuse auto, clear
reghdfe price weight, absorb(rep78) vce(cluster rep78)
* the HC2 + Bell-McCaffrey df machinery used in Phase 4 (stores e(adj_df))
sysuse auto, clear
regress price mpg weight, vce(hc2, dfadjust)
matrix list e(adj_df)
display as result "smoke: ok"
