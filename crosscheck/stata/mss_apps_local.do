* The Romer–Romer applications of Mei, Sheng & Shi (2026) with exact two-way fixed
* effects, by xtlp (whose tfe demeaning iterates to convergence; the authors' R code makes
* one sweep). An independent check of the exact numbers on the replication page.
* Inputs are the local, gitignored replication data under docs/references/, so the output
* (test/fixtures/local/mss_apps_xtlp.csv) is local too. Values are unscaled (RR is ×7).
do "$PLP_ROOT/crosscheck/stata/common.do"

local app "$PLP_ROOT/docs/references/code/panel-lp-replication/applications"
capture confirm file "`app'/empirical_RR_f4_lnunemp_1980.csv"
if _rc {
    display as text "MSS replication data absent (local-only): skipped"
    exit
}
capture mkdir "$PLP_ROOT/test/fixtures/local"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/local/mss_apps_xtlp.csv", write replace
file write `fh' "app,method,h,term,b,se,N" _n
foreach v in UNEMP LNGDP {
    local file = cond("`v'" == "UNEMP", "lnunemp", "lngdp")
    local name = cond("`v'" == "UNEMP", "RR_unemp", "RR_gdp")
    import delimited using "`app'/empirical_RR_f4_`file'_1980.csv", clear asdouble
    xtset countdums halfyear
    forvalues h = 0/10 {
        foreach m in fe spj {
            * the authors keep singletons; ase: no small-sample factor (their variance)
            quietly xtlp f`h'`=lower("`v'")' crisis l1`=lower("`v'")' l2`=lower("`v'")' ///
                l3`=lower("`v'")' l4`=lower("`v'")' l1crisis l2crisis l3crisis l4crisis, ///
                tfe m(`m') vce(cluster countdums, ase) keepsingletons
            file write `fh' "`name'," (upper("`m'")) ",`h',CRISIS," %25.17e (_b[crisis]) "," ///
                %25.17e (_se[crisis]) "," %12.0f (e(N)) _n
        }
    }
}
file close `fh'
display as result "mss_apps_local: done"
