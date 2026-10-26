# Phase 6b cross-check with fixest: lag-augmented panel LP-IV. At horizon h, p = min(h, 2)
# lags of y (exogenous) and of s_i ffr_t (endogenous, instrumented by the same lags of
# s_i mp_t), unit and time effects. Clustered by period without small-sample correction,
# this is the t-LAHR variance of Almuzara & Sancibrián's LP-IV (supplement A.4).
#     Rscript crosscheck/R/p6b_fixest_la.R
# Writes test/fixtures/r/p6b_fixest_la.csv; spec labels map to panellp/vcov calls in
# test/test_crosscheck_p6.jl.
suppressPackageStartupMessages({ library(fixest); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
setFixest_nthreads(1)

shiftcol <- function(d, v, k) {   # v at (id, t + k), NA where the cell is absent
  key <- d[, .(id, t = t - k, val = get(v))]
  key[d, on = .(id, t), val]
}
none <- ssc(K.adj = FALSE, G.adj = FALSE)
out <- list()
emit <- function(ds, spec, h, m, V) {
  b <- coef(m); se <- sqrt(diag(V))
  for (v in names(b)) out[[length(out) + 1]] <<- data.table(
    dataset = ds, spec = spec, h = h, term = v, b = b[[v]], se = se[[v]], nobs = nobs(m))
}
for (ds in c("fx_lpiv", "fx_lpiv_unbal")) {
  d <- fread(file.path(root, "test", "fixtures", "synthetic", paste0(ds, ".csv")))
  setkey(d, id, t)
  d[, `:=`(sffr = s * ffr, smp = s * mp)]
  for (l in 1:2) {
    d[, paste0("Ly", l) := shiftcol(d, "y", -l)]
    d[, paste0("Lsffr", l) := shiftcol(d, "sffr", -l)]
    d[, paste0("Lsmp", l) := shiftcol(d, "smp", -l)]
  }
  for (h in 0:4) {
    p <- min(h, 2)
    d[, Fy := shiftcol(d, "y", h)]
    lags <- function(v) if (p == 0) "" else paste0(" + ", paste0(v, 1:p, collapse = " + "))
    f <- as.formula(paste0("Fy ~ w", lags("Ly"), " | id + t | sffr", lags("Lsffr"),
                           " ~ smp", lags("Lsmp")))
    m <- feols(f, d, fixef.rm = "none")
    emit(ds, "la_twfe_cl_t0", h, m, vcov(m, cluster = ~t, ssc = none))
    emit(ds, "la_twfe_hetero0", h, m, vcov(m, vcov = "hetero", ssc = none))
  }
}
res <- rbindlist(out)
for (v in intersect(c("b", "se"), names(res))) set(res, j = v, value = sprintf("%.17g", res[[v]]))
fwrite(res, file.path(root, "test", "fixtures", "r", "p6b_fixest_la.csv"))
cat(sprintf("p6b_fixest_la: %d rows (fixest %s)\n", nrow(res), as.character(packageVersion("fixest"))))
