# Phase 6a cross-check with fixest: panel LP-IV, y_{i,t+h} on s_i ffr_t instrumented by
# s_i mp_t, with w and y_{i,t-1}, unit and time effects (2SLS on the fixed-effects design).
#     Rscript crosscheck/R/p6_fixest.R
# Writes test/fixtures/r/p6_fixest.csv; spec labels map to panellp/vcov calls in
# test/test_crosscheck_p6.jl. Every vcov is passed with its small-sample correction
# explicitly; `vcov_fix = FALSE` keeps two-way variances unrepaired (psd = :none).
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
  d[, `:=`(sffr = s * ffr, smp = s * mp, Ly = shiftcol(d, "y", -1))]
  for (h in 0:4) {
    d[, Fy := shiftcol(d, "y", h)]
    m <- feols(Fy ~ w + Ly | id + t | sffr ~ smp, d, fixef.rm = "none")
    emit(ds, "twfe_hetero0", h, m, vcov(m, vcov = "hetero", ssc = none))
    emit(ds, "twfe_cl_id0", h, m, vcov(m, cluster = ~id, ssc = none))
    emit(ds, "twfe_cl_t0", h, m, vcov(m, cluster = ~t, ssc = none))
    emit(ds, "twfe_cl_id_t0", h, m, vcov(m, cluster = ~id + t, ssc = none, vcov_fix = FALSE))
    emit(ds, "twfe_dk2_0", h, m, vcov(m, vcov = DK(2) ~ t, ssc = none))
  }
}
res <- rbindlist(out)
for (v in intersect(c("b", "se"), names(res))) set(res, j = v, value = sprintf("%.17g", res[[v]]))
fwrite(res, file.path(root, "test", "fixtures", "r", "p6_fixest.csv"))
cat(sprintf("p6_fixest: %d rows (fixest %s)\n", nrow(res), as.character(packageVersion("fixest"))))
