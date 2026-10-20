# Phase 3 cross-check with fixest: cumulative-difference OLS (Ugarte-Ruiz 2026),
# y_{i,t+h} − y_{i,t−1} on s_it and s_{i,t−1} (pooled, with an intercept).
#     Rscript crosscheck/R/p3_fixest.R
# Writes test/fixtures/r/p3_fixest.csv; spec labels map to panellp/vcov calls in
# test/test_crosscheck_p3.jl.
suppressPackageStartupMessages({ library(fixest); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
setFixest_nthreads(1)

shiftcol <- function(d, v, k) {   # v at (id, t + k), NA where the cell is absent
  key <- d[, .(id, t = t - k, val = get(v))]
  key[d, on = .(id, t), val]
}
out <- list()
emit <- function(ds, spec, h, m, V) {
  b <- coef(m); se <- sqrt(diag(V))
  for (v in setdiff(names(b), "(Intercept)")) out[[length(out) + 1]] <<- data.table(
    dataset = ds, spec = spec, h = h, term = v, b = b[[v]], se = se[[v]], nobs = nobs(m))
}
for (ds in c("xtarsim_s8y0", "xtarsim_s0y8", "xtarsim_s8y8", "fx_cdols", "fx_unbalanced")) {
  dir <- if (startsWith(ds, "fx_")) "synthetic" else "stata"
  d <- fread(file.path(root, "test", "fixtures", dir, paste0(ds, ".csv")))
  if (ds == "fx_unbalanced") d[, s := x]   # the unit-level shock x is the CD-OLS shock
  setkey(d, id, t)
  d[, Ly := shiftcol(d, "y", -1)]
  for (h in 0:9) {
    d[, ld := shiftcol(d, "y", h) - Ly]
    m <- feols(ld ~ s + l(s, 1), d, panel.id = ~id + t)
    emit(ds, "fixest_cdols_cl_id", h, m, vcov(m, cluster = ~id))
    emit(ds, "fixest_cdols_hetero", h, m, vcov(m, vcov = "hetero"))
    emit(ds, "fixest_cdols_cl_id_none", h, m,
         vcov(m, cluster = ~id, ssc = ssc(K.adj = FALSE, G.adj = FALSE)))
  }
}
res <- rbindlist(out)
for (v in intersect(c("b", "se"), names(res))) set(res, j = v, value = sprintf("%.17g", res[[v]]))
fwrite(res, file.path(root, "test", "fixtures", "r", "p3_fixest.csv"))
cat(sprintf("p3_fixest: %d rows (fixest %s)\n", nrow(res), as.character(packageVersion("fixest"))))
