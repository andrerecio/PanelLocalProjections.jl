# Phase 1 cross-check with fixest: FE panel LP coefficients and variances.
#     Rscript crosscheck/R/p1_fixest.R
# Writes test/fixtures/r/p1_fixest.csv; spec labels map to panellp/vcov calls in
# test/test_crosscheck_p1.jl.
suppressPackageStartupMessages({ library(fixest); library(data.table) })
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
setFixest_nthreads(1)

out <- list()
emit <- function(ds, spec, h, m, V) {
  b <- coef(m); se <- suppressWarnings(sqrt(diag(V)))   # NaN for a negative variance
  for (v in names(b)) out[[length(out) + 1]] <<- data.table(
    dataset = ds, spec = spec, h = h, term = v, b = b[[v]], se = se[[v]], nobs = nobs(m))
}
none <- ssc(K.adj = FALSE, G.adj = FALSE)
f2 <- function(h) as.formula(sprintf("f(y, %d) ~ x + l(x, 1) + l(y, 1) + w | id + t", h))
f1 <- function(h) as.formula(sprintf("f(y, %d) ~ x + l(x, 1) + l(y, 1) + w | id", h))
for (ds in c("fx_balanced", "fx_unbalanced")) {
  d <- fread(file.path(root, "test", "fixtures", "synthetic", paste0(ds, ".csv")))
  for (h in 0:6) {
    for (rm in c("singleton", "none")) {
      tag <- if (rm == "none") "_keepsing" else ""
      m1 <- feols(f1(h), d, panel.id = ~id + t, fixef.rm = rm, fixef.tol = 1e-11)
      m2 <- feols(f2(h), d, panel.id = ~id + t, fixef.rm = rm, fixef.tol = 1e-11)
      emit(ds, paste0("fixest_fe_cl_id_none", tag), h, m1, vcov(m1, cluster = ~id, ssc = none))
      emit(ds, paste0("fixest_fe_hetero_none", tag), h, m1, vcov(m1, vcov = "hetero", ssc = none))
      emit(ds, paste0("fixest_tfe_cl_id_none", tag), h, m2, vcov(m2, cluster = ~id, ssc = none))
      emit(ds, paste0("fixest_tfe_twoway_none", tag), h, m2, vcov(m2, cluster = ~id + t, ssc = none))
      emit(ds, paste0("fixest_tfe_cl_id_default", tag), h, m2, vcov(m2, cluster = ~id))
    }
  }
}
# Four units of fx_balanced: the two-way variance is not PSD, so fixest's default repair
# (vcov_fix = TRUE, an eigenvalue clip of V itself) fires; vcov_fix = FALSE leaves it.
d <- fread(file.path(root, "test", "fixtures", "synthetic", "fx_balanced.csv"))[id <= 4]
for (h in 0:6) {
  m <- feols(f1(h), d, panel.id = ~id + t, fixef.rm = "singleton", fixef.tol = 1e-11)
  emit("fx_balanced_id4", "fixest_fe_twoway_none_id4", h, m, vcov(m, cluster = ~id + t, ssc = none))
  emit("fx_balanced_id4", "fixest_fe_twoway_none_nofix_id4", h, m,
       vcov(m, cluster = ~id + t, ssc = none, vcov_fix = FALSE))
}
res <- rbindlist(out)
for (v in intersect(c("b", "se"), names(res))) set(res, j = v, value = sprintf("%.17g", res[[v]]))
fwrite(res, file.path(root, "test", "fixtures", "r", "p1_fixest.csv"))
cat(sprintf("p1_fixest: %d rows (fixest %s)\n", nrow(res), as.character(packageVersion("fixest"))))
