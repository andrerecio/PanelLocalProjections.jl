# Phases 1–2 cross-check with the Mei–Sheng–Shi replication function LP_panel
# (docs/references/code/panel-lp-replication/applications/LP_panel_all.r, MIT; local only).
#     Rscript crosscheck/R/p1p2_mss_lppanel.R
# LP_panel is run as in the replication notebook: one call per horizon (H = 0) on leads
# and lags pre-built by calendar joins on (id, t ± k), with controls passed explicitly.
# Writes test/fixtures/r/p1p2_mss_lppanel.csv (shock coefficient and SE).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
src <- file.path(root, "docs", "references", "code", "panel-lp-replication", "applications", "LP_panel_all.r")
if (!file.exists(src)) { cat("MSS reference code absent (local-only): skipped\n"); quit(status = 0) }
e <- new.env(); sys.source(src, envir = e)

shiftcol <- function(d, v, k) {   # v at (id, t + k), NA where the cell is absent
  key <- d[, .(id, t = t - k, val = get(v))]
  key[d, on = .(id, t), val]
}
specs <- list(
  mss_fe_unit            = list(method = "FE",  te = FALSE, tw = FALSE, robust = FALSE, eig = FALSE),
  mss_fe_unit_robust     = list(method = "FE",  te = FALSE, tw = FALSE, robust = TRUE,  eig = FALSE),
  mss_te_unit            = list(method = "FE",  te = TRUE,  tw = FALSE, robust = FALSE, eig = FALSE),
  mss_te_unit_robust     = list(method = "FE",  te = TRUE,  tw = FALSE, robust = TRUE,  eig = FALSE),
  mss_fe_tw              = list(method = "FE",  te = FALSE, tw = TRUE,  robust = FALSE, eig = FALSE),
  mss_fe_tw_robust_eig   = list(method = "FE",  te = FALSE, tw = TRUE,  robust = TRUE,  eig = TRUE),
  mss_te_tw              = list(method = "FE",  te = TRUE,  tw = TRUE,  robust = FALSE, eig = FALSE),
  mss_spj_unit           = list(method = "SPJ", te = FALSE, tw = FALSE, robust = FALSE, eig = FALSE),
  mss_spj_unit_robust    = list(method = "SPJ", te = FALSE, tw = FALSE, robust = TRUE,  eig = FALSE),
  mss_spj_te_unit        = list(method = "SPJ", te = TRUE,  tw = FALSE, robust = FALSE, eig = FALSE),
  mss_spj_tw             = list(method = "SPJ", te = FALSE, tw = TRUE,  robust = FALSE, eig = FALSE),
  mss_spj_tw_robust_eig  = list(method = "SPJ", te = FALSE, tw = TRUE,  robust = TRUE,  eig = TRUE),
  mss_spj_te_unit_robust = list(method = "SPJ", te = TRUE,  tw = FALSE, robust = TRUE,  eig = FALSE),
  mss_spj_te_tw          = list(method = "SPJ", te = TRUE,  tw = TRUE,  robust = FALSE, eig = FALSE),
  mss_spj_te_tw_robust_eig = list(method = "SPJ", te = TRUE, tw = TRUE, robust = TRUE, eig = TRUE))
out <- list()
# fx_balanced_id4 (units 1-4 of fx_balanced): the two-way variances are not PSD there, so
# eigens = T repairs them
for (ds in c("fx_balanced", "fx_unbalanced", "fx_balanced_id4")) {
  d <- fread(file.path(root, "test", "fixtures", "synthetic", paste0(sub("_id4$", "", ds), ".csv")))
  if (ds == "fx_balanced_id4") d <- d[id <= 4]
  setkey(d, id, t)
  d[, Lx := shiftcol(d, "x", -1)]
  d[, Ly := shiftcol(d, "y", -1)]
  for (h in 0:3) {
    d[, Fy := shiftcol(d, "y", h)]
    df <- as.data.frame(d[, .(id, t, Fy, x, Lx, Ly, w)])
    for (nm in names(specs)) {
      s <- specs[[nm]]
      r <- e$LP_panel(df, Y.name = "Fy", X.name = "x", c.name = c("Lx", "Ly", "w"),
                      id.name = "id", time.name = "t", method = s$method, te = s$te,
                      two_way_cluster = s$tw, robust = s$robust, eigens = s$eig)
      out[[length(out) + 1]] <- data.table(dataset = ds, spec = nm, h = h, term = "x",
                                           b = r$IRF[1, 1], se = r$se[1, 1])
    }
  }
}
res <- rbindlist(out)
for (v in intersect(c("b", "se"), names(res))) set(res, j = v, value = sprintf("%.17g", res[[v]]))
fwrite(res, file.path(root, "test", "fixtures", "r", "p1p2_mss_lppanel.csv"))
cat(sprintf("p1p2_mss_lppanel: %d rows\n", nrow(res)))
