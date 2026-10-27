# The four Mei–Sheng–Shi (2026) empirical applications, estimated exactly as in their
# replication notebook (applications/replication.ipynb): one LP_panel call per horizon on
# the pre-built outcome columns, FE and SPJ, with each application's options.
#     Rscript crosscheck/R/mss_apps_local.R
# Inputs are the local, gitignored replication data under docs/references/, so the output
# (test/fixtures/local/mss_apps_R.csv) is local too. Values are unscaled (RR is plotted ×7).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
app <- file.path(root, "docs", "references", "code", "panel-lp-replication", "applications")
if (!dir.exists(app)) { cat("MSS replication data absent (local-only): skipped\n"); quit(status = 0) }
e <- new.env(); sys.source(file.path(app, "LP_panel_all.r"), envir = e)

apps <- list(
  RR_gdp = list(file = "empirical_RR_f4_lngdp_1980.csv", X = "CRISIS",
    C = c(paste0("l", 1:4, "LNGDP"), paste0("l", 1:4, "CRISIS")), Y = paste0("f", 0:10, "LNGDP"),
    h = 0:10, te = TRUE, tw = FALSE, robust = FALSE, eig = FALSE, ymult = 1),
  RR_unemp = list(file = "empirical_RR_f4_lnunemp_1980.csv", X = "CRISIS",
    C = c(paste0("l", 1:4, "UNEMP"), paste0("l", 1:4, "CRISIS")), Y = paste0("f", 0:10, "UNEMP"),
    h = 0:10, te = TRUE, tw = FALSE, robust = FALSE, eig = FALSE, ymult = 1),
  BVX = list(file = "empirical_BVX_t1_y.csv", X = c("R_B", "R_N"),
    C = c(paste0("L", 1:3, "R_B"), paste0("L", 1:3, "R_N"), "D1y", paste0("L", 1:3, "D1y"),
          "D1d_y", paste0("L", 1:3, "D1d_y")), Y = paste0("Fd", 1:6, "y"),
    h = 1:6, te = FALSE, tw = TRUE, robust = TRUE, eig = TRUE, ymult = 100),
  MSV = list(file = "empirical_MSV_f2.csv", X = c("L0HHD_L1GDP", "L0NFD_L1GDP"),
    C = c(paste0("L", 0:4, "y"), paste0("L", 1:4, "HHD_L1GDP"), paste0("L", 1:4, "NFD_L1GDP")),
    Y = paste0("F", 1:10, "y"), h = 1:10, te = FALSE, tw = TRUE, robust = FALSE, eig = FALSE, ymult = 1),
  CS = list(file = "empirical_CS_f3.csv", X = "CRISIS",
    C = c(paste0("l", 1:4, "CRISIS"), paste0("l", 1:4, "GRRT_WB")), Y = paste0("cf", 1:10, "GRRT_WB"),
    h = 1:10, te = FALSE, tw = FALSE, robust = FALSE, eig = FALSE, ymult = 1))

out <- list()
for (nm in names(apps)) {
  a <- apps[[nm]]
  data <- read.csv(file.path(app, a$file))
  for (k in seq_along(a$Y)) {
    Y <- a$Y[k]
    data[, Y] <- a$ymult * data[, Y]
    for (method in c("FE", "SPJ")) {
      fit <- e$LP_panel(data, Y.name = Y, X.name = a$X, c.name = a$C, method = method,
                        te = a$te, two_way_cluster = a$tw, robust = a$robust, eigens = a$eig)
      for (j in seq_along(a$X)) out[[length(out) + 1]] <- data.table(
        app = nm, method = method, h = a$h[k], term = a$X[j],
        b = sprintf("%.17g", fit$IRF[1, j]), se = sprintf("%.17g", fit$se[1, j]))
    }
  }
}
res <- rbindlist(out)
dir.create(file.path(root, "test", "fixtures", "local"), showWarnings = FALSE)
fwrite(res, file.path(root, "test", "fixtures", "local", "mss_apps_R.csv"))
cat(sprintf("mss_apps_local: %d rows\n", nrow(res)))
