# Smoke check of the R side of the cross-checks (Phase 0 gate).
#     Rscript crosscheck/R/smoke.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script), "..", ".."))

cat(R.version.string, "\n")
for (p in c("fixest", "sandwich", "data.table", "pracma", "haven", "lubridate", "tidyverse")) {
  cat(sprintf("%-11s %s\n", p, as.character(packageVersion(p))))
}

refs <- file.path(root, "docs", "references", "code")
mss <- file.path(refs, "panel-lp-replication", "applications", "LP_panel_all.r")
as_r <- file.path(refs, "PanelLocalProjections", "R_code", "panel_LP.R")
if (file.exists(mss)) {
  e <- new.env(); sys.source(mss, envir = e)
  stopifnot(is.function(e$LP_panel)); cat("MSS LP_panel: ok\n")
} else cat("MSS reference code absent (local-only): skipped\n")
if (file.exists(as_r)) {
  # panel_LP.R calls require()/install.packages() at load time; every package it needs
  # is checked above, so nothing gets installed.
  e <- new.env(); suppressPackageStartupMessages(sys.source(as_r, envir = e))
  stopifnot(is.function(e$panel_LP)); cat("A&S panel_LP: ok\n")
} else cat("A&S reference code absent (local-only): skipped\n")
cat("smoke: ok\n")
