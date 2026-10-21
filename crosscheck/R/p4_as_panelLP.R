# Phase 4 cross-check with the Almuzara–Sancibrián reference code, panel_LP.R
# (github.com/TinchoAlmuzara/PanelLocalProjections, R_code/panel_LP.R; local copy in
# docs/references/code, which is not committed: the code has no licence). Their function is
# only *run* here, on our own synthetic panels; its outputs are committed as
# test/fixtures/r/p4_as_panelLP.csv (see test/fixtures/README.md).
#     Rscript crosscheck/R/p4_as_panelLP.R
# Spec labels map to panellp/vcov calls in test/test_crosscheck_p4.jl.
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args[grep("^--file=", args)])), "..", ".."))
src <- file.path(root, "docs", "references", "code", "PanelLocalProjections", "R_code", "panel_LP.R")
if (!file.exists(src)) stop("the Almuzara–Sancibrián code is not at ", src)
# panel_LP.R installs its packages when they are missing: refuse instead
for (p in c("tidyverse", "pracma", "lubridate", "fixest"))
  if (!requireNamespace(p, quietly = TRUE)) stop("R package ", p, " is not installed")
suppressPackageStartupMessages(suppressWarnings(source(src)))
fixest::setFixest_nthreads(1)

H <- 6
out <- list()
for (panel in c("fx_tlahr", "fx_tlahr_unbal")) {
  d <- read.csv(file.path(root, "test", "fixtures", "synthetic", paste0(panel, ".csv")))
  # numeric ids and periods (character ids make panel_LP.R zero-fill missing rows)
  id <- as.numeric(d$id); t <- as.numeric(d$t)
  specs <- list(
    s_twfe   = list(s = d$s, W = NULL, FE = cbind(id, t)),
    ones_unit = list(s = NULL, W = NULL, FE = cbind(id)),
    sv_twfe  = list(s = d$sv, W = NULL, FE = cbind(id, t)),
    s_twfe_w = list(s = d$s, W = cbind(d$sv), FE = cbind(id, t)),
    s2_twfe  = list(s = cbind(d$s, d$sv), W = NULL, FE = cbind(id, t)))
  for (sp in names(specs)) for (pmax in c(0, 2)) for (cumul in c(FALSE, TRUE))
    for (small in c(FALSE, TRUE)) {
      a <- specs[[sp]]
      r <- suppressWarnings(panel_LP(d$y, d$X, s = a$s, i_index = id, t_index = t,
        W = a$W, FE = a$FE, H = H, p_max = pmax, small_sample = small, cumulative = cumul))
      for (k in seq_len(ncol(r$estimate))) out[[length(out) + 1]] <- data.frame(
        panel = panel, spec = sp, pmax = pmax, small = as.integer(small),
        cumulative = as.integer(cumul), h = 0:H, k = k,
        estimate = sprintf("%.17g", r$estimate[, k]), se = sprintf("%.17g", r$SE[, k]),
        df = sprintf("%.17g", r$df[, k]), pval = sprintf("%.17g", r$pval[, k]))
    }
}
res <- do.call(rbind, out)
write.csv(res, file.path(root, "test", "fixtures", "r", "p4_as_panelLP.csv"), row.names = FALSE,
  quote = FALSE)
cat(sprintf("p4_as_panelLP: %d rows (R %s, fixest %s, pracma %s)\n", nrow(res),
  getRversion(), packageVersion("fixest"), packageVersion("pracma")))
