using TestItemRunner

const ROOT = dirname(@__DIR__)

# Local-only fixture sets (see test/fixtures/README.md). Items tagged with a set are
# skipped, with a note, when its sentinel file is absent (e.g. on CI).
const LOCAL_SETS = Dict(
    :local_mss => joinpath(ROOT, "docs", "references", "code", "panel-lp-replication",
    "applications", "empirical_RR_f4_lngdp_1980.csv"))
const MISSING_SETS = [k for (k, f) in LOCAL_SETS if !isfile(f)]
isempty(MISSING_SETS) ||
    @info "Local fixture sets absent; their test items are skipped" MISSING_SETS

# Paper-scale Monte Carlo runs only on request: PLP_MC_FULL=1.
const MC_FULL = get(ENV, "PLP_MC_FULL", "0") == "1"

@run_package_tests filter = ti -> !any(in(ti.tags), MISSING_SETS) &&
                                  (MC_FULL || !(:mc_full in ti.tags))
