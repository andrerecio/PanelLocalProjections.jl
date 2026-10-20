# Every struct and type alias of the package. Constructors and methods live in the files
# named in the comments.

# ---------------------------------------------------------------------------------------
# Panel index (panelindex.jl)
# ---------------------------------------------------------------------------------------

"""
    PanelIndex(id, time; frequency = nothing)

Index of a panel: unit codes (`1:N`, in sorted order of the id values), integer calendar
periods, and a dense lookup from (unit, period) to the data row.

- Integer `time` values are periods: consecutive integers are consecutive periods
  (e.g. years, or Stata `%tq` quarters).
- `Date`/`DateTime` values need `frequency` ∈ `(:year, :halfyear, :quarter, :month)`.
- Float values must be integral.

Missing ids or periods and duplicate (unit, period) pairs are errors. Gaps are allowed:
leads and lags into an absent cell are undefined (`NaN`).
"""
struct PanelIndex
    unit::Vector{Int32}          # unit code of each data row, 1:N
    period::Vector{Int}          # calendar period of each data row
    unitlabels::Vector{String}   # id value of each unit code
    pmin::Int                    # first and last period present
    pmax::Int
    frequency::Symbol            # :none for integer periods
    rowof::Matrix{Int32}         # rowof[i, p - pmin + 1] = data row, 0 if absent
end

# ---------------------------------------------------------------------------------------
# Formula terms (formula.jl). `@formula` yields StatsModels terms; `_rewrite` turns the
# package's helpers into these structs. They are evaluated by `_cols` on the panel, never
# by StatsModels.modelcols, so every shift goes through the PanelIndex. The user-facing
# ones are StatsModels terms, so that `~`, `+` and `&` compose them in formulas built
# programmatically.
# ---------------------------------------------------------------------------------------

"Left-hand side of a panel local projection: the response at horizon `h`."
abstract type PanelResponse <: StatsModels.AbstractTerm end

"`leads(y)`: y_{i,t+h}."
struct LeadsResponse{T} <: PanelResponse
    term::T
end

"`ldiff(y)`: y_{i,t+h} − y_{i,t−1}."
struct LDiffResponse{T} <: PanelResponse
    term::T
end

"`anchor(y, z)`: y_{i,t+h} − z_{i,t}."
struct AnchorResponse{T, S} <: PanelResponse
    term::T
    base::S
end

"`cumul(y)`: Σ_{k=0}^{h} y_{i,t+k}."
struct CumulResponse{T} <: PanelResponse
    term::T
end

"A plain response `y` (a pre-built outcome), allowed only with `horizon = 0`."
struct PlainResponse{T} <: PanelResponse
    term::T
end

"`lags(x, p)` or `lags(x, a:b)`: x_{i,t−ℓ} for each ℓ in the range."
struct LagsTerm{T} <: StatsModels.AbstractTerm
    term::T
    lags::UnitRange{Int}
end

"`firstdiff(x)`: x_{i,t} − x_{i,t−1}."
struct FirstDiffTerm{T} <: StatsModels.AbstractTerm
    term::T
end

"`fe(v)` or `fe(a) & fe(b)`: an absorbed fixed effect (an interaction for several vars)."
struct FETerm <: StatsModels.AbstractTerm
    vars::Vector{Symbol}
end

# Internal term tree for regressors and helper arguments (from StatsModels terms).
struct VarTerm
    sym::Symbol
end
struct ConstTerm
    value::Float64
end
struct FunTerm
    f::Function
    args::Vector{Any}
    name::String
end
struct ProductTerm
    terms::Vector{Any}
end

# Data needed to evaluate terms on the panel.
struct _EvalContext
    pidx::PanelIndex
    columns::Dict{Symbol, Vector{Float64}}
end

# ---------------------------------------------------------------------------------------
# Estimators (estimators.jl)
# ---------------------------------------------------------------------------------------

"Estimator of a panel local projection; dispatched on in `_fit_horizon`."
abstract type PanelLPEstimator end

"""
    WithinFE()

Fixed-effects (within/LSDV) estimator: OLS after absorbing the `fe()` terms of the
formula; pooled OLS with an intercept when there are none.
"""
struct WithinFE <: PanelLPEstimator end

"""
    SplitPanelJackknife(; split = :ceil)

Split-panel jackknife correction of the fixed-effects estimator (Mei, Sheng & Shi 2026,
eq. 16): θ̂ = 2θ̂_full − ½(θ̂_a + θ̂_b). The halves are formed per unit on each horizon's
estimation sample, from its rows in time order:

- `split = :ceil` puts the first ⌈n_i/2⌉ rows in half a, as the authors' R code and
  `xtlp` do;
- `split = :floor` puts the first ⌊n_i/2⌋ there, the paper's 𝒯_a = {t ≤ T_h/2}.

Every fixed effect, time effects included, is re-estimated inside each half. The variance
uses the full-sample bread, the SPJ residuals and the scores d = 2X̃ − X̃_half, which
reduce to the paper's d* on balanced panels with an even number of periods. Requires a
unit fixed effect, `fe(id)`. A horizon where a half is rank-deficient (a regressor
absorbed inside it) is `:spj_half_singular`, with NaN coefficients.
"""
struct SplitPanelJackknife <: PanelLPEstimator
    split::Symbol
    SplitPanelJackknife(split::Symbol) = new(_check_option(split, (:ceil, :floor), "split"))
end
SplitPanelJackknife(; split::Symbol = :ceil) = SplitPanelJackknife(split)

# ---------------------------------------------------------------------------------------
# Fixed-effect absorption (absorb.jl)
# ---------------------------------------------------------------------------------------

"How fixed effects are absorbed (within-transformed)."
abstract type AbsorptionMethod end

# Fixed-effect codes on one estimation sample, re-indexed to 1:L per dimension.
struct _FESample
    codes::Vector{Vector{Int32}}
    counts::Vector{Vector{Int}}
end

"""
    AlternatingProjections(; tol = 1e-12, maxiter = 10_000, accelerate = true)

Within transformation by alternating projections (Gauss–Seidel sweeps over the fixed-effect
dimensions), with an Irons–Tuck extrapolation after every two sweeps when `accelerate`. A
single dimension is demeaned exactly in one pass. Iteration stops when a sweep moves no
element by more than `tol · max|v₁|`, where `v₁` is the column after its first sweep
(which has removed its level): the rule does not depend on the units or the level of the
data. On slowly mixing panels the remaining error can exceed the last step by a factor
1/(1 − ρ), with ρ the convergence rate of the sweeps.
"""
struct AlternatingProjections <: AbsorptionMethod
    tol::Float64
    maxiter::Int
    accelerate::Bool
end
function AlternatingProjections(; tol::Real = 1e-12, maxiter::Integer = 10_000,
        accelerate::Bool = true)
    tol > 0 || throw(ArgumentError("tol must be positive"))
    maxiter >= 1 || throw(ArgumentError("maxiter must be at least 1"))
    return AlternatingProjections(Float64(tol), Int(maxiter), accelerate)
end

"""
    OneSweep()

One demeaning pass per fixed-effect dimension, in formula order. With `fe(id) + fe(t)` this
is the unit-then-time sweep of the Mei–Sheng–Shi replication code (`LP_panel_all.r`),
which is the exact two-way within transformation only on balanced samples. A
compatibility switch for reproducing published numbers; prefer `AlternatingProjections`.
"""
struct OneSweep <: AbsorptionMethod end

# ---------------------------------------------------------------------------------------
# Design (design.jl)
# ---------------------------------------------------------------------------------------

# Everything `panellp` needs that does not depend on the horizon, over all data rows.
struct PanelDesign{R <: PanelResponse}
    pidx::PanelIndex
    response::R
    responsename::String
    columns::Dict{Symbol, Vector{Float64}}   # numeric formula variables, missing → NaN
    X::Matrix{Float64}                       # right-hand side (NaN where undefined)
    coefnames::Vector{String}
    shockcols::Vector{Int}
    intercept::Bool
    rhsok::BitVector                         # rows where every RHS column is finite
    fe::Vector{Vector{Int32}}                # codes of each absorbed FE dimension (0 = missing)
    fenames::Vector{String}
    groups::Vector{Vector{Int32}}            # codes of the score-sum dimensions (0 = missing)
    groupnames::Vector{Symbol}               # :unit, :time, FE dims, declared clusters
    ngroups::Vector{Int}                     # global number of groups per dimension
    subset::BitVector
end

# ---------------------------------------------------------------------------------------
# Per-horizon fit (estimators.jl, scores.jl)
# ---------------------------------------------------------------------------------------

# Degrees-of-freedom bookkeeping of one horizon's fit.
struct DofInfo
    nobs::Int
    k::Int                    # estimated coefficients (intercept included when present)
    df_a::Int                 # rank of the absorbed fixed effects
    felevels::Vector{Int}     # levels of each FE dimension in the sample
    fecomponents::Matrix{Int} # [i, j], i < j: connected components of FE dims i and j
    ngroups::Vector{Int}      # groups present in the sample, per score-sum dimension
    nested::Matrix{Bool}      # nested[f, g]: FE dimension f nested in score dimension g
    singletons::Int           # rows dropped as singletons
end

# Sufficient statistics of the scores u_r = d_r ê_r of one horizon: every variance in
# vcov.jl is a function of these, so no design matrix is kept after estimation.
struct ScoreStats
    Mobs::Matrix{Float64}             # Σ_r u_r u_r'
    S::Vector{Matrix{Float64}}        # per score dimension: score sums by global group code
    pairs::Vector{Tuple{Int, Int}}    # pairs (a, b) of score dimensions, a < b
    Spair::Vector{Matrix{Float64}}    # score sums by intersection group (when :own)
    pairkind::Vector{Symbol}          # :obs (use Mobs), :first (= a), :second (= b), :own
end

struct HorizonFit
    h::Int
    # :ok; :empty (no observations); :singular (a shock is absorbed by the fixed effects or
    # collinear in this horizon's sample); :spj_half_singular (a split-panel half is
    # rank-deficient). Coefficients are NaN unless :ok.
    status::Symbol
    coef::Vector{Float64}     # NaN for omitted coefficients
    keep::BitVector           # coefficient estimated (not omitted)
    absorbed::BitVector       # omitted because the fixed effects (or the constant) absorb it
    bread::Matrix{Float64}    # (D'D)⁻¹ on the kept coefficients, NaN elsewhere
    # sd (n − 1) of the raw regressors in the sample, floored at 1e-3; 1 for the intercept
    xsd::Vector{Float64}
    esample::BitVector        # estimation sample over data rows
    timepresent::BitVector    # calendar periods (pmin:pmax) present in the sample
    dof::DofInfo
    scores::ScoreStats
    iterations::Int           # absorption sweeps (0 without FE)
    converged::Bool
end

struct PanelSummary
    id::Symbol
    time::Symbol
    frequency::Symbol
    nunits::Int
    nperiods::Int             # span of calendar periods, pmax - pmin + 1
    ntimes::Int               # distinct periods present in the data
    nrows::Int
    pmin::Int
    pmax::Int
end

"""
    PanelLocalProjection

Result of [`panellp`](@ref): one fit per horizon. Variances are computed afterwards with
`vcov(spec, r)`; see [`Robust`](@ref), [`Clustered`](@ref), [`DriscollKraaySE`](@ref).
"""
struct PanelLocalProjection{E <: PanelLPEstimator}
    estimator::E
    absorb::AbsorptionMethod
    formula::StatsModels.FormulaTerm
    response::String
    coefnames::Vector{String}
    intercept::Bool                  # the first coefficient is the intercept
    shocks::Vector{String}
    horizons::UnitRange{Int}
    fits::Vector{HorizonFit}
    groupnames::Vector{Symbol}
    ngroups::Vector{Int}
    fenames::Vector{String}
    panel::PanelSummary
    sample::Symbol
    singletons::Symbol
end

# ---------------------------------------------------------------------------------------
# Variance specifications and results (vcov.jl)
# ---------------------------------------------------------------------------------------

"Specification of a variance estimator, resolved per horizon by `vcov(spec, r)`."
abstract type PanelVcov end

"""
    Robust(; ssc = :stata)

Heteroskedasticity-robust (Eicker–Huber–White) variance. `ssc = :none` is HC0; `:stata`
is HC1 with the factor `n/(n − K − df_a)` of `regress, vce(robust)`, `reghdfe` and `xtlp
vce(robust)`, and t(n − K − df_a) critical values.
"""
struct Robust <: PanelVcov
    ssc::Symbol
    Robust(ssc::Symbol) = new(_check_option(ssc, (:none, :stata), "ssc"))
end

"""
    Clustered(dims...; ssc = :stata, psd = length(dims) == 2)

Cluster-robust variance, one- or two-way (Cameron–Gelbach–Miller). Each dimension is
`:unit`, `:time`, a fixed-effect or declared cluster column, or a tuple of columns for
their intersection (e.g. `(:sic3, :quarter)`).

`ssc` sets the small-sample correction:
- `:none`: no factor; normal critical values.
- `:stata`: `(n−1)/(n − a − K − df_a_eff) · G/(G−1)` with `G` the smallest number of
  clusters, absorbed effects nested in a cluster dimension not counted (`a = 1` if any
  is); t(G − 1) critical values. As `xtlp`, `reghdfe`, `xtreg, fe vce(cluster)`.
- `:mss`: `N/(N−1) · (n−1)/(n−K)` for `Clustered(:unit)`, or `G/(G−1) · (n−1)/(n−N−K)`
  with `G = min(N, T)` for `Clustered(:unit, :time)`, with N and T the numbers of units
  and periods in the data, as in the Mei–Sheng–Shi replication code with `robust = T`;
  normal critical values. Only these two specifications exist in that code.

`psd` repairs a two-way variance that is not positive semi-definite by clipping its
negative eigenvalues to zero:
- `:stata` (or `true`, the default for two dimensions) clips in the metric scaled by the
  standard deviations of the raw regressors (floored at 1e-3), as `reghdfe` and `xtlp`
  do. With an intercept, the slopes are clipped on their own and the intercept enters
  the metric with scale 1, as in `reghdfe`.
- `:raw` clips V itself, as the Mei–Sheng–Shi code with `eigens = T` and `fixest`'s
  default `vcov_fix = TRUE` (which floors the eigenvalues at 1e-16 instead of 0).
- `:none` (or `false`) leaves V unchanged; a negative variance then has a NaN standard
  error.
"""
struct Clustered <: PanelVcov
    dims::Vector{Symbol}
    ssc::Symbol
    psd::Symbol
    function Clustered(dims::Vector{Symbol}, ssc::Symbol, psd::Symbol)
        1 <= length(dims) <= 2 ||
            throw(ArgumentError("Clustered takes one or two dimensions; got $(length(dims))"))
        length(dims) == 2 && dims[1] == dims[2] &&
            throw(ArgumentError("the two cluster dimensions are the same"))
        return new(dims, _check_option(ssc, (:none, :stata, :mss), "ssc"),
            _check_option(psd, (:none, :stata, :raw), "psd"))
    end
end

"""
    DriscollKraaySE(bw = :nw94; ssc = :stata, timelags = :calendar)

Driscoll–Kraay variance: Bartlett-weighted long-run variance of the cross-sectional sums
of the scores, ordered by calendar period.

`bw` sets the number of lags m (weights `1 − j/(m+1)`):
- `:nw94`: `⌊4(T/100)^{2/9}⌋` (`xtscc`, `xtlp`).
- `:as2026`: `⌈4(T/100)^{2/9}⌉ − 1`, Almuzara–Sancibrián's kernel bandwidth M = ⌈⋯⌉
  with weights `1 − j/M`.
- `:mss2026`: `⌊T^{1/4}⌋` (Mei–Sheng–Shi, appendix B.2.2).
- an integer.

`ssc`:
- `:none`: no factor; normal critical values.
- `:stata`: `(n−1)/(n − K − df_a^{DK}) · T/(T−1)`, as `xtlp`'s default; t(T − 1).
- `:xtscc`: `(n−1)/(n − K − 1{FE}) · T/(T−1)`, as `xtscc` and `xtlp …, nodfadj`; t(T − 1).

`timelags = :position` lags over the periods present in the sample (as `xtlp`/`xtscc`)
instead of calendar periods.
"""
struct DriscollKraaySE <: PanelVcov
    rule::Symbol      # :nw94, :as2026, :mss2026 or :fixed
    lags::Int         # used when rule == :fixed
    ssc::Symbol
    timelags::Symbol
    function DriscollKraaySE(rule::Symbol, lags::Int, ssc::Symbol, timelags::Symbol)
        _check_option(rule, (:nw94, :as2026, :mss2026, :fixed), "bw")
        lags >= 0 || throw(ArgumentError("the number of lags must be non-negative"))
        return new(rule, lags, _check_option(ssc, (:none, :stata, :xtscc), "ssc"),
            _check_option(timelags, (:calendar, :position), "timelags"))
    end
end

# Per-horizon details of a variance computation.
struct VcovInfo
    factor::Float64       # small-sample factor applied
    groups::Int           # smallest number of clusters (0 if not clustered)
    periods::Int          # time periods entering a time aggregation (0 otherwise)
    lags::Int             # Driscoll–Kraay lags (0 otherwise)
    psdfixed::Bool        # negative eigenvalues were clipped
end

"""
    PanelLPVcov

Variances of a [`PanelLocalProjection`](@ref) for one specification: a matrix per horizon
(NaN for omitted coefficients, and where the variance is undefined), the degrees of
freedom of the critical values per coefficient (`Inf` = normal), and per-horizon details.
"""
struct PanelLPVcov{S <: PanelVcov}
    spec::S
    coefnames::Vector{String}
    horizons::UnitRange{Int}
    V::Vector{Matrix{Float64}}
    df::Vector{Vector{Float64}}
    info::Vector{VcovInfo}
    fingerprint::UInt64       # of the result it was computed from (checked by summarize)
end

# ---------------------------------------------------------------------------------------
# Summaries (summary.jl)
# ---------------------------------------------------------------------------------------

"""
    IRFSummary

Impulse response of one coefficient across horizons, with standard errors and confidence
bounds; a Tables.jl table (`DataFrame(s)` works).
"""
struct IRFSummary
    term::String
    response::String
    spec::String
    level::Float64
    scale::Float64
    horizon::Vector{Int}
    coef::Vector{Float64}
    se::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    df::Vector{Float64}
    nobs::Vector{Int}
end
