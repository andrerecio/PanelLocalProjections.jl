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

"""
    CDOLS(; validate = true)

Cumulative-difference OLS (Ugarte-Ruiz 2026): pooled OLS, with an intercept (unless
`0 +`) and without unit effects, of the long difference y_{i,t+h} − y_{i,t−1} on the
shock and its first lag. The long difference removes the unit effect, so there is no
within transformation and no Nickell bias.

With `validate = true` the formula must follow the method:
- the response is `ldiff(y)`, or `cumul(g)` of a growth rate `g`;
- each shock is dated t, right after the base period t − 1 of the long difference (t − ℓ
  with a lagged response, `ldiff(lags(y, ℓ))`), and its first lag is a regressor, however
  it is written (`lags(s, 1)`, `lags(s, 0:1)`, …);
- the outcome enters the right-hand side only lagged beyond the response's start and,
  when the response differences it (`ldiff(y)`, `cumul(firstdiff(y))`), only in
  differences, `lags(firstdiff(y), 1:p)`: a level of y carries the unit effect;
- no fixed effect is nested in the units, in the data or in any horizon's sample (time
  effects are allowed); fixed-effect groups with a single unit give a warning.

The paper also gives every control the lags of the shock (it is not checked).
`validate = false` skips the checks, e.g. for a pre-built long-difference column.
"""
struct CDOLS <: PanelLPEstimator
    validate::Bool
end
CDOLS(; validate::Bool = true) = CDOLS(validate)

"Lag-length rule of [`LagAugmented`](@ref)."
abstract type LagRule end

"""
    ICLags(; criterion = :aic, pmax = nothing, vars = :outcome, sample = :perp)

Almuzara–Sancibrián (2026) rule p(h) = max{1, min(h, p̂)} (eq. 15; supplement B.7), with
p̂ chosen by an information criterion on the *synthetic* time series (fn 11; supplement
B.3–B.4). For each p = 0, …, `pmax`, the outcome is regressed at h = 0 on the formula's
right-hand side and p of its own lags, with the formula's fixed effects; the residuals are
averaged across units with the first shock's exposure weights, ξ̄_t(p), and

    IC(p) = log σ̂²(p) + k(p) · g,   σ̂²(p) = Σ_t ω_t ξ̄_t(p)² / Σ_t ω_t,

with ω_t the synthetic weight of period t (constant on balanced panels with a
time-invariant exposure, where σ̂² is the mean of ξ̄_t²), k(p) the lags candidate p adds to
the regression (none for a lag already among the regressors, or spanned by them), and
g = 2/T (`:aic`) or log(T)/T (`:bic`), T the periods of the p = 0 sample. Ties go to the
smaller p, except between candidates that are the same regression (their lags are already
in the formula): p̂ is then the larger, so the formula's lags count towards it.

- `pmax`: largest candidate; default max{1, ⌈(T − H)^{1/3}⌉}, T the distinct periods in the
  data and H the largest horizon.
- `vars = :outcome_and_shock` also adds p lags of each shock to candidate p (each counted
  in k(p)).
- `sample = :common` fits every candidate on the `pmax`-lag sample.
"""
struct ICLags <: LagRule
    criterion::Symbol
    pmax::Union{Nothing, Int}
    vars::Symbol
    sample::Symbol
    function ICLags(criterion::Symbol, pmax::Union{Nothing, Int}, vars::Symbol, sample::Symbol)
        pmax === nothing || pmax >= 1 ||
            throw(ArgumentError("ICLags: pmax must be at least 1; got $pmax"))
        return new(_check_option(criterion, (:aic, :bic), "criterion"), pmax,
            _check_option(vars, (:outcome, :outcome_and_shock), "vars"),
            _check_option(sample, (:perp, :common), "sample"))
    end
end

"""
    CappedLags(pmax = nothing)

p(h) = min(h, pmax): no lag at h = 0. The rule of the Almuzara–Sancibrián code
(`panel_LP.m`/`panel_LP.R`), whose default `pmax` is ⌈(T − H)^{1/3}⌉ (T the distinct
periods in the data, H the largest horizon). A compatibility switch: the paper recommends
[`ICLags`](@ref).
"""
struct CappedLags <: LagRule
    pmax::Union{Nothing, Int}
    function CappedLags(pmax::Union{Nothing, Integer} = nothing)
        pmax === nothing || pmax >= 0 ||
            throw(ArgumentError("CappedLags: pmax must be non-negative; got $pmax"))
        return new(pmax === nothing ? nothing : Int(pmax))
    end
end

"""
    ConstantLags(p)

The same number of lags p at every horizon ("plain" t-LAHR, p = 1, in Almuzara &
Sancibrián's Table 1).
"""
struct ConstantLags <: LagRule
    p::Int
    function ConstantLags(p::Integer)
        p >= 0 || throw(ArgumentError("ConstantLags: p must be non-negative; got $p"))
        return new(Int(p))
    end
end

"""
    LagAugmented(rule = ICLags())

Lag-augmented fixed-effects estimator of Almuzara & Sancibrián (2026): at horizon h, the
regression adds p(h) lags of the outcome and of every shock column to the formula
(`rule` sets p(h)), so that the scores are serially uncorrelated and inference clustered
by period is valid ([`TLAHR`](@ref)). The lagged outcome is y for `leads(y)` and
`anchor(y, z)`, Δy for `ldiff(y)` and `leads(firstdiff(y))`, and g for `cumul(g)`.
Lags of a shock lag the whole column (s_{i,t−ℓ}X_{t−ℓ} for an exposure that varies over
time). A lag already among the regressors, whatever its name (columns are compared
numerically: `s & L1.X` is the lag of the shock `s & L0.X`), is not added again; one that
the formula's columns span (Δy lags with `lags(y, 1:2)` in the formula) is left out of the
regression. Requires unit fixed effects, `fe(id)`, or finer ones nested in the units (e.g.
`fe(id) & fe(state)`).
"""
struct LagAugmented{R <: LagRule} <: PanelLPEstimator
    rule::R
end
LagAugmented() = LagAugmented(ICLags())

# LagAugmented resolved for one call: the number of lags p(h) at each horizon.
struct _ResolvedLags <: PanelLPEstimator
    lags::Vector{Int}
    horizons::UnitRange{Int}
end

"""
    LagSelection

The lag rule of a [`LagAugmented`](@ref) fit: the rule, its largest lag, the selected
order, and for [`ICLags`](@ref) the criterion, residual variance and number of synthetic
periods for each candidate p = 0, …, `pmax`.
"""
struct LagSelection
    rule::String
    pmax::Int
    p::Int
    ic::Vector{Float64}
    sigma2::Vector{Float64}
    periods::Vector{Int}
end

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
    # per column: the aggregate series it multiplies, per calendar period of the span (NaN
    # where undefined), and its unit-level exposure, per data row; both empty when the
    # column is not an aggregate series times an exposure (`_exposure`)
    aggseries::Vector{Vector{Float64}}
    exposures::Vector{Vector{Float64}}
    # per column: the column whose synthetic regression its t-LAHR refinement uses (0: none)
    owner::Vector{Int}
    timefe::Vector{Int}                      # FE dimensions nested in periods
    lagorder::Vector{Int}                    # per column: 0, or ℓ for an added lag ℓ
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

# Synthetic time series of one coefficient that multiplies an aggregate series (a shock,
# or another such regressor; Almuzara & Sancibrián 2026, Remark 1; supplement B.2): per
# period present in the sample, weighted cross-sectional averages
# Σ_{i∈t} w_i z_it / Σ_{i∈t} w_i² of the within-transformed regressors, response and
# residuals. The weights are its exposure residualised per period (on an intercept when a
# time-type fixed effect is absorbed, and on the exposures of the other regressors of the
# same aggregate series); `Xraw` uses the raw exposure, as the authors' code.
struct SyntheticSeries
    col::Int                      # coefficient index
    weight::Vector{Float64}       # Σ_{i∈t} w_i² (0: the period has no weight)
    X::Matrix{Float64}            # periods × coefficients (NaN for omitted ones)
    y::Vector{Float64}
    resid::Vector{Float64}
    Xraw::Matrix{Float64}
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
    active::BitVector         # columns in this horizon's regression (LagAugmented: p(h) lags)
    nlags::Int                # lags added by LagAugmented (0 otherwise)
    # one per kept coefficient that owns a synthetic regression (`PanelDesign.owner`), rows
    # = periods of `timepresent`
    synthetic::Vector{SyntheticSeries}
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
    lags::Union{Nothing, LagSelection}   # LagAugmented only
    owner::Vector{Int}                   # per coefficient: see PanelDesign
    addedlag::Vector{Int}                # per coefficient: ℓ for a lag added by LagAugmented
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

"""
    TLAHR(; refine = true, weights = :residual)

t-LAHR inference of Almuzara & Sancibrián (2026): the variance clustered by period, with
no small-sample factor, of a lag-augmented regression ([`LagAugmented`](@ref)), whose
scores are serially uncorrelated.

- `refine = false`: CR0 by period (eq. 11), normal critical values. Equal to
  `Clustered(:time; ssc = :none)`; with `WithinFE` and no lags, the paper's "t-HR".
- `refine = true` (the paper's recommendation, supplement B.2): HC2 with the leverage of
  the synthetic time-series regression, and Student-t critical values with the
  Bell–McCaffrey/Imbens–Kolesár degrees of freedom ν of each coefficient. Needs every
  shock to be an aggregate series times a unit-level exposure (`s & X`, or `X` alone).

`weights` sets the synthetic regressions. `:residual` (the paper's) averages with the
exposure residualised per period, on an intercept when a time-type fixed effect is
absorbed (the paper's ŝ_i) and on the exposures of the other regressors of the same
aggregate series, and weights each period by ω_t = Σ_i w_i² (fn 24: with every regressor
an aggregate series times the exposure, the HC2 variance is then CR2 clustered by period).
`:raw` averages with the raw exposure and leaves the periods unweighted, as the authors'
code (`panel_LP.m`/`.R`). They agree on balanced panels with time effects and a
time-invariant exposure. A shock and its lags use the shock's synthetic regression,
another aggregate series times an exposure its own, every other coefficient the first
shock's.
"""
struct TLAHR <: PanelVcov
    refine::Bool
    weights::Symbol
    function TLAHR(refine::Bool, weights::Symbol)
        new(refine,
            _check_option(weights, (:residual, :raw), "weights"))
    end
end

"""
    CosineHAR(B = :llsw)

Equal-weighted cosine (EWC) variance of the scores summed by period (Lazarus, Lewis, Stock
& Watson 2018): Ω = B⁻¹ Σ_{j=1}^{B} Λ_j Λ_j', Λ_j = √2 Σ_t cos(πj(t − ½)/T) S_t, over the
calendar periods from the first to the last of the sample (absent ones contribute zero);
Student-t(B) critical values; no small-sample factor. `B = :llsw` is ⌊0.41 T^{2/3}⌋, as
Stata's `harreg`; an integer fixes it. A comparator of [`TLAHR`](@ref) (the paper's "t-HAR").
"""
struct CosineHAR <: PanelVcov
    rule::Symbol      # :llsw or :fixed
    B::Int
    function CosineHAR(rule::Symbol, B::Int)
        _check_option(rule, (:llsw, :fixed), "B")
        (rule === :llsw || B >= 1) ||
            throw(ArgumentError("CosineHAR: B must be at least 1; got $B"))
        return new(rule, B)
    end
end

# Per-horizon details of a variance computation.
struct VcovInfo
    factor::Float64       # small-sample factor applied
    groups::Int           # smallest number of clusters (0 if not clustered)
    periods::Int          # time periods entering a time aggregation (0 otherwise)
    lags::Int             # Driscoll–Kraay lags; CosineHAR: cosine terms B; 0 otherwise
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
