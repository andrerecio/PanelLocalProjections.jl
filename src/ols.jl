# Least squares on the within-transformed design of one horizon.
#
# Approach as in Regress.jl (G. Ragusa; src/ols_solver.jl: pivoted QR with a relative rank
# tolerance), with two differences: columns absorbed by the fixed effects are screened
# against their *pre-absorption* centered norm before any rescaling (rescaling would turn
# round-off into a unit-norm column that QR accepts), and an absorbed or collinear shock
# makes the horizon singular (NaN) instead of getting a zero coefficient.

const _ABSORBED_TOL = 1e-8   # ‖x̃‖² ≤ tol · ‖x − x̄‖² (squared norms): the column is absorbed
const _RANK_TOL = 1e-9       # |R_jj| ≤ tol · |R_11| on unit-norm columns: collinear

# Pivoted QR of the columns `kidx` of `X`, each scaled to unit norm; returns the norms too.
function _eqr(X::AbstractMatrix{Float64}, kidx::Vector{Int})
    Xk = X[:, kidx]
    norms = [LinearAlgebra.norm(view(Xk, :, j)) for j in axes(Xk, 2)]
    Xk ./= norms'
    return LinearAlgebra.qr!(Xk, LinearAlgebra.ColumnNorm()), norms
end

# Per column of the raw (pre-absorption) regressors: the centred and the raw sums of
# squares, what the absorbed-column screen and the PSD metric need, so that the raw
# columns need not be kept once absorbed.
function _colstats(X::AbstractMatrix{Float64})
    css = zeros(size(X, 2))
    ss = zeros(size(X, 2))
    for j in axes(X, 2)
        x = view(X, :, j)
        μ = Statistics.mean(x)
        css[j] = sum(v -> abs2(v - μ), x)
        ss[j] = sum(abs2, x)
    end
    return (; n = size(X, 1), css, ss)
end

# Columns absorbed by the fixed effects (or the constant), judged against the
# pre-absorption centred norm, before any rescaling. `raw`: the raw columns, or their
# `_colstats`.
function _absorbed_cols(Xt::AbstractMatrix{Float64}, raw::AbstractMatrix{Float64},
        hasfe::Bool, intercept::Bool)
    return _absorbed_cols(Xt, _colstats(raw), hasfe, intercept)
end
function _absorbed_cols(Xt::AbstractMatrix{Float64}, raw::NamedTuple, hasfe::Bool,
        intercept::Bool)
    absorbed = falses(size(Xt, 2))
    for j in axes(Xt, 2)
        css = raw.css[j]
        if hasfe
            # any fixed effect spans the constant
            absorbed[j] = css == 0 || sum(abs2, view(Xt, :, j)) <= _ABSORBED_TOL * css
        elseif !(intercept && j == 1)
            absorbed[j] = raw.ss[j] == 0 || (intercept && css == 0)
        end
    end
    return absorbed
end

"""
    _ols(Xt, yt, raw, hasfe, intercept, shockcols, added = Int[]; absorbed) -> NamedTuple

OLS of `yt` on `Xt` (within-transformed when `hasfe`; `raw`: the raw columns or their
`_colstats`). `absorbed` (by default `_absorbed_cols(Xt, raw, …)`) is passed when `Xt` is
not the absorbed data itself (the second stage of 2SLS, on fitted values). Returns the
coefficients and bread
`(X̃'X̃)⁻¹` on the kept columns, the kept-column mask, the columns omitted as absorbed by
the fixed effects or the constant, the residuals, and `status`: `:ok`, or `:singular`
when a shock column is absorbed or collinear (the other fields are then empty). Among
collinear columns, those in `added` (lags added by `LagAugmented`) are omitted first.
"""
function _ols(Xt::AbstractMatrix{Float64}, yt::Vector{Float64}, raw, hasfe::Bool,
        intercept::Bool, shockcols::Vector{Int}, added::Vector{Int} = Int[];
        absorbed::BitVector = _absorbed_cols(Xt, raw, hasfe, intercept))
    keep = .!absorbed
    singular = (; β = Float64[], B = zeros(0, 0), keep, absorbed, resid = Float64[],
        status = :singular)
    any(j -> absorbed[j], shockcols) && return singular
    kidx = findall(keep)
    isempty(kidx) && return singular
    F, norms = _eqr(Xt, kidx)
    d = abs.(LinearAlgebra.diag(F.R))
    rank = count(>(_RANK_TOL * d[1]), d)
    if rank < length(kidx)
        if any(j -> keep[j], added)
            _dropcollinear!(keep, Xt, added)
        else
            for j in kidx[F.p[(rank + 1):end]]
                keep[j] = false
            end
        end
        any(j -> !keep[j], shockcols) && return singular
        kidx = findall(keep)
        F, norms = _eqr(Xt, kidx)
    end
    β = (F \ yt) ./ norms
    resid = yt - view(Xt, :, kidx) * β
    # (X̃k'X̃k)⁻¹ = D⁻¹ P (R'R)⁻¹ P' D⁻¹, D = diag(norms)
    Rinv = inv(LinearAlgebra.UpperTriangular(F.R))
    A = Rinv * Rinv'
    p = F.p
    B = similar(A)
    B[p, p] = A
    B ./= norms .* norms'
    return (; β, B, keep, absorbed, resid, status = :ok)
end

# Collinear columns among the kept ones, omitting those in `added` rather than the others:
# the pivoted QR of the other columns first; then that of the `added` columns (unit norm)
# residualised on them, at the same tolerance.
function _dropcollinear!(keep::BitVector, Xt::AbstractMatrix{Float64}, added::Vector{Int})
    isadded = falses(length(keep))
    isadded[added] .= true
    fidx = findall(keep .& .!isadded)
    if !isempty(fidx)
        F, _ = _eqr(Xt, fidx)
        d = abs.(LinearAlgebra.diag(F.R))
        for j in fidx[F.p[(count(>(_RANK_TOL * d[1]), d) + 1):end]]
            keep[j] = false
        end
        fidx = findall(keep .& .!isadded)
    end
    aidx = findall(keep .& isadded)
    R = Xt[:, aidx]
    for c in axes(R, 2)
        nc = LinearAlgebra.norm(view(R, :, c))
        nc > 0 && (R[:, c] ./= nc)
    end
    if !isempty(fidx)
        # in the coordinates of the Q of the kept columns, the residuals on them are the
        # rows after the first length(fidx)
        F, _ = _eqr(Xt, fidx)
        R = (F.Q' * R)[(length(fidx) + 1):end, :]
    end
    Fa = LinearAlgebra.qr!(R, LinearAlgebra.ColumnNorm())
    d = abs.(LinearAlgebra.diag(Fa.R))
    for j in aidx[Fa.p[(count(>(_RANK_TOL), d) + 1):end]]
        keep[j] = false
    end
    return keep
end
