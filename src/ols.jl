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
function _eqr(X::Matrix{Float64}, kidx::Vector{Int})
    Xk = X[:, kidx]
    norms = [LinearAlgebra.norm(view(Xk, :, j)) for j in axes(Xk, 2)]
    Xk ./= norms'
    return LinearAlgebra.qr!(Xk, LinearAlgebra.ColumnNorm()), norms
end

# Columns absorbed by the fixed effects (or the constant), judged against the
# pre-absorption centred norm, before any rescaling.
function _absorbed_cols(Xt::AbstractMatrix{Float64}, Xraw::AbstractMatrix{Float64},
        hasfe::Bool, intercept::Bool)
    absorbed = falses(size(Xt, 2))
    for j in axes(Xt, 2)
        x = view(Xraw, :, j)
        μ = Statistics.mean(x)
        css = sum(v -> abs2(v - μ), x)
        if hasfe
            # any fixed effect spans the constant
            absorbed[j] = css == 0 || sum(abs2, view(Xt, :, j)) <= _ABSORBED_TOL * css
        elseif !(intercept && j == 1)
            absorbed[j] = sum(abs2, x) == 0 || (intercept && css == 0)
        end
    end
    return absorbed
end

"""
    _ols(Xt, yt, Xraw, hasfe, intercept, shockcols) -> NamedTuple

OLS of `yt` on `Xt` (within-transformed when `hasfe`). Returns the coefficients and bread
`(X̃'X̃)⁻¹` on the kept columns, the kept-column mask, the columns omitted as absorbed by
the fixed effects or the constant, the residuals, and `status`: `:ok`, or `:singular`
when a shock column is absorbed or collinear (the other fields are then empty).
"""
function _ols(Xt::Matrix{Float64}, yt::Vector{Float64}, Xraw::AbstractMatrix{Float64},
        hasfe::Bool, intercept::Bool, shockcols::Vector{Int})
    absorbed = _absorbed_cols(Xt, Xraw, hasfe, intercept)
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
        for j in kidx[F.p[(rank + 1):end]]
            keep[j] = false
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
