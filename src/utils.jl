# Small helpers shared across files.

# Float64 copy of a numeric column, with `missing` (and `nothing`) mapped to NaN.
function _nanfloat(v::AbstractVector)
    out = Vector{Float64}(undef, length(v))
    @inbounds for i in eachindex(v, out)
        x = v[i]
        out[i] = (x === missing || x === nothing) ? NaN : Float64(x)
    end
    return out
end

# Numeric content (Real, missing or nothing), whatever the element type (e.g. `Any`).
function _isnumericcol(v::AbstractVector)
    nonmissingtype(eltype(v)) <: Union{Nothing, Real} && return true
    return all(x -> x === missing || x === nothing || x isa Real, v)
end

# A missing value of a fixed-effect or cluster variable: `missing`, `nothing` or NaN, as
# for the formula variables (whose missing values become NaN).
_isna(x) = x === missing || x === nothing || (x isa AbstractFloat && isnan(x))

# Group codes 1:G in sorted order of the levels; 0 marks a missing value.
function _encode(v::AbstractVector)
    levels = sort!(unique(x for x in v if !_isna(x)))
    code = Dict(l => Int32(i) for (i, l) in enumerate(levels))
    return Int32[_isna(x) ? Int32(0) : code[x] for x in v], length(levels)
end

# Combine several code vectors into one (an interaction such as sector × quarter).
# 0 in any component gives 0.
function _combine_codes(codes::Vector{Vector{Int32}})
    length(codes) == 1 && return codes[1], Int(maximum(codes[1]; init = Int32(0)))
    n = length(codes[1])
    out = similar(codes[1])
    acc = copy(codes[1])
    G = 0
    for k in 2:length(codes)
        c = codes[k]
        d = Dict{NTuple{2, Int32}, Int32}()
        G = 0
        @inbounds for r in 1:n
            if acc[r] == 0 || c[r] == 0
                out[r] = 0
            else
                key = (acc[r], c[r])
                g = get(d, key, Int32(0))
                if g == 0
                    G += 1
                    g = d[key] = Int32(G)
                end
                out[r] = g
            end
        end
        acc .= out
    end
    return acc, G
end

# Row-wise Kronecker product of column blocks: the columns of an interaction term, with
# the first block varying slowest (StatsModels' `&` ordering).
function _rowkron(A::Matrix{Float64}, B::Matrix{Float64})
    n = size(A, 1)
    out = Matrix{Float64}(undef, n, size(A, 2) * size(B, 2))
    c = 0
    for j in axes(A, 2), k in axes(B, 2)

        c += 1
        @inbounds for r in 1:n
            out[r, c] = A[r, j] * B[r, k]
        end
    end
    return out
end

function _check_option(value::Symbol, allowed, name::AbstractString)
    value in allowed ||
        throw(ArgumentError("$name must be one of $(join(repr.(allowed), ", ")); got $(repr(value))"))
    return value
end

# All-finite test of row r of a matrix.
@inline function _rowfinite(X::Matrix{Float64}, r::Int)
    @inbounds for j in axes(X, 2)
        isfinite(X[r, j]) || return false
    end
    return true
end
