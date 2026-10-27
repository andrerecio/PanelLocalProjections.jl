# Minimal CSV reading for the US state example (no package dependency).

# Fields of one CSV line (RFC 4180: quoted fields may contain commas; "" is a quote).
function csvfields(line::AbstractString)
    out = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq
            if c == '"'
                j = nextind(line, i)
                if j <= lastindex(line) && line[j] == '"'
                    write(buf, '"')
                    i = j
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(out, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    push!(out, String(take!(buf)))
    return out
end


"""
    readtable(path) -> NamedTuple

A CSV file with a header row as a column table; columns whose fields all parse as integers
become `Vector{Int}`, as numbers `Vector{Float64}` (`NaN` for empty or non-numeric fields
in a numeric column), the rest `Vector{String}`.
"""
function readtable(path)
    lines = filter(!isempty, readlines(path))
    hdr = Symbol.(csvfields(lines[1]))
    rows = [csvfields(l) for l in lines[2:end]]
    function column(j)
        v = [r[j] for r in rows]
        all(x -> tryparse(Int, x) !== nothing, v) && return parse.(Int, v)
        f = [something(tryparse(Float64, x), NaN) for x in v]
        return count(isnan, f) < length(f) ? f : v
    end
    return NamedTuple{Tuple(hdr)}(Tuple(column(j) for j in eachindex(hdr)))
end
