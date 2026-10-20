# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`, so this file is formatted
# by hand (as in LocalProjections.jl).
#! format: off

@testmodule PLPSupport begin
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs
    const SYNTHETIC = joinpath(@__DIR__, "fixtures", "synthetic")
    fixture(name::Symbol) = read_csv(joinpath(SYNTHETIC, string(name, ".csv")))
end

@testitem "synthetic fixtures regenerate" tags = [:core, :fixtures] setup = [PLPSupport] begin
    using .PLPSupport: PLPDGPs, fixture
    same(a, b) = (ismissing(a) && ismissing(b)) ||
                 (!ismissing(a) && !ismissing(b) && (a == b || isapprox(a, b; rtol = 1e-14, atol = 1e-12)))
    bad = String[]
    for name in keys(PLPDGPs.FIXTURES)
        disk = fixture(name)
        mem = PLPDGPs.make_fixture(name)
        @test [Symbol(k) for (k, _) in mem] == collect(keys(disk))
        for (k, v) in mem
            all(same.(v, disk[Symbol(k)])) || push!(bad, string(name, ".", k))
        end
    end
    @test isempty(bad)
end

@testitem "fixture design properties" tags = [:core, :fixtures] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    b = fixture(:fx_balanced)
    @test length(b.id) == 30 * 60
    u = fixture(:fx_unbalanced)
    counts = Dict{Int, Int}()
    foreach(i -> counts[i] = get(counts, i, 0) + 1, u.id)
    @test counts[38] == 2 && counts[39] == 1 && counts[40] == 1
    @test !any(==(35), u.t)
    @test eltype(u.w) == Union{Missing, Float64} && count(ismissing, u.w) > 0
    d = fixture(:fx_dated)
    @test eltype(d.date) <: PLPSupport.PLPDGPs.Date
end
