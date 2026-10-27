# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "PanelIndex shifts equal a (unit, period) lookup" tags = [:core, :panel] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    const PLP = PanelLocalProjections
    u = fixture(:fx_unbalanced)
    p = PLP.PanelIndex(u.id, u.t)
    lookup = Dict((u.id[r], u.t[r]) => r for r in eachindex(u.id))
    for k in -6:6
        oracle = Int32[get(lookup, (u.id[r], u.t[r] + k), 0) for r in eachindex(u.id)]
        @test PLP._shift(p, k) == oracle
    end
    y = Float64.(u.y)
    ly = PLP._lag(p, y, 2)
    @test all(r -> (j = get(lookup, (u.id[r], u.t[r] - 2), 0); j == 0 ? isnan(ly[r]) : ly[r] == y[j]),
        eachindex(y))
    @test PLP.nunits(p) == 40 && PLP.nspan(p) == 70 && PLP.nrows(p) == length(u.id)
end

@testitem "Date periods" tags = [:core, :panel] setup = [PLPSupport] begin
    using .PLPSupport: fixture
    using Dates
    const PLP = PanelLocalProjections
    for (f, d1, d2) in ((:year, Date(2003, 7, 9), Date(2004, 1, 1)),
                        (:halfyear, Date(2003, 6, 30), Date(2003, 7, 1)),
                        (:quarter, Date(2003, 3, 31), Date(2003, 4, 1)),
                        (:month, Date(2003, 1, 31), Date(2003, 2, 1)))
        @test PLP._period_code(d2, f) == PLP._period_code(d1, f) + 1
        c = PLP._period_code(d1, f)
        @test PLP._period_code(PLP._period_date(c, f), f) == c
    end
    d = fixture(:fx_dated)
    p = PLP.PanelIndex(d.id, d.date; frequency = :quarter)
    @test PLP.nspan(p) == 40 && PLP.nrows(p) == 239
    @test PLP._period_date(p.pmin, :quarter) == Date(2001, 1, 1)
    # unit 2 has no 2005Q3: its 2005Q2 lead is undefined, as is every unit's last lead
    @test count(==(0), PLP._shift(p, 1)) == 6 + 1
    @test_throws ArgumentError PLP.PanelIndex(d.id, d.date)            # no frequency
    @test_throws ArgumentError PLP.PanelIndex(d.id, d.date; frequency = :week)
end

@testitem "PanelIndex errors and warnings" tags = [:core, :panel] begin
    const PLP = PanelLocalProjections
    @test_throws ArgumentError PLP.PanelIndex([1, 1, 2], [1, 1, 1])     # duplicate (id, t)
    @test_throws ArgumentError PLP.PanelIndex([1, missing], [1, 2])     # missing id
    @test_throws ArgumentError PLP.PanelIndex([1, 2], [1, missing])     # missing time
    @test_throws ArgumentError PLP.PanelIndex([1, 2], [1.0, 1.5])       # non-integral time
    @test_throws ArgumentError PLP.PanelIndex([1, 2], [1, 2]; frequency = :year)
    @test_throws ArgumentError PLP.PanelIndex(1:3, [1, 2, 10^10])       # span too large
    # unsupported time types: an error, never an endless recursion
    using Dates
    for tm in (["2001", "2002"], [:a, :b], ['a', 'b'], [Year(2001), Year(2002)], Any[1, "2"])
        @test_throws ArgumentError PLP.PanelIndex([1, 1], tm)
    end
    @test PLP.PanelIndex([1, 1], Any[1, 2]).period == [1, 2]
    # frequency as a string; duplicates are reported by period
    dts = [Date(2001, 1, 5), Date(2001, 2, 5), Date(2001, 4, 5)]
    @test PLP.PanelIndex([1, 1, 1], dts; frequency = "month").period ==
          PLP.PanelIndex([1, 1, 1], dts; frequency = :month).period
    @test_throws ArgumentError PLP.PanelIndex([1, 1], dts[1:2]; frequency = "week")
    err = try
        PLP.PanelIndex([1, 1, 1], dts; frequency = :quarter)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("(1, 2001Q1)", err.msg)
    p = PLP.PanelIndex(["b", "a", "b"], [1.0, 2.0, 2.0])
    @test p.unitlabels == ["a", "b"] && p.unit == Int32[2, 1, 2] && p.period == [1, 2, 2]
    # YYYYMM codes: most of the span is never observed
    @test_logs (:warn, r"absent for every unit") PLP.PanelIndex(ones(Int, 24),
        [200001 + 100 * (k ÷ 12) + k % 12 for k in 0:23])
end
