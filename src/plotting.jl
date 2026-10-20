# Plot recipes. RecipesBase (Plots.jl) recipes live here; Makie methods for `irfplot` and
# `irfplot!` are defined in ext/PanelLocalProjectionsMakieExt.jl.

"""
    irfplot(s::IRFSummary; kw...)
    irfplot!(ax, s::IRFSummary; kw...)

Impulse response with its confidence band, drawn with Makie. Load Makie (e.g. `using
CairoMakie`) to enable these methods. With Plots.jl, `plot(s)` works directly.
"""
function irfplot end
@doc (@doc irfplot) function irfplot! end

@recipe function f(s::IRFSummary)
    xguide --> "horizon"
    yguide --> s.response
    title --> s.term
    legend --> false
    @series begin
        seriestype := :path
        ribbon := (s.coef .- s.lower, s.upper .- s.coef)
        fillalpha --> 0.25
        linewidth --> 2
        s.horizon, s.coef
    end
    @series begin
        seriestype := :hline
        linestyle := :dash
        linecolor := :black
        linewidth := 1
        primary := false
        [0.0]
    end
end
