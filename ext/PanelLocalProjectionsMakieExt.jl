module PanelLocalProjectionsMakieExt

# Makie methods for `irfplot`/`irfplot!` (ported in spirit from LocalProjections.jl's
# LocalProjectionsMakieExt; here they draw an IRFSummary).

using PanelLocalProjections: PanelLocalProjections, IRFSummary
import PanelLocalProjections: irfplot, irfplot!
using Makie: Makie

"""
    irfplot!(ax, s::IRFSummary; color = :steelblue, bandalpha = 0.25, zeroline = true)

Draw the impulse response `s` with its confidence band on the Makie axis `ax`.
"""
function irfplot!(ax::Makie.Axis, s::IRFSummary; color = :steelblue, bandalpha = 0.25,
        zeroline::Bool = true)
    zeroline && Makie.hlines!(ax, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
    Makie.band!(ax, s.horizon, s.lower, s.upper; color = (color, bandalpha))
    Makie.lines!(ax, s.horizon, s.coef; color = color, linewidth = 2)
    return ax
end

"""
    irfplot(s::IRFSummary; figure = (;), axis = (;), kw...) -> Makie.Figure

New figure with the impulse response `s`; `figure` and `axis` are passed to
`Makie.Figure` and `Makie.Axis`, the other keywords to [`irfplot!`](@ref).
"""
function irfplot(s::IRFSummary; figure = (;), axis = (;), kw...)
    fig = Makie.Figure(; figure...)
    pct = round(100 * s.level; digits = 1)
    ax = Makie.Axis(fig[1, 1]; xlabel = "horizon", ylabel = s.response,
        title = "$(s.term) ($(pct)% band)", axis...)
    irfplot!(ax, s; kw...)
    return fig
end

end # module
