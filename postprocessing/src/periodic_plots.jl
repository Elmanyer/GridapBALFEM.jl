# ==============================================================
#  periodic_plots.jl — figures and animations of a periodic-box run (see periodic.jl)
#
#  State:   plot_state (all fields at one instant) · plot_state_wiggle (layers drawn inside the column) · animate_run (time marching, frame selection)
#  Modes:   plot_mode_timeseries · plot_mode_cycle · plot_mode_phasors · plot_mode_amplitudes
#           plot_mode_profiles · animate_mode
#  Space–time: plot_layers_hovmoller · plot_steady_split · plot_column · plot_komega
#  Energy:  plot_energy · plot_energy_spectrum · plot_growth_overview
#  All:     periodic_report (every figure + an animation into one folder)
# ==============================================================

const LAYER_COLORS = [:royalblue, :seagreen, :darkorange, :purple, :brown, :olive, :magenta]
const _M = (left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, top_margin = 2Plots.mm)
# diverging map for signed fields: red = positive (+x), blue = negative. Plots' :RdBu runs red → blue, low → high.
const SIGNED_CMAP = cgrad(:RdBu; rev = true)

layer_color(j) = LAYER_COLORS[mod1(j, length(LAYER_COLORS))]
function layer_name(r::PeriodicRun, j::Int)
    σ = _node_sigma(r)[j]
    j == 1 ? "u1 (bed)" : j == nlayers(r) ? "u$j (surface)" : @sprintf("u%d (σ=%.3f)", j, σ)
end
"σ position of each vertical node (mode j ↔ node j)."
function _node_sigma(r::PeriodicRun)
    b = r.basis; σ = zeros(b.Nσ)
    for (gid, loc) in zip(b.gnodes, b.local_nodes), (g, s) in zip(gid, loc)
        σ[g] = s
    end
    return σ
end
_rms(v) = vec(sqrt.(mean(v .^ 2; dims = 2)))
_symlim(a) = (a = a == 0 || !isfinite(a) ? 1e-12 : a; (-a, a))

# ---- state at one instant ----------------------------------------------------------------------------

_poslog(v, lo) = [x > lo ? x : NaN for x in v]          # values at or below the axis floor are not drawn
_kin(S) = S.mean .+ S.shear .+ S.w

"""Fixed log limits of the two spectrum panels over a set of snapshots, `(surface, velocity)` (used by `animate_run`)."""
function spectrum_limits(r::PeriodicRun, its; avg::Symbol = :period)
    S0 = state_spectrum(r, 1; avg); Ss = [state_spectrum(r, it; avg) for it in its]
    lim(v0, hi) = (minimum(v0[v0 .> 0]) / 10, 3hi)
    return (lim(S0.pot, maximum(maximum(S.pot) for S in Ss)), lim(_kin(S0), maximum(maximum(_kin(S)) for S in Ss)))
end

# totals of a state_spectrum: E/E0 (with w), surface share, the x-uniform flow's share, the disturbance's
# share and its dominant kind
function _spectrum_summary(S, S0, nmin)
    E = sum(S.pot) + sum(_kin(S)); E0 = sum(S0.pot) + sum(_kin(S0))
    dη = S.nη .≥ nmin; du = S.nu .≥ nmin
    kinds = ("surface" => sum(S.pot[dη]), "depth-mean" => sum(S.mean[du]), "shear" => sum(S.shear[du]), "w" => sum(S.w[du]))
    Ed = sum(last, kinds); kd = argmax(last, kinds)
    return (E = E / E0, pot = sum(S.pot) / E, n0 = (S.pot[1] + _kin(S)[1]) / E, dist = Ed / E, kind = first(kd),
            kindfrac = Ed > 0 ? last(kd) / Ed : 0.0)
end
# the velocity panel's legend, drawn in another (frameless) subplot
function _velocity_legend!(p, r)
    nan = [NaN]
    plot!(p, nan, nan; c = :gray, ls = :dash, lw = 1, label = "kinetic at t = 0")
    plot!(p, nan, nan; c = :gray, lw = 3.5, alpha = 0.45, label = "kinetic total")
    for (c, lb) in ((:royalblue, "depth-mean"), (:darkorange, "shear"), (:seagreen, "w"))
        plot!(p, nan, nan; c, lw = 1.2, m = :circle, ms = 1.5, msw = 0, label = lb)
    end
    scatter!(p, nan, nan; m = :star5, ms = 7, c = :red, msw = 0, label = "n = 0: x-uniform flow")
    plot!(p, nan, nan; c = :purple, ls = :dot, lw = 1, label = "n = $(r.meta.ncell): 1 λ per cell")
    plot!(p, nan, nan; c = :gold, alpha = 0.4, lw = 8, label = "band of the part shown")
    return p
end
_summary_line(Σ, nmin) = @sprintf("E/E0 = %.3f (with w) · n≥%d: %.2g %%, %.0f %% %s", Σ.E, nmin, 100Σ.dist, 100Σ.kindfrac, Σ.kind)

# shade the band of `part`, mark the carrier/disturbance divide and n = ncell (one wavelength per cell)
function _mark_bands!(p, r, part, nmin, ncarrier, nmax)
    band = part === :total ? nothing : part === :carrier ? (1, ncarrier) : part isa Symbol ? (nmin, nmax) : extrema(part)
    band === nothing || band[1] > nmax || vspan!(p, [band[1] - 0.5, min(band[2], nmax) + 0.5]; c = :gold, alpha = 0.2, lw = 0, label = "")
    vline!(p, [ncarrier + 0.5]; c = :gray, lw = 0.6, label = "")
    r.meta.ncell < nmax && vline!(p, [r.meta.ncell]; c = :purple, ls = :dot, lw = 1, label = "n = $(r.meta.ncell): 1 λ per cell")
    return p
end

# beside η: surface energy ½gη² per η mode, n = 0 … ncell·p_eta/2 (from the nodal values: no images)
function _surface_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims, kw...)
    n = collect(S.nη); nmax = last(n)
    yl = something(ylims, (minimum(S0.pot[S0.pot .> 0]) / 10, 3maximum(S.pot)))
    p = plot(; kw..., yscale = :log10, ylims = yl, xlims = (-0.5, nmax + 0.5), xlabel = "η mode n", ylabel = "energy (m³/s²)",
             legend = :topright, legendfontsize = 7, titlefontsize = 9,
             title = @sprintf("surface energy per η mode (Q%d: n ≤ %d)", r.meta.p_eta, nmax))
    _mark_bands!(p, r, part, nmin, ncarrier, nmax)
    plot!(p, n, _poslog(S0.pot, yl[1]); c = :gray, ls = :dash, lw = 1, label = "t = 0")
    plot!(p, n, _poslog(S.pot, yl[1]); c = :black, lw = 1.5, m = :circle, ms = 2, msw = 0, label = "surface ½gη²")
    return p
end

# beside the water column: kinetic energy per velocity mode, n = 0 … ncell·p_u/2, split into depth-mean,
# shear and w (from the nodal values: no images)
function _velocity_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims, title, legend = :outerbottom, kw...)
    n = collect(S.nu); nmax = last(n); K = _kin(S); K0 = _kin(S0)
    yl = something(ylims, (minimum(K0[K0 .> 0]) / 10, 3maximum(K)))
    p = plot(; xlabel = "velocity mode n", kw..., yscale = :log10, ylims = yl, xlims = (-0.5, nmax + 0.5), ylabel = "energy (m³/s²)",
             title, titlefontsize = 9, legend, legendfontsize = 7, legend_columns = 2)
    _mark_bands!(p, r, part, nmin, ncarrier, nmax)
    plot!(p, n, _poslog(K0, yl[1]); c = :gray, ls = :dash, lw = 1, label = "kinetic at t = 0")
    plot!(p, n, _poslog(K, yl[1]); c = :gray, lw = 3.5, alpha = 0.45, label = "kinetic total")
    for (v, c, lb) in ((S.mean, :royalblue, "depth-mean"), (S.shear, :darkorange, "shear"), (S.w, :seagreen, "w"))
        plot!(p, n, _poslog(v, yl[1]); c, lw = 1.2, m = :circle, ms = 1.5, msw = 0, label = lb)
    end
    scatter!(p, [0], _poslog(K[1:1], yl[1]); m = :star5, ms = 7, c = :red, msw = 0, label = "n = 0: x-uniform flow")
    return p
end

"""
    plot_state(r, t; part=:total, normalise=:auto, streamlines=true, zw=(0, d), avg=:period,
               lims=nothing, elims=nothing, timeline=nothing, nmin=5, ncarrier=4) -> Plot

Everything at one instant `t` (or snapshot index `it = t::Int`), for one `part` of the solution
(`:total`, `:carrier`, `:disturbance`, `:steady`, `:oscillating`, or a harmonic list such as `[16]`):
* top: the free surface η(x), and beside it the surface energy ½gη² per η mode (n ≤ ncell·p_eta/2);
* middle: every layer velocity u_j(x) and the depth average ū (black, dashed), and beside it the energy
  totals: E/E0 (with w), the surface share, the x-uniform flow (n = 0), the disturbance's share and kind;
* bottom: the velocity inside the water column, u(x,σ) (colour; red = +x) with streamlines, and beside it
  the kinetic energy per velocity mode (n ≤ ncell·p_u/2), split into depth-mean, shear and w.
The spectra (`state_spectrum`) are of the whole solution, from the nodal values (each mode once, no
images), against t = 0 (dashed); the band of `part` is shaded. `avg = :period` (default) averages them
over one wave period centred on t, `:instant` does not.
`normalise` (default: on for the disturbance parts) divides every panel by its own maximum, printed
in the titles. `lims = (η, u)` and `elims = (surface, velocity)` (spectra) fix the axes instead (used by
`animate_run`).
`timeline = (t, y)` adds a strip with the curve y(t) and the current time marked.
"""
function plot_state(r::PeriodicRun, t::Real; part = :total, normalise = :auto, streamlines::Bool = true,
                    zw = (0.0, r.meta.d), avg::Symbol = :period, lims = nothing, elims = nothing, timeline = nothing, timeline_label = "",
                    nmin::Int = 5, ncarrier::Int = 4, size = (1350, 950))
    it = t isa Integer ? t : tindex(r, t); N = nlayers(r); d = r.meta.d
    η, U, ub = part_state(r, it, part; nmin, ncarrier)
    nrm = normalise === :auto ? !(part in (:total, :carrier)) : normalise
    aη = maximum(abs.(η)); au = maximum(abs.(U))
    sη, su = nrm ? (aη, au) : (1.0, 1.0)
    ylη = lims === nothing ? _symlim(1.1 * aη / sη) : _symlim(lims[1]); ylu = lims === nothing ? _symlim(1.1 * au / su) : _symlim(lims[2])
    plab = part_label(part; nmin, ncarrier)

    p1 = plot(r.x, η ./ sη; c = :black, lw = 2, legend = false, ylims = ylη, xlims = extrema(r.x), ylabel = "η" * (nrm ? " / max" : " (m)"),
              title = @sprintf("t = %.2f s   ·   %s   ·   free surface (max |η| %.2e m)", r.t[it], plab, aη), titlefontsize = 10)
    p2 = plot(; legend = :outerright, ylims = ylu, xlims = extrema(r.x), ylabel = "u" * (nrm ? " / max" : " (m/s)"),
              title = @sprintf("layer velocities (max |u_j| %.2e m/s)", au), titlefontsize = 10, legendfontsize = 7)
    for j in 1:N
        plot!(p2, r.x, U[j, :] ./ su; c = layer_color(j), lw = 1.8, label = layer_name(r, j))
    end
    plot!(p2, r.x, ub ./ su; c = :black, lw = 2.2, ls = :dash, label = "ū")
    σs, u, ψ = column_field(r, U)
    z = σs .* d; iz = findall(zz -> zw[1] ≤ zz ≤ zw[2], z)
    ac = maximum(abs.(u[iz, :]))
    p3 = heatmap(r.x, z[iz], u[iz, :] ./ (nrm || lims === nothing ? ac : lims[2]); c = SIGNED_CMAP, clims = (-1, 1), colorbar = false,
                 xlabel = "x (m)", ylabel = "height above bed (m)", xlims = extrema(r.x), ylims = (z[iz[1]], z[iz[end]]),
                 title = @sprintf("inside the water column: u(x,σ) (colour, red = +x) and streamlines (max |u| %.2e m/s)", ac),
                 titlefontsize = 10)
    if streamlines
        b = maximum(abs.(ψ[iz, :]))
        b > 0 && contour!(p3, r.x, z[iz], ψ[iz, :]; levels = collect(range(-b, b; length = 12))[2:end-1], c = :black, lw = 0.6, colorbar = false)
    end
    for s in _node_sigma(r)[2:end-1]
        hline!(p3, [s * d]; c = :gray, ls = :dash, lw = 0.8, label = "")
    end
    S = state_spectrum(r, it; avg); S0 = state_spectrum(r, 1; avg); Σ = _spectrum_summary(S, S0, nmin); el = something(elims, (nothing, nothing))
    ps = _surface_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims = el[1])
    pu = _velocity_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims = el[2], legend = false,
                                  title = @sprintf("kinetic energy per velocity mode (Q%d: n ≤ %d)", r.meta.p_u, last(S.nu)))
    ptxt = plot(; framestyle = :none, grid = false, xlims = (0, 1), ylims = (0, 1), legend = :bottom, legendfontsize = 7,
                legend_columns = 2, foreground_color_legend = nothing)
    annotate!(ptxt, 0.02, 0.8, text(@sprintf("energy (%s)\nE/E0 = %.3f  (with w)  ·  surface %.0f %%\nn = 0 (x-uniform)  %.2g %%\nn ≥ %d  %.2g %%, of which %.0f %% %s",
                                         avg === :period ? "period mean" : "instant", Σ.E, 100Σ.pot, 100Σ.n0, nmin, 100Σ.dist,
                                         100Σ.kindfrac, Σ.kind), 9, :left))
    _velocity_legend!(ptxt, r)
    rows = [plot(a, b; layout = grid(1, 2; widths = [0.72, 0.28])) for (a, b) in ((p1, ps), (p2, ptxt), (p3, pu))]
    top = plot(rows...; layout = grid(3, 1; heights = [0.2, 0.3, 0.5]))
    timeline === nothing && return plot(top; size = size, _M...)
    tl, yl = timeline
    p5 = plot(tl, yl; yscale = :log10, c = :black, lw = 1.5, legend = false, xlabel = "t (s)", ylabel = timeline_label,
              titlefontsize = 9, xlims = (r.t[1], r.t[end]))
    vline!(p5, [r.t[it]]; c = :red, lw = 1.5)
    return plot(top, p5; layout = grid(2, 1; heights = [0.84, 0.16]), size = size .+ (0, 150), _M...)
end

"""
    plot_state_wiggle(r, t; part=:total, normalise=:auto, streamlines=true, amp=0.45, avg=:period,
                      elims=nothing, nmin=5, ncarrier=4) -> Plot

`plot_state` with the layer-velocity and water-column panels merged into one (σ, u) panel:
* top: the free surface η(x);
* middle: u(x,σ) (colour, red = +x) with streamlines, against **σ on the left axis**; on top of it each layer
  velocity u_j(x) is drawn **on its own node**, as the curve σ_j + s·u_j(x) (black, white halo) about a dotted
  baseline at σ_j. All layers share one scale s: the largest |u_j| fills `amp` of the smallest node gap. The
  **right axis reads u_j in m/s** (u/max when normalised) about each baseline: every node carries its own 0
  and ±U ticks. The bed and surface curves swing outside [0, 1], into small margins;
* right: beside η the surface energy per η mode, beside the merged panel the kinetic energy per velocity
  mode, as in `plot_state` (the totals in the second title).
Keywords as in `plot_state`.
"""
function plot_state_wiggle(r::PeriodicRun, t::Real; part = :total, normalise = :auto, streamlines::Bool = true,
                           amp::Real = 0.45, avg::Symbol = :period, elims = nothing, nmin::Int = 5, ncarrier::Int = 4,
                           size = (1350, 850))
    it = t isa Integer ? t : tindex(r, t); N = nlayers(r)
    η, U, ub = part_state(r, it, part; nmin, ncarrier)
    nrm = normalise === :auto ? !(part in (:total, :carrier)) : normalise
    aη = maximum(abs.(η)); au = maximum(abs.(U)); su = nrm ? au : 1.0
    σn = _node_sigma(r)
    plab = part_label(part; nmin, ncarrier)

    # one common scale for every layer curve: the largest |u_j| fills `amp` of the smallest node gap
    gap = minimum(diff(σn)); s = amp * gap / (au / su)          # σ per (m/s) [or per unit of u/max]
    pad = 1.08 * s * au / su; ylσ = (-pad, 1 + pad)

    p1 = plot(r.x, η ./ (nrm ? aη : 1.0); c = :black, lw = 2, legend = false, xlims = extrema(r.x),
              ylims = _symlim(1.1 * (nrm ? 1.0 : aη)), ylabel = "η" * (nrm ? " / max" : " (m)"),
              title = @sprintf("t = %.2f s   ·   %s   ·   free surface (max |η| %.2e m)", r.t[it], plab, aη), titlefontsize = 10)

    # merged panel: u(x,σ) colour + streamlines, layer curves on their own σ baselines
    σs, u, ψ = column_field(r, U)
    ac = maximum(abs.(u))
    p2 = heatmap(r.x, σs, u ./ ac; c = SIGNED_CMAP, clims = (-1, 1), colorbar = false, xlims = extrema(r.x), ylims = ylσ,
                 xlabel = "x (m)", ylabel = "σ", yticks = 0:0.2:1, legend = false,
                 title = @sprintf("u(x,σ) (colour, red = +x; max %.2e m/s), streamlines, and each u_j(x) on its node σ_j (max |u_j| %.2e m/s)", ac, au),
                 titlefontsize = 10)
    if streamlines
        b = maximum(abs.(ψ))
        b > 0 && contour!(p2, r.x, σs, ψ; levels = collect(range(-b, b; length = 12))[2:end-1], c = :gray30, lw = 0.5, colorbar = false)
    end
    hline!(p2, [0.0, 1.0]; c = :black, lw = 1)                               # bed and surface
    for j in 1:N
        hline!(p2, [σn[j]]; c = :black, lw = 0.6, ls = :dot)                  # baseline u_j = 0
        y = σn[j] .+ s .* U[j, :] ./ su
        plot!(p2, r.x, y; c = :white, lw = 4.5)                               # halo
        plot!(p2, r.x, y; c = :black, lw = 1.8)
        annotate!(p2, r.x[1] + 0.01 * r.meta.Lx, σn[j] + (j == N ? -0.035 : 0.035), text(layer_name(r, j), 7, :left, :black))
    end
    # right axis: 0 and ±U ticks at every baseline, U a round number below max |u_j|
    Umax = au / su; U0 = 10.0^floor(log10(0.8Umax)); Ut = U0 * (0.8Umax ÷ U0)
    tk = Float64[]; tl = String[]
    for j in 1:N, v in (-Ut, 0.0, Ut)
        push!(tk, σn[j] + s * v); push!(tl, v == 0 ? "0" : @sprintf("%+.2g", v))
    end
    pr = twinx(p2)
    plot!(pr, [NaN], [NaN]; ylims = ylσ, yticks = (tk, tl), xlims = extrema(r.x), xticks = false, legend = false,
          ylabel = "u_j" * (nrm ? " / max" : " (m/s)"), ytickfontsize = 7)

    S = state_spectrum(r, it; avg); S0 = state_spectrum(r, 1; avg); el = something(elims, (nothing, nothing))
    ps = _surface_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims = el[1], left_margin = 14Plots.mm, bottom_margin = 4Plots.mm)
    pu = _velocity_spectrum_panel(r, S, S0; part, nmin, ncarrier, ylims = el[2], left_margin = 14Plots.mm, xlabel = "",
                                  title = @sprintf("kinetic energy per velocity mode (Q%d: n ≤ %d)\n%s", r.meta.p_u, last(S.nu),
                                                   _summary_line(_spectrum_summary(S, S0, nmin), nmin)))
    # margins per panel, set before the rows are assembled (a margin given to the outer plot would
    # override these on every subplot)
    plot!(p1; left_margin = 5Plots.mm); plot!(p2; left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, right_margin = 4Plots.mm)
    row1 = plot(p1, ps; layout = grid(1, 2; widths = [0.70, 0.30]))
    row2 = plot(p2, pu; layout = grid(1, 2; widths = [0.70, 0.30]))
    return plot(row1, row2; layout = grid(2, 1; heights = [0.22, 0.78]), size = size, top_margin = 2Plots.mm)
end

"""
    animate_run(r; frames=select_frames(r; dt=1.0), part=:total, out="<run>/viz/state_<part>.mp4",
                fps=12, timeline=:disturbance, kw...) -> out

Time marching: one `plot_state` per frame. Choose the frames with `select_frames` (time range,
step, explicit times, number of frames, or a stroboscopic phase). For `:total` and `:carrier` the
axes are fixed over the frames; for the disturbance parts every frame is normalised (see
`plot_state`). The spectrum panel keeps one log scale over all frames. `timeline`: `:disturbance` (rms of the surface-layer disturbance), `:energy` (monitor
energy ratio) or `nothing`. `out` ending in `.gif` writes a GIF (large), otherwise `.mp4`.
"""
function animate_run(r::PeriodicRun; frames = select_frames(r; dt = 1.0), part = :total, out = nothing,
                     fps::Int = 12, timeline = :disturbance, nmin::Int = 5, ncarrier::Int = 4, kw...)
    tag = part isa Symbol ? String(part) : "n" * join(string.(collect(part)[[1, end]]), "-")
    out = something(out, joinpath(mkpath(joinpath(r.dir, "viz")), "state_$tag.mp4"))
    lims = nothing
    if part in (:total, :carrier)
        rowsη = part_rows(r, r.η, frames, part; nmin, ncarrier)
        rowsU = [part_rows(r, r.U[j, :, :], frames, part; nmin, ncarrier) for j in 1:nlayers(r)]
        lims = (1.1 * maximum(abs.(rowsη)), 1.1 * maximum(maximum(abs.(u)) for u in rowsU))
    end
    elims = spectrum_limits(r, frames)
    tl = timeline === :disturbance ? ((r.t, _rms(xband(r.U[nlayers(r), :, :], nmin:nharmonics(r)))), "rms u_surf, n≥$nmin") :
         timeline === :energy && r.diag !== nothing ? ((r.diag["t"], abs.(r.diag["energy_ratio"])), "E / E0") : (nothing, "")
    anim = @animate for it in frames
        plot_state(r, it; part, lims, elims, timeline = tl[1], timeline_label = tl[2], nmin, ncarrier, kw...)
    end
    endswith(out, ".gif") ? gif(anim, out; fps) : mp4(anim, out; fps)
    return out
end

"Animation of one harmonic (or a set) alone: `animate_run` with `part = [n]` (normalised per frame)."
animate_mode(r::PeriodicRun, n; frames = select_frames(r; dt = 1.0), kw...) =
    animate_run(r; frames, part = n isa Integer ? [n] : n, kw...)

# ---- modes ---------------------------------------------------------------------------------------------

"""
    plot_mode_timeseries(r, ns; tw, labels) -> Plot

Each harmonic n seen at x = 0, one line per layer plus ū, over the window `tw` (pick a few
periods); grey lines one wave period apart. Shows whether a harmonic oscillates (and how fast) or
stands still, and whether the layers move together or against each other.
"""
function plot_mode_timeseries(r::PeriodicRun, ns; tw = (r.t[end] - 4r.meta.T, r.t[end]), labels = string.(ns))
    its = twindow(r, tw); ub = depth_average(r); N = nlayers(r)
    ps = map(enumerate(zip(ns, labels))) do (k, (n, lb))
        ys = [real.(xharmonic(r.U[j, its, :], n)) for j in 1:N]; yb = real.(xharmonic(ub[its, :], n))
        a = maximum(maximum(abs.(y)) for y in ys)
        p = plot(; xlabel = "t (s)", ylabel = "u / max", title = "n = $n: $lb", titlefontsize = 10,
                 ylims = (-1.1, 1.1), legend = k == 1 ? :outerright : false)
        for j in 1:N
            plot!(p, r.t[its], ys[j] ./ a; lw = 2, c = layer_color(j), m = :circle, ms = 2, label = layer_name(r, j))
        end
        plot!(p, r.t[its], yb ./ a; lw = 2.5, c = :black, ls = :dash, label = "ū")
        vline!(p, tw[1] .+ r.meta.T .* (0:floor(Int, (tw[2] - tw[1]) / r.meta.T)); c = :gray, lw = 0.4, label = "")
        p
    end
    plot(ps...; layout = (length(ps), 1), size = (1100, 240 * length(ps)), plot_titlefontsize = 11, _M...,
         plot_title = @sprintf("Harmonics seen at x = 0, t ∈ [%.1f, %.1f] s", r.t[its[1]], r.t[its[end]]))
end

"""
    plot_mode_cycle(r, n; t0=end − T, cycles=1) -> Plot

The vertical profile of harmonic n at every snapshot over `cycles` wave periods from `t0`, coloured
by time (the profile is taken at x = 0, as the real part of the complex amplitude). A steady mode
draws the same curve every time (only growing); an oscillating one swings through ± (a surface wave
keeps its shape and flips sign; a mixture changes shape). Right: the same at the x where the
harmonic is largest at t0.
"""
function plot_mode_cycle(r::PeriodicRun, n::Int; t0 = r.t[end] - r.meta.T, cycles = 1)
    its = twindow(r, (t0, t0 + cycles * r.meta.T)); N = nlayers(r)
    C = [xharmonic(r.U[j, its, :], n) for j in 1:N]           # C[j][i]
    σf = range(0, 1; length = 101); Φf = [phi(r.basis, j, σ) for σ in σf, j in 1:N]
    θmax = -angle(C[N][1])                                      # x where the surface value peaks at t0
    a = maximum(abs(C[j][i]) for j in 1:N, i in eachindex(its))
    cols = cgrad(:viridis, max(length(its), 2); categorical = true)
    ps = map((("at x = 0", 0.0), ("at the x of the surface maximum", θmax))) do (lb, θ)
        p = plot(; xlabel = "u / max", ylabel = "height (m)", title = "n = $n, " * lb, titlefontsize = 10,
                 xlims = (-1.1, 1.1), ylims = (0, r.meta.d), legend = false)
        vline!(p, [0]; c = :gray, lw = 0.5)
        for (i, it) in enumerate(its)
            u = [real(C[j][i] * cis(θ)) for j in 1:N] ./ a
            plot!(p, Φf * u, σf .* r.meta.d; c = cols[i], lw = 1.6)
            scatter!(p, u, _node_sigma(r) .* r.meta.d; c = cols[i], ms = 3, msw = 0)
        end
        p
    end
    plot(ps...; layout = (1, 2), size = (900, 480), plot_titlefontsize = 10, _M...,
         plot_title = @sprintf("Profile of harmonic %d at every snapshot, t ∈ [%.1f, %.1f] s (dark → light = early → late)", n, r.t[its[1]], r.t[its[end]]))
end

"""
    plot_mode_phasors(r, ns; tw) -> Plot

The complex amplitude û_j(t) of each layer of harmonic n in the complex plane, over the window `tw`,
each sample divided by the surface layer's |û| at that time (this removes the growth and keeps the
phase). A steady pattern sits at fixed points; a travelling or oscillating one goes round circles
(once per period of its frequency). Layers on opposite sides of 0 move against each other.
"""
function plot_mode_phasors(r::PeriodicRun, ns; tw = (r.t[end] - 2r.meta.T, r.t[end]), labels = string.(ns))
    its = twindow(r, tw); N = nlayers(r); ub = depth_average(r)
    ps = map(enumerate(zip(ns, labels))) do (k, (n, lb))
        C = [xharmonic(r.U[j, its, :], n) for j in 1:N]; s = abs.(C[N])
        p = plot(; aspect_ratio = :equal, xlims = (-1.2, 1.2), ylims = (-1.2, 1.2), xlabel = "Re", ylabel = "Im",
                 title = "n = $n: $lb", titlefontsize = 10, legend = k == 1 ? :bottomright : false, legendfontsize = 7)
        θ = range(0, 2π; length = 100); plot!(p, cos.(θ), sin.(θ); c = :gray, lw = 0.5, label = "")
        for j in 1:N
            z = C[j] ./ s
            plot!(p, real.(z), imag.(z); c = layer_color(j), lw = 1, label = layer_name(r, j))
            scatter!(p, [real(z[end])], [imag(z[end])]; c = layer_color(j), ms = 5, label = "")
        end
        zb = xharmonic(ub[its, :], n) ./ s
        plot!(p, real.(zb), imag.(zb); c = :black, lw = 1.5, ls = :dash, label = "ū")
        p
    end
    plot(ps...; layout = (1, length(ps)), size = (340 * length(ps), 400), plot_titlefontsize = 10, _M...,
         plot_title = @sprintf("Phasors û_j(t) / |û_surface(t)|, t ∈ [%.1f, %.1f] s (dot = last sample)", r.t[its[1]], r.t[its[end]]))
end

"Amplitude |û_n| of every layer, of ū and of η√(g/d) over the whole run, one panel per harmonic."
function plot_mode_amplitudes(r::PeriodicRun, ns; labels = string.(ns))
    ub = depth_average(r); c = sqrt(r.meta.g / r.meta.d); N = nlayers(r)
    ps = map(enumerate(zip(ns, labels))) do (k, (n, lb))
        p = plot(; yscale = :log10, xlabel = "t (s)", ylabel = k == 1 ? "|û| (m/s)" : "", title = "n = $n: $lb",
                 titlefontsize = 10, legend = k == 1 ? :bottomright : false, legendfontsize = 7)
        for j in 1:N
            plot!(p, r.t, abs.(xharmonic(r.U[j, :, :], n)) .+ 1e-16; lw = 1.5, c = layer_color(j), label = layer_name(r, j))
        end
        plot!(p, r.t, abs.(xharmonic(ub, n)) .+ 1e-16; lw = 2, c = :black, label = "ū")
        plot!(p, r.t, c .* abs.(xharmonic(r.η, n)) .+ 1e-16; lw = 1.2, c = :black, ls = :dot, label = "η·√(g/d)")
        p
    end
    plot(ps...; layout = (1, length(ps)), size = (330 * length(ps), 420), plot_titlefontsize = 11, _M...,
         plot_title = "Amplitude of each harmonic, per layer")
end

"Vertical profile of each harmonic in `ns` at time t, at the x where its surface value is largest."
function plot_mode_profiles(r::PeriodicRun, t::Real, ns; labels = string.(ns))
    it = tindex(r, t); N = nlayers(r); σf = range(0, 1; length = 101); Φf = [phi(r.basis, j, σ) for σ in σf, j in 1:N]
    ps = map(zip(ns, labels)) do (n, lb)
        C = [xharmonic(r.U[j, it:it, :], n)[1] for j in 1:N]; θ = -angle(C[N])
        u = real.(C .* cis(θ)); a = maximum(abs.(Φf * u)); ū = sum(r.D .* u)
        p = plot(Φf * u ./ a, σf .* r.meta.d; lw = 2.5, c = :black, label = "", xlims = (-1.15, 1.15),
                 xlabel = "u / max", ylabel = "height (m)", title = "n = $n: $lb", titlefontsize = 10, legend = :bottomleft)
        scatter!(p, u ./ a, _node_sigma(r) .* r.meta.d; c = [layer_color(j) for j in 1:N], ms = 6, label = "")
        vline!(p, [ū / a]; c = :red, ls = :dash, lw = 2, label = "depth average")
        vline!(p, [0]; c = :gray, lw = 0.5, label = "")
        p
    end
    plot(ps...; layout = (1, length(ps)), size = (320 * length(ps), 450), plot_titlefontsize = 11, _M...,
         plot_title = @sprintf("Vertical profile of each harmonic at t = %.1f s (dots = the layer unknowns)", r.t[it]))
end

# ---- space–time ----------------------------------------------------------------------------------------

"""
    plot_layers_hovmoller(r; tw, part=:disturbance, row_normalise=true) -> Plot

x–t maps of every layer, of ū and of η√(g/d) for a `part` over the window `tw`. With
`row_normalise`, every time row is divided by its largest value over all panels (relative sizes are
kept; growth is removed). Vertical stripes = standing still, slanted = travelling.
"""
function plot_layers_hovmoller(r::PeriodicRun; tw = (r.t[end] - 15, r.t[end]), part = :disturbance,
                               row_normalise::Bool = true, nmin::Int = 5, ncarrier::Int = 4)
    its = twindow(r, tw); N = nlayers(r); c = sqrt(r.meta.g / r.meta.d)
    F = [part_rows(r, r.U[j, :, :], its, part; nmin, ncarrier) for j in 1:N]
    push!(F, sum(r.D[j] .* F[j] for j in 1:N), c .* part_rows(r, r.η, its, part; nmin, ncarrier))
    names = vcat([layer_name(r, j) for j in 1:N], ["ū (depth average)", "η·√(g/d)"])
    sc = row_normalise ? maximum(reduce(hcat, [maximum(abs.(f); dims = 2) for f in F]); dims = 2) : maximum(maximum(abs.(f)) for f in F)
    ps = [heatmap(r.x, r.t[its], f ./ sc; c = SIGNED_CMAP, clims = (-1, 1), colorbar = false, title = nm, titlefontsize = 9,
                  xlabel = "x (m)", ylabel = k == 1 ? "t (s)" : "") for (k, (f, nm)) in enumerate(zip(F, names))]
    plot(ps...; layout = (1, length(ps)), size = (290 * length(ps), 500), plot_titlefontsize = 11, _M...,
         plot_title = @sprintf("%s in space and time, t ∈ [%.1f, %.1f] s%s", part_label(part; nmin, ncarrier),
                               r.t[its[1]], r.t[its[end]], row_normalise ? " (normalised per time row)" : ""))
end

"The disturbance of one layer split into its steady part (mean over a period) and the rest, as x–t maps."
function plot_steady_split(r::PeriodicRun; tw = (r.t[end] - 15, r.t[end]), layer::Int = nlayers(r), nmin::Int = 5)
    its = twindow(r, tw)
    F = [part_rows(r, r.U[layer, :, :], its, p; nmin) for p in (:disturbance, :steady, :oscillating)]
    sc = maximum(abs.(F[1]); dims = 2)
    ttl = ["disturbance (n ≥ $nmin)", "steady part (mean over one period)", "oscillating part (the rest)"]
    ps = [heatmap(r.x, r.t[its], f ./ sc; c = SIGNED_CMAP, clims = (-1, 1), colorbar = false, title = tt, titlefontsize = 10,
                  xlabel = "x (m)", ylabel = k == 1 ? "t (s)" : "") for (k, (f, tt)) in enumerate(zip(F, ttl))]
    plot(ps...; layout = (1, 3), size = (1250, 500), plot_titlefontsize = 11, _M...,
         plot_title = layer_name(r, layer) * @sprintf(", t ∈ [%.1f, %.1f] s (normalised per time row)", r.t[its[1]], r.t[its[end]]))
end

"""
    plot_column(r, t; parts=[:carrier, :steady, :oscillating], xw=(0, Lx), zw=(0, d)) -> Plot

The velocity inside the water column (colour) and streamlines for several parts at time t, side
by side; `xw`/`zw` zoom (the carrier is always drawn over the whole box).
"""
function plot_column(r::PeriodicRun, t::Real; parts = [:carrier, :steady, :oscillating], xw = (0.0, r.meta.Lx),
                     zw = (0.0, r.meta.d), nmin::Int = 5, ncarrier::Int = 4)
    it = tindex(r, t); d = r.meta.d
    ps = map(parts) do part
        _, U, _ = part_state(r, it, part; nmin, ncarrier)
        ix = part === :carrier ? eachindex(r.x) : findall(x -> xw[1] ≤ x ≤ xw[2], r.x)
        σs, u, ψ = column_field(r, U[:, ix]); z = σs .* d; iz = findall(zz -> zw[1] ≤ zz ≤ zw[2], z)
        a = maximum(abs.(u[iz, :])); b = maximum(abs.(ψ[iz, :]))
        p = heatmap(r.x[ix], z[iz], u[iz, :] ./ a; c = SIGNED_CMAP, clims = (-1, 1), colorbar = false, xlabel = "x (m)",
                    ylabel = "height above bed (m)", title = part_label(part; nmin, ncarrier) * @sprintf("\nmax |u| %.1e m/s", a),
                    titlefontsize = 9, aspect_ratio = part === :carrier ? :auto : :equal,
                    xlims = (r.x[ix[1]], r.x[ix[end]]), ylims = (z[iz[1]], z[iz[end]]), legend = false)
        b > 0 && contour!(p, r.x[ix], z[iz], ψ[iz, :]; levels = collect(range(-b, b; length = 16))[2:end-1], c = :black, lw = 0.7, colorbar = false)
        for s in _node_sigma(r)[2:end-1]
            hline!(p, [s * d]; c = :gray, ls = :dash, lw = 1)
        end
        p
    end
    plot(ps...; layout = (1, length(ps)), size = (480 * length(ps), 560), plot_titlefontsize = 10, _M..., left_margin = 9Plots.mm,
         plot_title = @sprintf("Inside the water column at t = %.1f s: colour u (red = +x), black streamlines, dashed layer interfaces", r.t[it]))
end

"""
    plot_komega(r; tw, layer=Nσ, nmax=50, marks=[]) -> Plot

Power of layer `layer` in wavenumber (box harmonic n) and frequency ω over the window `tw` (Hann in
time); ω > 0 travels towards +x. White: the model's surface-wave branch ±ω(k); red: ω = 0.
`marks = [(n, ω/ω0, "label"), …]` circles points (e.g. the members of a triad).
⚠ ω is folded into ±π/Δt_snap (±4ω0 for T/8 sampling).
"""
function plot_komega(r::PeriodicRun; tw = (r.t[1], r.t[end]), layer::Int = nlayers(r), nmax::Int = 50, marks = [])
    its = twindow(r, tw); nt = length(its); dts = r.t[2] - r.t[1]
    v = r.U[layer, its, :]; v = v .- mean(v; dims = 1); Ns = size(v, 2)
    F = fft(v .* (0.5 .- 0.5 .* cos.(2π .* (0:nt-1) ./ nt)), (1, 2))
    ωobs = -2π .* fftfreq(nt, 1 / dts); ps = sortperm(ωobs); ns = -min(nmax, Ns ÷ 2 - 1):min(nmax, Ns ÷ 2 - 1)
    P = abs2.(F[ps, mod.(ns, Ns) .+ 1])        # e^{i(k_n x − ωt)} sits at temporal bin −ω
    ωk, _ = model_dispersion(r); k1 = 2π / r.meta.Lx
    p = heatmap(ns, ωobs[ps] ./ ω0(r), log10.(P ./ maximum(P) .+ 1e-14); c = :viridis, clims = (-10, 0),
                xlabel = "box harmonic n (wavenumber = n × carrier wavenumber)", ylabel = "ω / ω0",
                title = @sprintf("Power of %s, t ∈ [%.1f, %.1f] s (colour: log10)", layer_name(r, layer), r.t[its[1]], r.t[its[end]]),
                titlefontsize = 10, size = (1050, 620), _M..., background_color_legend = :gray25,
                foreground_color_legend = nothing, legendfontcolor = :white)
    nn = range(first(ns), last(ns); length = 801)
    plot!(p, nn, [sign(n) * ωk(abs(n) * k1) for n in nn] ./ ω0(r); c = :white, lw = 1, ls = :dash, label = "surface-wave branch")
    plot!(p, nn, [-sign(n) * ωk(abs(n) * k1) for n in nn] ./ ω0(r); c = :white, lw = 1, ls = :dash, label = "")
    hline!(p, [0]; c = :red, lw = 0.6, label = "ω = 0 (steady)")
    for (n, w, lb) in marks
        θ = range(0, 2π; length = 60)
        plot!(p, n .+ 1.6 .* cos.(θ), w .+ 0.22 .* sin.(θ); c = :red, lw = 2, label = "")
        annotate!(p, n, w + 0.42, text(lb, 9, :red))
    end
    return p
end

# ---- energy --------------------------------------------------------------------------------------------

"""
    plot_energy(r; bands=default_bands(r), nmin=5) -> Plot

(a) total energy over E0: the solver's monitor (diagnostics.csv: ½gη² + ½H uᵀMu, which leaves out
the kinetic energy of w), the same from the snapshots linearised (H → d, a check), and the
snapshots' linear energy WITH ½d³(∂ₓu)ᵀB(∂ₓu), the energy the linear layered model conserves (its
dispersion is ω² = g d k² Dᵀ(A + k²d²B)⁻¹D); the axis is clipped to the pre-break-down range;
(b) energy per band of harmonics over time (log, including w);
(c) the disturbance (n ≥ nmin) energy split into surface ½gη², depth-mean kinetic ½dū², vertical-
shear kinetic ½d(uᵀAu − ū²) and vertical kinetic ½d³(∂ₓu)ᵀB(∂ₓu); (d) the same as fractions.
"""
function plot_energy(r::PeriodicRun; bands = default_bands(r), nmin::Int = 5)
    E = harmonic_energy(r); hor = E.pot .+ E.mean .+ E.shear; tot = hor .+ E.w
    es = vec(sum(hor; dims = 2)); ew = vec(sum(tot; dims = 2))
    curves = Any[(E.t, es ./ es[1], "snapshots, no w (linearised monitor)", :red, :dash), (E.t, ew ./ ew[1], "snapshots, with w (linear model energy)", :blue, :solid)]
    md = ""
    if r.diag !== nothing
        pushfirst!(curves, (r.diag["t"], r.diag["energy_ratio"], "solver monitor (no w)", :black, :solid))
        md = @sprintf(", max |mass drift| %.1e m^2", maximum(abs.(r.diag["mass_drift"])))
    end
    ok = [c[1] .≤ end_of_linear_growth(r; nmin) + 20 for c in curves]          # clip the final blow-up
    lo, hi = extrema(reduce(vcat, [c[2][k] for (c, k) in zip(curves, ok)])); pad = 0.1 * (hi - lo) + 1e-4
    pa = plot(; xlabel = "t (s)", ylabel = "E / E0", title = "(a) total energy" * md, titlefontsize = 10, legend = :bottomleft,
              ylims = (lo - pad, hi + pad), legendfontsize = 7)
    for (t, y, lb, c, ls) in curves
        plot!(pa, t, y; c = c, ls = ls, lw = 1.3, label = lb)
    end
    pb = plot(; yscale = :log10, xlabel = "t (s)", ylabel = "energy (m^3/s^2 per m^2)", title = "(b) energy per band of harmonics",
              titlefontsize = 10, legend = :outerright, legendfontsize = 7)
    for (k, (nm, rg)) in enumerate(bands)
        plot!(pb, E.t, vec(sum(tot[:, rg]; dims = 2)) .+ 1e-30; lw = 1.8, label = nm, c = k)
    end
    rg = nmin:nharmonics(r)
    parts = [("surface gη^2/2", E.pot), ("depth-mean kinetic dū^2/2", E.mean), ("vertical-shear kinetic", E.shear), ("vertical kinetic (w)", E.w)]
    cols = [:black, :royalblue, :darkorange, :seagreen]
    pc = plot(; yscale = :log10, xlabel = "t (s)", ylabel = "energy", title = "(c) disturbance (n ≥ $nmin): where its energy is",
              titlefontsize = 10, legend = :outerright, legendfontsize = 7)
    S = [vec(sum(e[:, rg]; dims = 2)) for (_, e) in parts]
    for ((nm, _), s, c) in zip(parts, S, cols)
        plot!(pc, E.t, s .+ 1e-30; lw = 1.8, c = c, label = nm)
    end
    ST = sum(S)
    pd = plot(; xlabel = "t (s)", ylabel = "fraction", ylims = (0, 1), title = "(d) the same, as fractions of the disturbance energy",
              titlefontsize = 10, legend = :outerright, legendfontsize = 7)
    for ((nm, _), s, c) in zip(parts, S, cols)
        plot!(pd, E.t, s ./ ST; lw = 1.8, c = c, label = nm)
    end
    plot(pa, pb, pc, pd; layout = (2, 2), size = (1500, 850), plot_titlefontsize = 11, _M...,
         plot_title = r.meta.name * ": energy")
end

"""
    plot_energy_history(r; modes=:groups, avg=:period, relative=false, tw=(t_first, t_last),
                        nmin=5, ncarrier=4) -> Plot

Energy against time (log scale), from `energy_history` (nodal modes: each mode once, no images):
* `modes = :all`: one line per mode n = 1 … ncell·p_u/2, coloured by n (colour bar), and n = 0 (the
  x-uniform flow) red dashed;
* `modes = :groups`: total, carrier, bound harmonics, disturbance and n = 0;
* `modes = [14, 16, 32]`: those modes.
`relative = true` divides by the total energy at t = 0. `avg = :period` smooths the exchange between the
surface and velocity energy within each wave period; `:instant` shows the snapshots.
"""
function plot_energy_history(r::PeriodicRun; modes = :groups, avg::Symbol = :period, relative::Bool = false,
                             tw = (r.t[1], r.t[end]), nmin::Int = 5, ncarrier::Int = 4, size = (1100, 600))
    H = energy_history(r; modes, avg, nmin, ncarrier); its = twindow(r, tw); t = H.t[its]
    y = H.E[its, :] ./ (relative ? H.E0 : 1.0)
    y[.!(y .> 0)] .= NaN                                        # log axis: drop zeros
    what = modes === :all ? "every mode" : modes === :groups ? "groups of modes" : "selected modes"
    p = plot(; yscale = :log10, xlabel = "t (s)", ylabel = relative ? "energy / total at t = 0" : "energy (m³/s²)",
             xlims = extrema(t), legend = modes === :all ? :bottomright : :outerright, legendfontsize = 8, titlefontsize = 11, size, _M...,
             title = @sprintf("%s: energy of %s (%s)", r.meta.name, what, avg === :period ? "period mean" : "instant"))
    if modes === :all
        nmax = last(H.n)
        for (k, n) in enumerate(H.n)
            n == 0 && continue
            plot!(p, t, y[:, k]; line_z = fill(n, length(t)), c = :turbo, clims = (1, nmax), lw = 1, label = "",
                  colorbar_title = "mode n")
        end
        plot!(p, t, y[:, 1]; c = :red, ls = :dash, lw = 1.5, label = "n = 0 (x-uniform flow)")
    elseif modes === :groups
        for (k, (c, lw, ls)) in enumerate(((:gray, 3.0, :solid), (:black, 2.0, :solid), (:darkorange, 1.8, :solid),
                                           (:crimson, 2.0, :solid), (:purple, 1.5, :dash)))
            plot!(p, t, y[:, k]; c, lw, ls, label = H.labels[k])
        end
    else
        for k in eachindex(H.n)
            plot!(p, t, y[:, k]; lw = 1.8, c = k, label = H.labels[k])
        end
    end
    return p
end

"""
    plot_energy_spectrum(r; times=5 evenly spaced, nmax=nh, split_at=last of times) -> Plot

(a) energy per harmonic n (surface + kinetic + vertical kinetic) at several times: how the spectrum
fills up; (b) at one time (`split_at`), the composition of each harmonic's energy.
"""
function plot_energy_spectrum(r::PeriodicRun; times = range(r.t[1], r.t[end]; length = 5), nmax::Int = nharmonics(r),
                              split_at = last(times))
    its = unique([tindex(r, t) for t in times]); E = harmonic_energy(r; its = its); n = 1:nmax
    cols = cgrad(:viridis, max(length(its), 2); categorical = true)
    pa = plot(; yscale = :log10, xlabel = "box harmonic n", ylabel = "energy", title = "(a) energy spectrum at several times",
              titlefontsize = 10, legend = :topright, legendfontsize = 7)
    for (i, it) in enumerate(its)
        plot!(pa, n, (E.pot .+ E.mean .+ E.shear .+ E.w)[i, n] .+ 1e-30; c = cols[i], lw = 1.6, m = :circle, ms = 2, msw = 0,
              label = @sprintf("t = %.0f s", r.t[it]))
    end
    vline!(pa, [r.meta.ncell]; c = :gray, ls = :dash, label = "1 wavelength per cell")
    i = argmin(abs.(r.t[its] .- split_at))
    pb = plot(; yscale = :log10, xlabel = "box harmonic n", ylabel = "energy", titlefontsize = 10, legend = :topright, legendfontsize = 7,
              title = @sprintf("(b) composition at t = %.0f s", r.t[its[i]]))
    for (nm, e, c) in (("surface gη^2/2", E.pot, :black), ("depth-mean kinetic", E.mean, :royalblue),
                       ("vertical-shear kinetic", E.shear, :darkorange), ("vertical kinetic (w)", E.w, :seagreen))
        plot!(pb, n, e[i, n] .+ 1e-30; c = c, lw = 1.6, m = :circle, ms = 2, msw = 0, label = nm)
    end
    plot(pa, pb; layout = (1, 2), size = (1450, 520), plot_titlefontsize = 11, _M..., plot_title = r.meta.name * ": energy per harmonic")
end

"Size (rms over the box) of the disturbance n ≥ nmin in every layer, in ū and in η√(g/d), over the run."
function plot_growth_overview(r::PeriodicRun; nmin::Int = 5)
    nh = nharmonics(r); ub = depth_average(r); c = sqrt(r.meta.g / r.meta.d)
    p = plot(; yscale = :log10, xlabel = "t (s)", ylabel = "rms over the box (m/s)", legend = :bottomright,
             title = r.meta.name * ": size of the disturbance (n ≥ $nmin) per layer", size = (950, 520), _M...)
    for j in 1:nlayers(r)
        plot!(p, r.t, _rms(xband(r.U[j, :, :], nmin:nh)) .+ 1e-30; lw = 2, c = layer_color(j), label = layer_name(r, j))
    end
    plot!(p, r.t, _rms(xband(ub, nmin:nh)) .+ 1e-30; lw = 2.5, c = :black, label = "ū (depth average)")
    plot!(p, r.t, c .* _rms(xband(r.η, nmin:nh)) .+ 1e-30; lw = 1.5, c = :black, ls = :dot, label = "η·√(g/d)")
    plot!(p, r.t, _rms(xband(ub, 1:1)); lw = 1.5, c = :gray, ls = :dash, label = "for scale: the carrier's ū (n = 1)")
    return p
end

# ---- everything at once --------------------------------------------------------------------------------

"""
    periodic_report(r; outdir="<run>/viz", animate=true, frames=select_frames(r; dt=1.0),
                    parts_anim=[:total, :disturbance], nmin=5, io=stdout) -> outdir

Every figure of this file for one run, with the windows, harmonics and triad chosen from the data
(`end_of_linear_growth`, `growing_harmonics`, `pick_triad`), plus the animations of `parts_anim`.
Prints the growth table and the choices made.
"""
function periodic_report(r::PeriodicRun; outdir = joinpath(r.dir, "viz"), animate::Bool = true,
                         frames = select_frames(r; dt = 1.0), parts_anim = [:total, :disturbance], nmin::Int = 5, io = stdout)
    mkpath(outdir); T = r.meta.T
    tl = end_of_linear_growth(r; nmin); tw = (max(r.t[1], tl - 60), tl)
    rows = growing_harmonics(r; tw, top = 12, io)
    growing_harmonics(r; tw, top = 8, io, sortby = :amp)      # the same table, largest first
    triad = pick_triad(rows; nmin)
    ns = triad === nothing ? [1, 2, rows[1].n, r.meta.ncell] : [1, 2, triad[1], triad[2], r.meta.ncell]
    lb = triad === nothing ? ["carrier", "2nd harmonic", "fastest growing", "1 wavelength per cell"] :
         ["carrier", "pump (2nd harmonic)", "steady member", "fast member", "1 wavelength per cell"]
    @printf(io, "end of the clean growth t = %.1f s; harmonics %s; triad %s\n", tl, string(ns), string(triad))
    sv(p, f) = (savefig(p, joinpath(outdir, f)); println(io, "  ", f))
    sv(plot_growth_overview(r; nmin), "01_growth_overview.png")
    sv(plot_energy(r; nmin), "02_energy.png")
    sv(plot_energy_spectrum(r), "03_energy_spectrum.png")
    sv(plot_energy_history(r; modes = :groups, relative = true, nmin), "15_energy_history_groups.png")
    sv(plot_energy_history(r; modes = :all), "16_energy_history_modes.png")
    sv(plot_mode_amplitudes(r, ns; labels = lb), "04_mode_amplitudes.png")
    sv(plot_mode_timeseries(r, ns; tw = (tl - 4T - 5, tl - 5), labels = lb), "05_mode_timeseries.png")
    sv(plot_mode_phasors(r, ns; tw = (tl - 2T - 5, tl - 5), labels = lb), "06_mode_phasors.png")
    for (k, n) in enumerate(ns[3:end])
        sv(plot_mode_cycle(r, n; t0 = tl - T - 5), @sprintf("07_mode_cycle_n%d.png", n))
    end
    sv(plot_mode_profiles(r, tl - 5, ns; labels = lb), "08_mode_profiles.png")
    sv(plot_state(r, tl - 5; part = :total), "09_state_total.png")
    sv(plot_state(r, tl - 5; part = :disturbance), "10_state_disturbance.png")
    sv(plot_column(r, tl - 5; xw = (r.meta.Lx / 2, 3r.meta.Lx / 4), zw = (0.55r.meta.d, r.meta.d), nmin), "11_column.png")
    sv(plot_layers_hovmoller(r; tw = (tl - 15, tl), nmin), "12_hovmoller_layers.png")
    sv(plot_steady_split(r; tw = (tl - 15, tl), nmin), "13_steady_split.png")
    marks = triad === nothing ? [] : [(2, 2.0, "pump"), (triad[1], 0.0, "steady"), (-triad[2], 2.0, "fast")]
    sv(plot_komega(r; tw, marks), "14_k_omega.png")
    if animate
        for part in parts_anim
            println(io, "  ", basename(animate_run(r; frames, part, out = joinpath(outdir, "state_$(part).mp4"), nmin)))
        end
    end
    return outdir
end
