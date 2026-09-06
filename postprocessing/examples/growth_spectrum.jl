# ==============================================================
#  growth_spectrum.jl — rank an instability by what GROWS, never by what is LARGEST
#
#  WHY THIS EXISTS AS A SCRIPT. The 2026-09-05 diagnosis of the nonlinear grid-scale
#  instability (building_files/NONLINEAR_INSTABILITY.md §5b) was done ad hoc, and its
#  FIRST attempt returned the wrong answer convincingly: ranking the spectrum by
#  AMPLITUDE above an arbitrary cut found the carrier's own harmonics at λ ≈ 1.5–1.8 m,
#  a wavelength that does not scale with dx — which reads as evidence AGAINST a grid
#  mode. Ranking the same data by GAIN found the real mode at λ ≈ 2–3·dx.
#
#  An unstable mode is small in absolute terms for most of its life. That is exactly
#  why it goes unnoticed until it dominates.
#
#  Two guards are built in because both were violated in that analysis:
#    * NO band cut. The full resolved range is reported, binned. A cut at 0.95·k_Nyq
#      once placed a reported "peak" exactly on the window edge — and when an extremum
#      sits at the edge of the search window, the WINDOW is the finding.
#    * The carrier and its harmonics are reported SEPARATELY from the gain ranking,
#      so physical steepening is never mistaken for the instability (or damped by a
#      "cure" that flattens both).
#
#  USAGE
#    julia --project=postprocessing postprocessing/examples/growth_spectrum.jl \
#          <run_dir> [t_early] [t_late]
#
#    Compare two runs (e.g. two element pairings, or two meshes) by running it on
#    each and reading the `gain` column — that column, not eta_max, is what a cure
#    has to change.
# ==============================================================

using GridapBALFEMPost
using FFTW, Printf, Statistics

length(ARGS) >= 1 || error("usage: growth_spectrum.jl <run_dir> [t_early] [t_late]")
const RUN = ARGS[1]

sim = load_simulation(RUN; fields = ["eta"])
times = sim.times
length(times) >= 2 || error("growth_spectrum: need ≥2 snapshots in $RUN")

#  Extract η along the flume centreline: take the unique x-nodes and average over y
#  (a 1-D case is y-invariant, so this is exact; in 2-D it is a transect average and
#  should be read as such).
pts = sim.points
xs  = pts[:, 1]
xu  = sort(unique(round.(xs; digits = 9)))
idx = [findall(abs.(xs .- x) .< 1e-9) for x in xu]
etam = sim.fields["eta"].data                     # [n_points × n_times]
line(it) = [mean(etam[ii, it]) for ii in idx]

#  Snapshot selection: default to the first and last available, but prefer an EARLY
#  time after the ramp and a LATE time before any abort.
t_early = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : times[max(1, cld(length(times), 8))]
t_late  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : times[end]
ie = argmin(abs.(times .- t_early));  il = argmin(abs.(times .- t_late))
ie == il && error("growth_spectrum: early and late snapshots coincide (t=$(times[ie]))")

dx    = length(xu) > 1 ? xu[2] - xu[1] : NaN
n     = length(xu)
k_nyq = π / dx

#  Detrend and window before the FFT: a non-periodic transect otherwise leaks the
#  carrier across the whole spectrum, which is precisely the band the mode lives in.
function spec(v)
    v = v .- mean(v)
    w = 0.5 .* (1 .- cos.(2π .* (0:n-1) ./ (n - 1)))       # Hann
    F = abs.(rfft(v .* w))
    k = 2π .* (0:length(F)-1) ./ (n * dx)
    return k, F
end

k, Fe = spec(line(ie))
_, Fl = spec(line(il))

@printf("\n  run        : %s\n", RUN)
@printf("  x-nodes    : %d   dx = %.4f m   k_Nyq = %.3f rad/m\n", n, dx, k_nyq)
@printf("  snapshots  : t_early = %.3f s  →  t_late = %.3f s  (Δt = %.3f s)\n\n",
        times[ie], times[il], times[il] - times[ie])

#  ---- gain, binned over the FULL resolved range (no cut) ----------------------
#  Binning matters: a per-mode gain is noisy, and the mode of interest is broad.
const NBIN = 12
edges = range(0.0, k_nyq; length = NBIN + 1)
println("  k/k_Nyq      λ/dx       A_early      A_late        GAIN")
println("  " * "-"^58)
best_gain = -Inf; best_bin = 0
for b in 1:NBIN
    sel = findall(x -> edges[b] <= x < edges[b+1], k)
    isempty(sel) && continue
    ae = sqrt(sum(Fe[sel].^2));  al = sqrt(sum(Fl[sel].^2))
    ae < 1e-14 && continue
    gn = al / ae
    kc = 0.5*(edges[b] + edges[b+1])
    @printf("  %6.3f     %6.2f    %10.3e   %10.3e   %10.1f×\n",
            kc/k_nyq, 2π/kc/dx, ae, al, gn)
    gn > best_gain && (best_gain = gn; best_bin = b)
end

kc = 0.5*(edges[best_bin] + edges[best_bin+1])
@printf("\n  PEAK GAIN  %.1f×  at k/k_Nyq = %.3f  (λ ≈ %.2f·dx = %.3f m)\n",
        best_gain, kc/k_nyq, 2π/kc/dx, 2π/kc)
if best_bin == NBIN || best_bin == 1
    println("  ⚠ THE PEAK SITS ON THE EDGE OF THE SEARCH WINDOW — the window is the")
    println("    finding, not the peak. Re-read with a different snapshot pair before")
    println("    quoting this number.")
end

#  ---- the carrier and its harmonics, reported SEPARATELY ----------------------
#  These grow by ordinary nonlinear steepening and are PHYSICAL. A cure that flattens
#  them along with the grid-scale mode has damped the physics, and this is where that
#  shows up.
kc_car = k[argmax(Fe)]                      # carrier from the EARLY spectrum
@printf("\n  carrier k_c = %.3f rad/m (λ = %.2f m).  Harmonic bands n·k_c ± 0.4:\n",
        kc_car, 2π/kc_car)
println("     n     A_early      A_late        gain")
for nh in 1:5
    sel = findall(x -> abs(x - nh*kc_car) < 0.4, k)
    isempty(sel) && continue
    ae = sqrt(sum(Fe[sel].^2));  al = sqrt(sum(Fl[sel].^2))
    @printf("    %2d   %10.3e   %10.3e   %8.1f×\n", nh, ae, al, ae > 0 ? al/ae : NaN)
end
println()
