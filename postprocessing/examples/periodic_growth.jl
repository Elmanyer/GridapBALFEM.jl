# ==============================================================
#  periodic_growth.jl — growth rates of the wavenumber bands of an x-PERIODIC run
#
#  Companion of examples/local_1d/run_periodic_1d.jl. On an x-periodic box of exactly
#  n_λ wavelengths the Fourier transform along x is EXACT (no windowing, no leakage), so
#  the energy of every box harmonic can be tracked through time and a growth rate fitted
#  to each band. That is the Floquet-type rate the frozen eigen-analysis could not give:
#  it includes the carrier's full orbital cycle, so compression/extension phases average
#  out exactly as they do in the real dynamics.
#
#  ⚠ READ THE RIGHT BAND. The initial state is the LINEAR eigenmode, so the carrier sheds
#  bounded free harmonics at low n; that is physics, not instability. The stability
#  verdict is the HIGH band (upper half of the node spectrum) and, separately, the
#  sub-element band above the node Nyquist (C⁰ kink content). Rank by GAIN, never by
#  amplitude (CLAUDE.md rule 38f).
#
#  ⚠ NEEDS SUB-CELL VTK. Vertex-only snapshots alias every sub-element mode away; the
#  periodic driver writes BALFEM_VTK_NSUB = 2·p_u samples per cell by default.
#
#  USAGE
#    julia --project=postprocessing postprocessing/examples/periodic_growth.jl <run_dir> \
#          [t_fit_start] [t_fit_end]
#  Writes <run_dir>/band_energy.csv and prints one σ per band (energy growth rate, 1/s;
#  amplitude grows at σ/2).
# ==============================================================

include(joinpath(@__DIR__, "..", "src", "GridapBALFEMPost.jl"))
using .GridapBALFEMPost
using FFTW, Printf, Statistics

length(ARGS) >= 1 || error("usage: periodic_growth.jl <run_dir> [t_fit_start] [t_fit_end]")
RUN = ARGS[1]

sim   = load_simulation(RUN; fields = :all, regularize = false)
times = sim.times
nt    = length(times)
nt >= 3 || error("periodic_growth: need ≥ 3 snapshots in $RUN")

#  the surface-layer x-velocity: highest vertical node index among u<k>x
ukeys = sort(filter(k -> occursin(r"^u\d+x$", k), collect(keys(sim.fields))))
utop  = ukeys[argmax(parse.(Int, [m.match for m in match.(r"\d+", ukeys)]))]

#  unique x samples along the flume (y-invariant case: average over y), dropping the
#  periodic duplicate at x = Lx
xs  = sim.points[:, 1]
tol = 1e-9 * max(1.0, maximum(abs, xs))
xu  = sort(unique(round.(xs; digits = 9)))
Lx  = xu[end] - xu[1]
xu  = xu[1:end-1]                                   # x = Lx ≡ x = 0
idx = [findall(abs.(xs .- x) .< 1e-8) for x in xu]
N   = length(xu)
dxs = Lx / N
maximum(abs, diff(xu) .- dxs) < 1e-6 * dxs ||
    error("periodic_growth: samples are not uniform in x (nonuniform sub-cell output?)")
line(F, it) = [mean(F[ii, it]) for ii in idx]

#  harmonic index n ↔ k = 2πn/Lx; spectral energy |F̂_n|² (one-sided, n = 0 … N/2)
spec(v) = abs2.(rfft(v) ./ N)

#  bands, in BOX harmonics. The node Nyquist is read from the run (cells/λ and p_u are
#  passed as env so the script needs no solver import); defaults match the driver.
nlam  = parse(Int, get(ENV, "PG_NLAMBDA", "1"))
ncell = parse(Int, get(ENV, "PG_NCELL",   "16"))
pu    = parse(Int, get(ENV, "PG_PU",      "2"))
n_nodeNyq = pu * ncell * nlam ÷ 2                   # node spacing λ/(p·ncell)
n_top     = N ÷ 2
bands = [("carrier (n=nλ)",            nlam:nlam),
         ("low harmonics (2–4)·k0",    (2nlam):(4nlam)),
         ("mid band",                  (4nlam + 1):(n_nodeNyq ÷ 2)),             # may be empty on coarse meshes
         ("HIGH band (node spectrum)", (max(4nlam, n_nodeNyq ÷ 2) + 1):n_nodeNyq),  # disjoint from the bands above
         ("sub-element (> node Nyq)",  (n_nodeNyq + 1):n_top)]

E = Dict(f => zeros(length(bands), nt) for f in ("eta", utop))
for it in 1:nt, f in ("eta", utop)
    S = spec(line(sim.fields[f], it))
    for (b, (_, r)) in enumerate(bands)
        rr = intersect(r, 1:(length(S) - 1)) .+ 1   # harmonic n sits at index n+1
        E[f][b, it] = isempty(rr) ? NaN : sum(S[rr])
    end
end

open(joinpath(RUN, "band_energy.csv"), "w") do io
    println(io, "t," * join(["eta_b$b" for b in eachindex(bands)], ",") * "," *
                join(["u_b$b" for b in eachindex(bands)], ","))
    for it in 1:nt
        println(io, join([times[it]; E["eta"][:, it]; E[utop][:, it]], ","))
    end
end

#  least-squares slope of ln E over the fit window
t0 = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : times[max(1, cld(nt, 5))]
t1 = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : times[end]
w  = findall(t -> t0 <= t <= t1, times)
length(w) >= 3 || error("periodic_growth: fewer than 3 snapshots in [$t0, $t1]")
function slope(y)
    ok = [i for i in w if isfinite(y[i]) && y[i] > 0]
    length(ok) < 3 && return NaN
    tt = times[ok]; ly = log.(y[ok])
    tm = mean(tt); return sum((tt .- tm) .* (ly .- mean(ly))) / sum((tt .- tm) .^ 2)
end

@printf("\nperiodic_growth: %s\n", RUN)
@printf("  %d samples along x (Δx_sample = %.4f m), %d snapshots, fit window t ∈ [%.2f, %.2f] s\n",
        N, dxs, nt, times[w[1]], times[w[end]])
@printf("  bands in box harmonics (k = 2πn/Lx); node Nyquist n = %d; surface velocity field %s\n\n",
        n_nodeNyq, utop)
@printf("  %-28s %-12s %12s %12s %12s %12s\n", "band", "n range", "σ_E(η) 1/s", "σ_E(u) 1/s",
        "E(η) end/start", "E(u) end/start")
for (b, (name, r)) in enumerate(bands)
    ge = E["eta"][b, w[end]] / E["eta"][b, w[1]]
    gu = E[utop][b, w[end]] / E[utop][b, w[1]]
    @printf("  %-28s %-12s %+12.4e %+12.4e %12.3e %12.3e\n", name,
            "$(first(r))–$(last(r))", slope(E["eta"][b, :]), slope(E[utop][b, :]), ge, gu)
end
println("\n  Reading: σ_E ≈ 0 (|σ_E| ≪ 0.1/s) in the HIGH and sub-element bands ⇒ interior stable at")
println("  this resolution. σ_E > 0 there, rising on the dx ladder ⇒ grid-scale instability.")
