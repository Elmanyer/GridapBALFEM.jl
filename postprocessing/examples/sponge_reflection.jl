# ==============================================================
#  sponge_reflection.jl — reflection coefficient of the outflow sponge, per case of the sweep
#
#  Companion of run/local/run_1d_sponge_sweep.sh (CH8_CAMPAIGN_TODO.md §7). For every finished case
#  of the stage's table (<sweep>/_logs/_cases_<stage>.tsv):
#    1. the complex amplitude of η at ω, a(x) = (2/N) Σ_j η(x, t_j) e^{iωt_j}, over the LAST TEN
#       PERIODS (N = 10T/Δt_save samples: a whole number of periods, so no leakage from the
#       carrier; Δt_save = 0.1 s divides every T of the sweep);
#    2. the least-squares fit a(x) ≈ a_I e^{ikx} + a_R e^{−ikx} over the fixed 30 m stretch
#       [x_fit0, x_fit1] of the interior (the fit of LINEAR_PRODUCTION_RESULTS.md). k is fitted too,
#       by a scan of ±5 % around Airy (the discrete model's k differs from Airy by ≲ 0.5 %, and a
#       wrong k reads as reflection);
#    3. Kr = |a_R| / |a_I|, the relative rms residual of the fit, the generated amplitude
#       |a_I| / A, and |a| at the last node (what reaches the wall), relative to |a_I|.
#  Sign: with the e^{+iωt} kernel a wave cos(kx − ωt) travelling in +x maps to e^{+ikx}.
#
#  ⚠ The residual is part of the result. A residual comparable to Kr means the two-wave model does
#  not describe the stretch (unsteady window, a third wave, ...), and that Kr is not a measurement.
#  ⚠ Floor cases (W = 48, μ = 4) measure what the generation and the relaxation zone reflect on
#  their own: a matrix Kr is only meaningful well above the floor at the same T.
#
#  USAGE
#    julia --project=postprocessing postprocessing/examples/sponge_reflection.jl <sweep_dir> [A|B]
#  Writes <sweep_dir>/kr_<stage>.csv and kr_<stage>_T<T>.png (heat map of Kr(μ_max, W) per T) and
#  prints the tables.
# ==============================================================

include(joinpath(@__DIR__, "..", "src", "GridapBALFEMPost.jl"))
using .GridapBALFEMPost
using Printf, Statistics, LinearAlgebra
import .GridapBALFEMPost: Plots

length(ARGS) >= 1 || error("usage: sponge_reflection.jl <sweep_dir> [A|B]")
SWEEP = ARGS[1]
STAGE = length(ARGS) >= 2 ? ARGS[2] : "A"
const G, D, A_WAVE = 9.81, 3.5, 0.1

function airy_k(T)
    w2 = (2π / T)^2 / G
    k = max(w2, sqrt(w2 / D))
    for _ in 1:100
        k -= (k * tanh(k * D) - w2) / (tanh(k * D) + k * D / cosh(k * D)^2)
    end
    return k
end

function read_cases(path)
    lines = filter(!isempty, readlines(path))
    hdr = split(lines[1], '\t')
    return [Dict(zip(hdr, split(l, '\t'))) for l in lines[2:end]]
end

"η on the unique x nodes (y-averaged) for the last `nper` periods: (x, t, H[nx × nt])."
function last_periods(rundir, T, nper)
    pvd = filter(f -> endswith(f, ".pvd"), readdir(rundir; join = true))
    isempty(pvd) && error("no .pvd in $rundir")
    vtus, times = GridapBALFEMPost._parse_pvd(first(pvd))
    q = median(diff(times[max(1, end - 20):end]))
    N = round(Int, nper * T / q)
    abs(N * q - nper * T) < 1e-6 * T || error("save interval $q does not divide $nper T = $(nper*T)")
    length(vtus) >= N || error("only $(length(vtus)) snapshots, need $N")
    sel = (length(vtus) - N + 1):length(vtus)
    P, F = load_snapshot(vtus[sel[1]])
    xr  = round.(P[:, 1]; digits = 8)
    xu  = sort(unique(xr))
    idx = [findall(==(x), xr) for x in xu]
    H = zeros(length(xu), N)
    for (j, it) in enumerate(sel)
        _, Fj = load_snapshot(vtus[it])
        η = Fj["eta"]
        H[:, j] .= [mean(η[ii]) for ii in idx]
    end
    return xu, times[sel], H
end

"Two-wave fit at fixed k: (a_I, a_R, relative rms residual)."
function fit2(x, a, k)
    M = hcat(exp.(im * k .* x), exp.(-im * k .* x))
    c = M \ a
    r = a - M * c
    return c[1], c[2], norm(r) / sqrt(length(a)) / abs(c[1])
end

function analyse(c)
    T = parse(Float64, c["T"])
    rundir = joinpath(SWEEP, c["name"])
    x, t, H = last_periods(rundir, T, 10)
    ω = 2π / T
    a = (2 / length(t)) .* (H * exp.(im * ω .* t))
    x0, x1 = parse(Float64, c["x_fit0"]), parse(Float64, c["x_fit1"])
    m = (x .>= x0) .& (x .<= x1)
    k0 = airy_k(T)
    best = (Inf, k0, 0.0im, 0.0im)
    for k in range(0.95k0, 1.05k0; length = 2001)
        aI, aR, res = fit2(x[m], a[m], k)
        res < best[1] && (best = (res, k, aI, aR))
    end
    res, k, aI, aR = best
    return (name = c["name"], T = T, W = parse(Float64, c["W"]), mu = parse(Float64, c["mu"]),
            relax_w = parse(Float64, c["relax_w"]), role = c["role"],
            Kr = abs(aR) / abs(aI), resid = res, k_fit = k, k_airy = k0,
            gen = abs(aI) / A_WAVE, a_end = abs(a[end]) / abs(aI), t_end = t[end])
end

cases = read_cases(joinpath(SWEEP, "_logs", "_cases_$STAGE.tsv"))
done(c) = (log = joinpath(SWEEP, "_logs", c["name"] * ".log");
           isfile(log) && occursin(r"flume_1d \[.*\] done:", read(log, String)))
R = []
for c in cases
    if !done(c)
        @printf("  skip %-28s (not finished)\n", c["name"])
        continue
    end
    try
        push!(R, analyse(c))
    catch e
        @printf("  FAIL %-28s %s\n", c["name"], sprint(showerror, e))
    end
end
isempty(R) && error("no finished case in $SWEEP (stage $STAGE)")

csv = joinpath(SWEEP, "kr_$STAGE.csv")
open(csv, "w") do io
    println(io, "name,role,T,W,mu,relax_w,Kr,resid,k_fit,k_airy,k_err,gen_amp,a_end,t_end")
    for r in R
        @printf(io, "%s,%s,%g,%g,%g,%g,%.5e,%.3e,%.6f,%.6f,%.3e,%.4f,%.3e,%.2f\n",
                r.name, r.role, r.T, r.W, r.mu, r.relax_w, r.Kr, r.resid, r.k_fit, r.k_airy,
                r.k_fit / r.k_airy - 1, r.gen, r.a_end, r.t_end)
    end
end
println("\nwrote $csv  ($(length(R)) of $(length(cases)) cases)\n")

for T in sort(unique(r.T for r in R))
    RT = filter(r -> r.T == T, R)
    floor = filter(r -> r.role == "floor", RT)
    @printf("T = %.1f s   (k_fit/k_Airy − 1 = %+.2f…%+.2f %%; floor Kr = %s)\n", T,
            100 * minimum(r.k_fit / r.k_airy - 1 for r in RT),
            100 * maximum(r.k_fit / r.k_airy - 1 for r in RT),
            isempty(floor) ? "—" : @sprintf("%.2f %% (resid %.1e)", 100 * floor[1].Kr, floor[1].resid))
    Ws  = sort(unique(r.W for r in RT if r.role != "floor"))
    mus = sort(unique(r.mu for r in RT if r.role != "floor"))
    isempty(Ws) && continue
    @printf("  Kr [%%]   %s\n", join([@sprintf("W=%-6g", W) for W in Ws], " "))
    K = fill(NaN, length(mus), length(Ws))
    for (i, μ) in enumerate(mus)
        row = String[]
        for (j, W) in enumerate(Ws)
            hit = filter(r -> r.W == W && r.mu == μ && r.role != "floor", RT)
            if isempty(hit)
                push!(row, "   —    ")
            else
                K[i, j] = 100 * hit[1].Kr
                push!(row, @sprintf("%6.2f%s ", K[i, j], hit[1].resid > 0.5 * hit[1].Kr ? "?" : " "))
            end
        end
        @printf("  μ=%-5g  %s\n", μ, join(row, " "))
    end
    println("  (? = fit residual > Kr/2: not a measurement)\n")
    if STAGE == "A"
        p = Plots.heatmap(string.(Int.(Ws)), string.(Int.(mus)), K;
                          xlabel = "sponge width W [m]", ylabel = "μ_max [1/s]",
                          title = @sprintf("Kr [%%] at T = %.1f s", T), c = :viridis,
                          clims = (0, max(1.0, maximum(filter(!isnan, K)))))
        Plots.savefig(p, joinpath(SWEEP, @sprintf("kr_%s_T%.1f.png", STAGE, T)))
    end
end

rep = filter(r -> r.role == "replica", R)
if !isempty(rep)
    println("Round-1 replica (W = 12, μ = 40) against Snellius (5 % at 2.0 s, 9 % at 2.2 s, 18 % at 2.9 s):")
    for r in sort(rep; by = r -> r.T)
        @printf("  T = %.1f s: Kr = %.2f %%\n", r.T, 100 * r.Kr)
    end
end
