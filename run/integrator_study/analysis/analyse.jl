# analyse.jl — growth analysis of the Snellius integrator study from the box harmonics (2026-11-23)
#
# THE SPLIT. At box harmonic n (k_n = 2πn/L) the linear P1LFE-2 operator, (A + k²d²B̃) u̇ = −ikg D η,
# has one surface-wave velocity profile s(k) = M⁻¹D (M = A + k²d²B̃) and Nσ − 1 = 2 steady
# vertical-structure (VS) profiles with Dᵀu = 0 (zero depth-average, ū = 0). They are M-orthogonal
# (sᵀ M v = Dᵀ v = 0), so every velocity vector splits uniquely:
#     u = α s + v,   α = Dᵀu / DᵀM⁻¹D  (∝ ū),   v = the VS part,   |v|_M = √(vᴴ M v)  [m/s]
# |v|_M is the VS amplitude (kinetic-energy norm, vertical velocity included). Its growth is what the
# 1-D campaign found under every scheme (CLAUDE.md §5: "every growing content is a vertical-structure
# velocity mode"); no integrator damps it at ω = 0 (rule 15).
#
# PER CASE (sampled at fixed numbers of periods so cases of different T compare):
#   carrier    |η̂₁|(t)/|η̂₁|(1 period)        (amplitude kept; rule 15's carrier cost)
#   VS         √Σ_n |v_n|²_M over bands: low n = 1–4 (n = 1 included: the carrier's own profile
#              distortion), grid = 5 … ncell·p_u/2 (node spectrum), sub = above the node Nyquist
#   surf       √Σ |η̂_n|² for n ≥ 5 (short surface waves: the burst signature)
#   σ          amplitude growth rate (1/s) of each, log-linear over the last 25 % of the record
#              (≥ 10 periods), and the same relative to the carrier
#
# Usage (julia-mcp, postprocessing env, repository root as cwd), after extract.jl:
#   include("run/integrator_study/analysis/analyse.jl")     # reads/writes $IS_OUT (summary.tsv first)
# Writes metrics.tsv (one row per case) and events.tsv in $IS_OUT.
using Serialization, Printf, Statistics, LinearAlgebra, DelimitedFiles

const IS_ROOT2 = normpath(joinpath(@__DIR__, "..", "..", ".."))
const AN_DIR = get(ENV, "IS_OUT", joinpath(IS_ROOT2, "output", "snellius", "integrator_study_analysis"))
const AN_g = 9.81; const AN_d = 3.5

"A, B̃, D of P1LFE-2 (hat basis on σ ∈ [0, 0.728, 1]); as temporal_spectrum.jl ts_linear_branch."
function an_vertical()
    s = [0.0, 0.728, 1.0]; ng = 4000
    σ = range(0, 1; length = ng + 1); hσ = 1 / ng
    φ = [[clamp(j == 1 ? (s[2] - x) / s[2] : j == 3 ? (x - s[2]) / (1 - s[2]) :
                 (x ≤ s[2] ? x / s[2] : (1 - x) / (1 - s[2])), 0, 1) for x in σ] for j in 1:3]
    w = [vcat(0.0, cumsum((f[1:end-1] .+ f[2:end]) ./ 2) .* hσ) for f in φ]
    q(a, b) = hσ * (sum(a .* b) - (a[1] * b[1] + a[end] * b[end]) / 2)
    A = [q(φ[i], φ[j]) for i in 1:3, j in 1:3]; B = [q(w[i], w[j]) for i in 1:3, j in 1:3]
    D = [q(φ[i], ones(ng + 1)) for i in 1:3]
    return A, B, D
end
const AN_A, AN_B, AN_Dv = an_vertical()

"VS amplitude |v_n|_M (per sample) of harmonic n, from the layer coefficients u₁…u₃."
function an_vs(r, n)
    k = 2π * n / r.L
    M = AN_A + k^2 * AN_d^2 * AN_B
    s = M \ AN_Dv; den = dot(AN_Dv, s)
    out = zeros(length(r.t))
    for i in eachindex(r.t)
        u = ComplexF64[r.C["u1"][i, n], r.C["u2"][i, n], r.C["u3"][i, n]]
        v = u .- (dot(AN_Dv, u) / den) .* s
        out[i] = sqrt(max(real(dot(v, M * v)), 0.0))
    end
    return out
end

an_rate(t, y) = (ok = (y .> 0) .& isfinite.(y); sum(ok) < 3 ? NaN :
                 (tt = t[ok]; ly = log.(y[ok]); cov(tt, ly) / var(tt)))

function an_case(case, meta)
    r = deserialize(joinpath(AN_DIR, "_coeffs", case * ".jls"))
    T = parse(Float64, meta["T"]); ncell = parse(Int, meta["ncell"])
    nh = size(r.C["eta"], 2); nnyq = min(3ncell ÷ 2, nh)           # node Nyquist (p_u = 3, one λ)
    vs = hcat((an_vs(r, n) for n in 1:nh)...)                       # nt × nh
    band(c) = vec(sqrt.(sum(abs2, c; dims = 2)))
    B = Dict("vs_low" => band(vs[:, 1:4]), "vs_grid" => band(vs[:, 5:nnyq]),
             "vs_sub" => nnyq < nh ? band(vs[:, nnyq+1:nh]) : zeros(length(r.t)),
             "vs_all" => band(vs), "surf" => band(r.C["eta"][:, 5:nh]),
             "car" => abs.(r.C["eta"][:, 1]))
    t = r.t; per = t ./ T
    i1 = argmin(abs.(per .- 1.0)); car0 = B["car"][i1]
    at(key, P) = (P > per[end] + 1e-9) ? NaN : B[key][argmin(abs.(per .- P))]
    # last-25 % window (at least 10 periods)
    w0 = max(per[end] - max(0.25per[end], 10.0), per[1])
    iw = findall(p -> p ≥ w0, per)
    σ = Dict(k => an_rate(t[iw], B[k][iw]) for k in keys(B))
    # fastest single VS harmonic over the window
    rates = [an_rate(t[iw], vs[iw, n]) for n in 2:nh]
    nfast = argmax(replace(rates, NaN => -Inf)) + 1
    return (r = r, B = B, per = per, car0 = car0, at = at, σ = σ, w0 = w0, nfast = nfast,
            σfast = rates[nfast - 1], vsfast_end = vs[end, nfast], nnyq = nnyq)
end

# ---- driver -------------------------------------------------------------------------------------
sm = readdlm(joinpath(AN_DIR, "summary.tsv"), '\t', String; header = true)
hdr = vec(sm[2]); rows = [Dict(hdr .=> vec(sm[1][i, :])) for i in 1:size(sm[1], 1)]
const AN_P = (1, 10, 25, 50, 75, 100, 125, 150, 175, 200)
open(joinpath(AN_DIR, "metrics.tsv"), "w") do io
    cols = ["name", "outcome", "per_end", "car_kept_end", "sig_car",
            ["vs_all_P$P" for P in AN_P]..., ["vs_grid_P$P" for P in AN_P]..., ["surf_P$P" for P in AN_P]...,
            "win_from_P", "sig_vs_low", "sig_vs_grid", "sig_vs_sub", "sig_vs_all", "sig_surf",
            "n_fast", "sig_fast", "vs_fast_end", "vs_all_over_car_end"]
    println(io, join(cols, '\t'))
    for m in rows
        isfile(joinpath(AN_DIR, "_coeffs", m["name"] * ".jls")) || continue
        a = an_case(m["name"], m)
        vals = Any[m["name"], m["outcome"], a.per[end], a.B["car"][end] / a.car0, a.σ["car"],
                   [a.at("vs_all", P) for P in AN_P]..., [a.at("vs_grid", P) for P in AN_P]...,
                   [a.at("surf", P) for P in AN_P]..., a.w0,
                   a.σ["vs_low"], a.σ["vs_grid"], a.σ["vs_sub"], a.σ["vs_all"], a.σ["surf"],
                   a.nfast, a.σfast, a.vsfast_end, a.B["vs_all"][end] / a.B["car"][end]]
        println(io, join([v isa AbstractFloat ? @sprintf("%.4g", v) : string(v) for v in vals], '\t'))
    end
end
println("metrics.tsv written")

# ---- events: burst onset and pre-burst growth, from period-rms series ----------------------------
# Per band, the rms over each period (removes the carrier-phase modulation). Baseline = median of
# periods 10–40 (or what exists). ONSET = first period whose rms exceeds 10× the baseline.
# Pre-onset growth: amplitude rate (1/s) of the period-rms over [P20, onset − 5 periods or end].
function an_events(case, meta)
    a = an_case(case, meta)
    T = parse(Float64, meta["T"])
    np = floor(Int, a.per[end])
    vs2 = Dict("vs_low2" => nothing)
    r = a.r; nh = size(r.C["eta"], 2)
    low2 = vec(sqrt.(sum(abs2, hcat((an_vs(r, n) for n in 2:4)...); dims = 2)))
    series = Dict("vs_low2" => low2, "vs_grid" => a.B["vs_grid"], "vs_sub" => a.B["vs_sub"],
                  "surf" => a.B["surf"], "car" => a.B["car"])
    prms(y) = [sqrt(mean(abs2, y[(a.per .> p - 1) .& (a.per .≤ p)])) for p in 1:np]
    out = Dict{String,Any}()
    for (k, y) in series
        s = prms(y); out[k] = s
    end
    onset = Dict{String,Float64}()
    for k in ("vs_grid", "vs_sub", "surf", "vs_low2")
        s = out[k]; length(s) < 12 && (onset[k] = NaN; continue)
        base = median(s[10:min(40, end)])
        i = findfirst(>(10base), s[11:end]); onset[k] = i === nothing ? Inf : i + 10
    end
    on = minimum(values(onset))
    p1 = isfinite(on) ? Int(on) - 5 : np
    rate(k) = p1 - 20 < 10 ? NaN : an_rate(collect(20:p1) .* T, out[k][20:p1])
    return (np = np, onset = onset, on = on, pre = Dict(k => rate(k) for k in keys(series)),
            car_end = out["car"][end] / out["car"][1])
end

open(joinpath(AN_DIR, "events.tsv"), "w") do io
    println(io, join(["name", "outcome", "per_end", "onset_grid", "onset_sub", "onset_surf", "onset_low2",
                      "onset", "pre_car", "pre_vs_low2", "pre_vs_grid", "pre_vs_sub", "pre_surf"], '\t'))
    for m in rows
        isfile(joinpath(AN_DIR, "_coeffs", m["name"] * ".jls")) || continue
        e = try an_events(m["name"], m) catch err; @warn m["name"] err; continue end
        f(x) = x isa AbstractFloat ? (isinf(x) ? "none" : @sprintf("%.4g", x)) : string(x)
        println(io, join(f.(Any[m["name"], m["outcome"], e.np, e.onset["vs_grid"], e.onset["vs_sub"],
                               e.onset["surf"], e.onset["vs_low2"], e.on, e.pre["car"], e.pre["vs_low2"],
                               e.pre["vs_grid"], e.pre["vs_sub"], e.pre["surf"]]), '\t'))
    end
end
println("events.tsv written")
