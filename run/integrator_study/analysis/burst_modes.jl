# burst_modes.jl — what grows just before each burst onset (2026-11-23)
# For a case with an onset (events.tsv), window = [onset − 15, onset + 2] periods: per box harmonic,
# the amplitude growth rate of the VS part |v_n|_M and the observed frequency of û_n(u₁) (phase fit;
# ω_obs aliased into (−4ω₀, 4ω₀] by the T/8 sampling, ω₀ = 2π/T of the case). Top 4 by rate.
# Usage: after analyse.jl in the same session: include(".../burst_modes.jl")
ev = readdlm(joinpath(AN_DIR, "events.tsv"), '\t', String; header = true)
eh = vec(ev[2]); erows = [Dict(eh .=> vec(ev[1][i, :])) for i in 1:size(ev[1], 1)]
unwrap(p) = (q = copy(p); for i in 2:length(q); q[i] -= 2π * round((q[i] - q[i-1]) / 2π); end; q)
open(joinpath(AN_DIR, "burst_modes.log"), "w") do io
    for e in erows
        e["onset"] in ("none", "NaN") && continue
        m = only(filter(x -> x["name"] == e["name"], rows))
        T = parse(Float64, m["T"]); ω0 = 2π / T
        on = parse(Float64, e["onset"])
        r = deserialize(joinpath(AN_DIR, "_coeffs", e["name"] * ".jls"))
        per = r.t ./ T; iw = findall(p -> on - 15 ≤ p ≤ on + 2, per)
        length(iw) < 20 && continue
        nh = size(r.C["eta"], 2); tt = r.t[iw]
        res = map(2:nh) do n
            v = an_vs(r, n)[iw]; σ = an_rate(tt, v)
            ph = unwrap(angle.(r.C["u1"][iw, n])); ω = -(cov(tt, ph) / var(tt))
            ηr = abs(r.C["eta"][iw[end], n]) / max(v[end], 1e-300)
            (n = n, σ = σ, ω = ω / ω0, amp = v[end], ηr = ηr)
        end
        sort!(res; by = q -> isnan(q.σ) ? Inf : -q.σ)
        @printf(io, "%-36s %-6s onset P%-4s  ", e["name"], e["outcome"], e["onset"])
        for q in res[1:4]
            @printf(io, " n%-3d σ%+.3f ω/ω₀%+.2f |v|%.0e η/v %.1e |", q.n, q.σ, q.ω, q.amp, q.ηr)
        end
        println(io)
    end
end
println("burst_modes.log written")
