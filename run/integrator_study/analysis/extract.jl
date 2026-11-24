# extract.jl — box harmonics of every case of the Snellius integrator study (2026-11-23)
#
# For each case: read the VTK snapshots of all segments in time order (later segments win on a
# repeated time; an unreadable file — the quota cut can truncate the last one — is skipped and
# counted), and store the complex box harmonics û_n(t), n = 1 … Ns/2, of η, ū, u₁, u₂, u₃
# (coefficient of e^{+ik_n x}; y-averaged; the periodic duplicate at x = L dropped), as in
# output/1D_stability_campaign/_logs/temporal_spectrum.jl. Cache: _coeffs/<case>.jls.
#
# Usage (julia-mcp, postprocessing env, repository root as cwd):
#   include("run/integrator_study/analysis/extract.jl")     # env IS_STUDY, IS_OUT (see below)
isdefined(Main, :GridapBALFEMPost) || include(joinpath(@__DIR__, "..", "..", "..", "postprocessing", "src", "GridapBALFEMPost.jl"))
#  IS_STUDY: the downloaded study tree; IS_OUT: where the caches and tables go (defaults: the
#  Snellius download of 2026-11-22 and its analysis folder, both under the gitignored output/)
const IS_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const IS_STUDY = get(ENV, "IS_STUDY", joinpath(IS_ROOT, "output", "snellius", "integrator_study"))
const IS_OUT = get(ENV, "IS_OUT", joinpath(IS_ROOT, "output", "snellius", "integrator_study_analysis"))
using .GridapBALFEMPost, FFTW, Printf, Statistics, Serialization

const EX_STUDY = IS_STUDY
const EX_CACHE = joinpath(IS_OUT, "_coeffs")
const EX_D = [0.364, 0.5, 0.136]                    # P1LFE-2 depth weights, c_bdy [0, 0.728, 1]

function ex_case(case)
    cache = joinpath(EX_CACHE, case * ".jls")
    isfile(cache) && return deserialize(cache)
    files = Dict{Float64,String}()
    for seg in sort(filter(s -> startswith(s, "seg"), readdir(joinpath(EX_STUDY, case))))
        dir = joinpath(EX_STUDY, case, seg)
        for f in readdir(dir)
            occursin(r"^sol_t_.*\.vtu$", f) || continue
            files[GridapBALFEMPost._time_from_name(f)] = joinpath(dir, f)
        end
    end
    ts = sort(collect(keys(files)))
    isempty(ts) && return nothing
    P, _ = GridapBALFEMPost.load_snapshot(files[ts[1]])
    xr = round.(P[:, 1]; digits = 8); xs = sort(unique(xr))[1:end-1]
    idx = [findall(==(x), xr) for x in xs]; Ns = length(xs); nh = Ns ÷ 2
    C = Dict(k => zeros(ComplexF64, length(ts), nh) for k in ("eta", "ubar", "u1", "u2", "u3"))
    keep = trues(length(ts)); bad = 0
    for (it, t) in enumerate(ts)
        F = try
            last(GridapBALFEMPost.load_snapshot(files[t]))
        catch
            keep[it] = false; bad += 1; continue
        end
        g(k) = [mean(F[k][i]) for i in idx]
        u = [g("u$(j)x") for j in 1:3]
        for (k, v) in ("eta" => g("eta"), "ubar" => sum(EX_D[j] .* u[j] for j in 1:3),
                       "u1" => u[1], "u2" => u[2], "u3" => u[3])
            C[k][it, :] = fft(v)[2:nh+1] ./ Ns
        end
    end
    r = (t = ts[keep], C = Dict(k => v[keep, :] for (k, v) in C), L = xs[end] + xs[2] - xs[1],
         Ns = Ns, bad = bad)
    serialize(cache, r)
    return r
end

mkpath(EX_CACHE)
cases = sort(filter(c -> startswith(c, "m") && isdir(joinpath(EX_STUDY, c)), readdir(EX_STUDY)))
t0 = time()
for (i, c) in enumerate(cases)
    r = ex_case(c)
    @printf("%3d/%d %-40s %5d snapshots  t_end %7.1f  bad %d  (%.0f s)\n", i, length(cases), c,
            r === nothing ? 0 : length(r.t), r === nothing ? NaN : r.t[end], r === nothing ? 0 : r.bad,
            time() - t0)
    flush(stdout)
end
