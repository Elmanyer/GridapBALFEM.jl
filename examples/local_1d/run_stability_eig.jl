# ==============================================================
#  run_stability_eig.jl — frozen-state eigen-analysis batch driver (method: stability_eig.jl)
#
#  Models `:off` (nl_pressure=false) and `:on` (all eight 𝓝 components, Class III by the broken
#  formulation), unstabilised and with the ghost penalty.
#
#  Runs as a detached process, not through julia-mcp: the nonlinear-pressure residual takes long
#  to compile, beyond the MCP tool's 30-min idle limit. One compile then serves the whole ladder.
#
#  OUTPUT (persistent disk, rule 47):  output/local_1d/stability_eig/
#     summary.csv           one row per case
#     eig_<tag>.csv         every eigenvalue: re, im, kdom, uy_fraction
#
#  RUN:  julia --project=. examples/local_1d/run_stability_eig.jl
# ==============================================================
using GridapBALFEM, LinearAlgebra, Printf
LinearAlgebra.BLAS.set_num_threads(parse(Int, get(ENV, "SE_BLAS", "4")))
include(joinpath(@__DIR__, "stability_eig.jl"))

const OUT = joinpath(@__DIR__, "..", "..", "output", "local_1d", "stability_eig")
mkpath(OUT)
const SUMMARY = joinpath(OUT, "summary.csv")
isfile(SUMMARY) || open(SUMMARY, "w") do io
    println(io, "model,pu,pe,gamma,stab,A,ncell,dx,ndof,sigma_max,omega_at_max,k_at_max_over_k0,",
                "lambda_over_dx_at_max,uy_at_max,rho,m_gate,ncolour,wall_s")
end

tagof(c) = @sprintf("%s_Q%dQ%d_g%g_A%.2f_n%d", c.model, c.pu, c.pe, c.gamma, c.A, c.ncell)
done_tags() = isfile(SUMMARY) ? Set(begin f = split(l, ','); @sprintf("%s_Q%sQ%s_g%g_A%.2f_n%s",
                  f[1], f[2], f[3], parse(Float64, f[4]), parse(Float64, f[6]), f[7]) end
                  for l in readlines(SUMMARY)[2:end]) : Set{String}()

case(model, pu, A, n; gamma = 0.0) = (model = model, pu = pu, pe = pu - 1, gamma = gamma, A = A, ncell = n)

cases = [
    # sanity: at rest every eigenvalue must be imaginary (no sponge, non-dissipative core)
    case(:on, 3, 0.0, 8),
    # THE LADDER: A = 0.10 (Reg03), dx refinement, both pairings, 𝓝 off vs on
    [case(m, 2, 0.10, n) for m in (:off, :on) for n in (8, 16, 32)]...,
    [case(m, 3, 0.10, n) for m in (:off, :on) for n in (8, 16, 32)]...,
    # amplitude
    [case(:on, 3, 0.15, n) for n in (8, 16, 32)]...,
    # the stabilised operator (ghost penalty, γ* = 0.01 and its window edges)
    [case(:on, 3, 0.10, n; gamma = g) for g in (0.0033, 0.01, 0.03) for n in (16, 32)]...,
]

println("stability_eig batch: $(length(cases)) cases -> $OUT"); flush(stdout)
done = done_tags()
for c in cases
    tag = tagof(c)
    tag in done && (println("skip (done) $tag"); continue)
    t0 = time()
    r = try
        stability_case(; pu = c.pu, pe = c.pe, model = c.model, A = c.A, ncell = c.ncell,
                         gamma = c.gamma)
    catch e
        println("FAILED $tag: ", sprint(showerror, e)); flush(stdout); continue
    end
    wall = time() - t0
    println(se_line(r), @sprintf("  [%.0f s]", wall)); flush(stdout)
    open(SUMMARY, "a") do io
        @printf(io, "%s,%d,%d,%g,%s,%.4f,%d,%.6f,%d,%.8e,%.6e,%.4f,%.4f,%.4f,%.6e,%.3e,%d,%.1f\n",
                r.model, r.pu, r.pe, r.gamma, r.stab, r.A, r.ncell, r.dx, r.ndof, r.σmax, r.ωmax,
                r.kmax / r.k0, r.kmax > 0 ? (2π / r.kmax) / r.dx : Inf, r.uymax, r.rho,
                r.m_gate, r.ncolour, wall)
    end
    open(joinpath(OUT, "eig_$tag.csv"), "w") do io
        println(io, "re,im,kdom,uy_fraction")
        for m in eachindex(r.λ)
            @printf(io, "%.10e,%.10e,%.8e,%.6f\n", real(r.λ[m]), imag(r.λ[m]), r.kdom[m], r.uyfr[m])
        end
    end
end
println("stability_eig batch: finished"); flush(stdout)
