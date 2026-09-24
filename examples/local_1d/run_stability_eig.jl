# ==============================================================
#  run_stability_eig.jl — OPTION C batch driver (see stability_eig.jl for the method)
#
#  Runs as a detached process, not through julia-mcp: the mixed :full residual alone takes
#  >28 min to compile, beyond the MCP tool's 30-min idle limit (two sessions were killed on
#  2026-09-24). One compile then serves the whole ladder.
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
    println(io, "model,pu,pe,paux,A,ncell,dx,ndof,sigma_max,omega_at_max,k_at_max_over_k0,",
                "lambda_over_dx_at_max,uy_at_max,rho,aux_res,mix_leak,m_gate,ncolour,wall_s")
end

tagof(c) = @sprintf("%s_Q%dQ%d_aux%d_A%.2f_n%d", c.model, c.pu, c.pe, c.paux, c.A, c.ncell)
done_tags() = isfile(SUMMARY) ? Set(begin f = split(l, ','); @sprintf("%s_Q%sQ%s_aux%s_A%.2f_n%s",
                  f[1], f[2], f[3], f[4], parse(Float64, f[5]), f[6]) end
                  for l in readlines(SUMMARY)[2:end]) : Set{String}()

case(model, pu, A, n; paux = pu) = (model = model, pu = pu, pe = pu - 1, paux = paux, A = A, ncell = n)

cases = [
    # sanity: the mixed layout must also be exactly non-dissipative at rest
    case(:mixed, 2, 0.0, 8),
    # THE MAIN LADDER: A = 0.10 (Reg03), dx refinement, both pairings, :native vs :full-mixed
    [case(m, 2, 0.10, n) for m in (:native, :mixed) for n in (8, 16, 32, 64)]...,
    [case(m, 3, 0.10, n) for m in (:native, :mixed) for n in (8, 16, 32)]...,
    # amplitude: the A = 0.15 mixed run died where A = 0.10 held
    [case(:mixed, 2, 0.15, n) for n in (8, 16, 32, 64)]...,
    [case(:mixed, 3, 0.15, n) for n in (8, 16, 32)]...,
    # test A counterpart: aux one order below u (both runs died EARLIER than aux = p_u)
    [case(:mixed, 2, 0.10, n; paux = 1) for n in (16, 32)]...,
    [case(:mixed, 3, 0.10, n; paux = 2) for n in (16, 32)]...,
    # :none reference at A = 0.10 (advection only, no N)
    [case(:none, 2, 0.10, n) for n in (16, 32)]...,
    # finest levels last (largest dense eigenproblems)
    case(:native, 3, 0.10, 64), case(:mixed, 3, 0.10, 64), case(:mixed, 2, 0.10, 128),
]

println("stability_eig batch: $(length(cases)) cases -> $OUT"); flush(stdout)
done = done_tags()
for c in cases
    tag = tagof(c)
    tag in done && (println("skip (done) $tag"); continue)
    t0 = time()
    r = try
        stability_case(; pu = c.pu, pe = c.pe, paux = c.paux, model = c.model, A = c.A, ncell = c.ncell)
    catch e
        println("FAILED $tag: ", sprint(showerror, e)); flush(stdout); continue
    end
    wall = time() - t0
    println(se_line(r), @sprintf("  [%.0f s]", wall)); flush(stdout)
    open(SUMMARY, "a") do io
        @printf(io, "%s,%d,%d,%d,%.4f,%d,%.6f,%d,%.8e,%.6e,%.4f,%.4f,%.4f,%.6e,%.3e,%.3e,%.3e,%d,%.1f\n",
                r.model, r.pu, r.pe, r.paux, r.A, r.ncell, r.dx, r.ndof, r.σmax, r.ωmax,
                r.kmax / r.k0, r.kmax > 0 ? (2π / r.kmax) / r.dx : Inf, r.uymax, r.rho,
                r.aux_res, r.mix_leak, r.m_gate, r.ncolour, wall)
    end
    open(joinpath(OUT, "eig_$tag.csv"), "w") do io
        println(io, "re,im,kdom,uy_fraction")
        for m in eachindex(r.λ)
            @printf(io, "%.10e,%.10e,%.8e,%.6f\n", real(r.λ[m]), imag(r.λ[m]), r.kdom[m], r.uyfr[m])
        end
    end
end
println("stability_eig batch: finished"); flush(stdout)
