# ==============================================================
#  run_skew_mms.jl — V5/V6 of the skew-advection verification ladder
#                    (building_files/SKEW_SYMMETRIC_ADVECTION_PLAN.md §3)
#
#  THE QUESTION. The energy-consistent correction adds ½∫𝒞(𝗠·U)·W with 𝒞 the
#  POINTWISE strong residual of the continuity row. 𝒞 ≡ 0 for the exact solution,
#  so the scheme stays consistent and the MMS forcing is unchanged — but 𝒞 evaluated
#  on the DISCRETE solution is O(h^{p_h}) through ∇·(Hū), one order below the
#  velocity error itself. The classical skew/Temam analysis retains optimal rates,
#  but that is an ENERGY-norm result and these are L² rates.
#
#  ⚠ SO A DROP OF p_u FROM 4 TO 3 IS A PLAUSIBLE OUTCOME, NOT A BUG. Decide by the
#  rule in the plan §3, and read the pairwise SEQUENCE, never the fitted slope
#  (CLAUDE.md rule 33): a rising sequence with the error still falling is
#  PRE-ASYMPTOTIC; a flat sequence with a stalled error is not.
#
#  ⚠ EVERY MMS NUMBER RECORDED BEFORE THIS BRANCH IS VOID FOR skew=true. This
#  script therefore runs each model TWICE — off and on — in the same process, so the
#  comparison is in-harness and one-variable rather than against a remembered value.
#  The study tag carries "SKEW", so the CSV cannot silently mix them.
#
#  Models covered: the four NONLINEAR :none/:native models (3,4,5,6). The linear
#  models 1-2 assemble no advection block and are unaffected by construction —
#  running them would report a false confirmation, not a check. Models 7-8 (:full)
#  are excluded here: their floor is an amplitude-sensitive quasi-Newton limit
#  (CLAUDE.md §5) and mixing that with a residual change measures neither.
#
#  RUN:  julia --project=. examples/local_mms/run_skew_mms.jl
#        BALFEM_SKEW_MMS_MODELS=M3,M5 julia --project=. examples/local_mms/run_skew_mms.jl
#        BALFEM_SKEW_MMS_MODE=transient  (V6; default is the spatial ladder)
# ==============================================================

using GridapBALFEM
using Printf, DelimitedFiles

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
genv(k, d) = get(ENV, k, d)

const MODE   = Symbol(genv("BALFEM_SKEW_MMS_MODE", "static"))
const LEVELS = parse(Int, genv("BALFEM_SKEW_MMS_LEVELS", MODE === :static ? "4" : "4"))
const NX0    = parse(Int, genv("BALFEM_SKEW_MMS_NX0", "8"))
const P_U    = parse(Int, genv("BALFEM_SKEW_MMS_PU", "3"))     # Q3/Q2 — the verified pairing

all_models = [
    (tag="M3", regime=:nonlinear, nlp=:none,   flat=true ),
    (tag="M4", regime=:nonlinear, nlp=:none,   flat=false),
    (tag="M5", regime=:nonlinear, nlp=:native, flat=true ),
    (tag="M6", regime=:nonlinear, nlp=:native, flat=false),
]
SEL = genv("BALFEM_SKEW_MMS_MODELS", "all")
models = SEL == "all" ? all_models :
         let want = strip.(split(SEL, ","))
             bad = setdiff(want, [m.tag for m in all_models])
             isempty(bad) || error("run_skew_mms: unknown model tag(s): $(join(bad, ", "))")
             [m for m in all_models if m.tag in want]
         end

outdir = genv("BALFEM_OUTDIR", joinpath(ROOT, "output", "local", "skew_mms"))
mkpath(outdir)

println("=" ^ 78)
println("  run_skew_mms.jl — MMS convergence, skew-advection OFF vs ON")
@printf("  mode=%s  Q%d/Q%d  levels=%d (nx0=%d)  models=%s\n",
        MODE, P_U, P_U-1, LEVELS, NX0, SEL)
println("=" ^ 78)
flush(stdout)

rows = Any[]
for m in models, skew in (false, true)
    @printf("\n>>> %s  %s / %s / %s   skew=%s\n", m.tag, m.regime,
            m.flat ? "flat" : "varbed", m.nlp, skew)
    flush(stdout)
    r = run_conv_study(; p_u = P_U, domain = :d1, mode = MODE, levels = LEVELS,
                         nx0 = NX0, M = 2, p_vert = 1,
                         regime = m.regime, nl_pressure = m.nlp,
                         flat_bed = m.flat, a_b = m.flat ? 0.0 : 0.2,
                         skew_advection = skew,
                         #  a_eta at the campaign default: these are :none/:native
                         #  models, whose quasi-Newton gap is not the amplitude cliff
                         #  the :full pair hits.
                         nl_iter = 200, nl_tol = 1e-14)
    push!(rows, (m.tag, skew, r.tag, r.fit_eta, r.opt_eta, r.fit_u, r.opt_u,
                 r.pw_eta[end], r.pw_u[end], r.e_eta[end], r.e_u[end],
                 join(round.(r.pw_eta, digits=3), " "),
                 join(round.(r.pw_u,   digits=3), " ")))
    flush(stdout)
end

csv = joinpath(outdir, "skew_mms_$(MODE)_Q$(P_U).csv")
open(csv, "w") do io
    println(io, "model,skew,tag,fit_eta,opt_eta,fit_u,opt_u,pw_eta_last,pw_u_last,",
                "e_eta_fine,e_u_fine,pw_eta_seq,pw_u_seq")
    for r in rows
        println(io, join(r, ","))
    end
end

println()
println("=" ^ 78)
println("  SUMMARY — the comparison that matters is WITHIN each model, across skew")
println("=" ^ 78)
@printf("  %-5s %-6s  %-9s %-9s  %-11s %-11s\n",
        "model", "skew", "pw_eta", "pw_u", "e_eta(fine)", "e_u(fine)")
for r in rows
    @printf("  %-5s %-6s  %9.3f %9.3f  %11.3e %11.3e\n", r[1], r[2], r[8], r[9], r[10], r[11])
end
println()
println("  Full pairwise sequences (read these, not the fitted slopes):")
for r in rows
    @printf("    %-5s skew=%-5s  eta: [%s]   u: [%s]\n", r[1], r[2], r[12], r[13])
end
println("\n  wrote $csv")
