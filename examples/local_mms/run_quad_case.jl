# ==============================================================================
#  run_quad_case.jl — ONE quadrature-probe study, so the probe can be run in
#  parallel processes instead of one sequential script.
#
#  USAGE: julia --project=. examples/local_mms/run_quad_case.jl <model> <quad_extra>
#
#  Parameters are copied EXACTLY from the Campaign C job (run_phaseB_shard.jl):
#  p_u=3, :d1, levels=5, nx0=4, M=2, p_vert=1, a_eta=0.8, nl_tol default 1e-14.
#  ONLY `quad_extra` and `model` vary — everything else is pinned so a difference
#  between rows can only come from those two axes (rule 38c).
#
#  WHY A DOSE-RESPONSE AND NOT A SINGLE TREATMENT. The nonlinear advection
#  integrand H*(u.grad u)*phi is degree 10 per direction at Q3/Q2; Gridap's
#  default Measure(8) is exact to 9 (measured), quad_extra=2 -> exact to 11,
#  quad_extra=4 -> exact to 13. If under-integration is the cause, the rate must
#  recover ALREADY at q=2 (the first degree that covers the integrand) and then
#  STAY recovered at q=4. A single treatment cannot show that shape, and a
#  monotone dose-response is much harder to fake than one point.
#
#  MODEL KEY (matches MODELS in run_phaseB_shard.jl):
#     3 = nonlinear / flat     / :none      <- where the degradation was measured
#     4 = nonlinear / VARIABLE / :none      <- does any fix survive a sloping bed?
# ==============================================================================
using GridapBALFEM, Printf, Dates

const MODEL = parse(Int, ARGS[1])
const QEX   = parse(Int, ARGS[2])

MODELS = Dict(
 1=>(regime=:linear,    flat_bed=true,  nlp=:none),
 2=>(regime=:linear,    flat_bed=false, nlp=:none),
 3=>(regime=:nonlinear, flat_bed=true,  nlp=:none),
 4=>(regime=:nonlinear, flat_bed=false, nlp=:none),
 5=>(regime=:nonlinear, flat_bed=true,  nlp=:native),
 6=>(regime=:nonlinear, flat_bed=false, nlp=:native))
mo = MODELS[MODEL]

OUT = joinpath(@__DIR__, "..", "..", "output", "local", "mms", "quadrature_test")
mkpath(OUT)
CSV = joinpath(OUT, @sprintf("result_m%d_q%d.csv", MODEL, QEX))

@printf("[quad] model %d (%s / %s / %s)  quad_extra=%d   start %s\n",
        MODEL, mo.regime, mo.flat_bed ? "flat" : "varbed", mo.nlp, QEX,
        Dates.format(now(), "HH:MM:SS")); flush(stdout)

t0 = time()
r = run_conv_study(p_u = 3, domain = :d1, mode = :static, levels = 5, nx0 = 4,
                   M = 2, p_vert = 1, a_eta = 0.8,
                   regime = mo.regime, flat_bed = mo.flat_bed, nl_pressure = mo.nlp,
                   a_b = mo.flat_bed ? 0.0 : 0.2,
                   quad_extra = QEX,
                   output_dir = joinpath(OUT, @sprintf("m%d_q%d", MODEL, QEX)),
                   diag_every = 5, diag_csv = true,
                   monitor_factory = () -> SolverMonitor(),
                   verbose = true)
el = time() - t0

ratios = [r.e_eta[i]/r.e_eta[i+1] for i in 1:length(r.e_eta)-1]
open(CSV, "w") do io
    println(io, "model,regime,flat_bed,nl_pressure,quad_extra,p_eta_last,p_u_last,",
                "fit_eta,fit_u,e_eta_list,e_u_list,pw_eta_list,pw_u_list,ratio_list,seconds")
    @printf(io, "%d,%s,%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%s,%s,%s,%s,%s,%.1f\n",
            MODEL, mo.regime, mo.flat_bed, mo.nlp, QEX,
            r.pw_eta[end], r.pw_u[end], r.fit_eta, r.fit_u,
            join([@sprintf("%.6e",x) for x in r.e_eta], ";"),
            join([@sprintf("%.6e",x) for x in r.e_u],   ";"),
            join([@sprintf("%.4f",x) for x in r.pw_eta], ";"),
            join([@sprintf("%.4f",x) for x in r.pw_u],   ";"),
            join([@sprintf("%.3f",x) for x in ratios],  ";"), el)
end
@printf("[quad] model %d q=%d DONE in %.1f min | p_eta(last)=%.4f (opt 3) | ratios %s\n",
        MODEL, QEX, el/60, r.pw_eta[end],
        join([@sprintf("%.2f",x) for x in ratios], " "))
flush(stdout)
