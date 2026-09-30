# ==============================================================
#  run_broken_mms_check.jl — does the C⁰-IP penalty preserve the MMS orders? (T8)
#
#  markdown_files/BROKEN_FORMULATION_PLAN.md §T8. The penalty J_h is CONSISTENT (it vanishes
#  on the exact, smooth solution), so the analytic MMS spatial orders must be unchanged by it:
#  Q3/Q2 ⇒ p_η = 3, p_u = 4. Same ladder as test_mms_convergence.jl G6 (Lx=1.7, Ly=1.1, d=1,
#  nx = 6,12,24, dt = 2e-4, T = 0.02), for the linear model (Model 1) and the :native
#  nonlinear model (Model 5, flat bed), each at γ = 0 (control, same batch) and at the
#  CONSTRAINT-regime γ_u = γ_η = 0.3 used by the screening batch.
#
#  RUN:  julia --project=. examples/local_mms/run_broken_mms_check.jl
#  OUT:  output/local/broken_mms/mms_cip.csv (+ stdout)
# ==============================================================
using GridapBALFEM, Printf

const OUT = joinpath(@__DIR__, "..", "..", "output", "local", "broken_mms")
mkpath(OUT)
const CASE = (Lx = 1.7, Ly = 1.1, d = 1.0, M = 2, p_vert = 1, p_u = 3, p_eta = 2)
const GAMMA = parse(Float64, get(ENV, "MMS_CIP_GAMMA", "0.3"))

rows = String[]
for (lbl, regime, nlp) in (("linear", :linear, :none), ("native", :nonlinear, :native))
    for γ in (0.0, GAMMA)
        @printf("\n=== %s  γ_u = γ_η = %g ===\n", lbl, γ); flush(stdout)
        sp = run_mms_refinement(:space; levels = 3, nx0 = 6, ny0 = 4, dt0 = 2e-4, T_final = 0.02,
                                regime = regime, nl_pressure = nlp, CASE...,
                                cip_gamma_u = γ, cip_gamma_eta = γ)
        @printf("  → p_eta = %.3f   p_u = %.3f\n", sp.p_eta, sp.p_u); flush(stdout)
        push!(rows, @sprintf("%s,%g,%.6e,%.6e,%.4f,%.4f", lbl, γ, sp.e_eta[end], sp.e_u[end],
                             sp.p_eta, sp.p_u))
    end
end
open(joinpath(OUT, "mms_cip.csv"), "w") do io
    println(io, "model,gamma,e_eta_fine,e_u_fine,p_eta,p_u")
    foreach(r -> println(io, r), rows)
end
println("\nwritten ", joinpath(OUT, "mms_cip.csv"))
