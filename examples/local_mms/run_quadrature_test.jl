# ==============================================================================
#  run_quadrature_test.jl — is the nonlinear order reduction a QUADRATURE crime?
#
#  THE OBSERVATION (Campaign C, 2026-09-11, 1-D Q3/Q2, a_eta=0.8, nl_tol=1e-14):
#      linear      e_eta ratios  7.95  7.99  8.00  8.00   -> optimal 3rd order
#      nonlinear   e_eta ratios  7.84  7.58  6.76  5.46   -> decaying toward 2nd
#  The error still FALLS, so this is an order reduction, not an algebraic floor.
#
#  ALREADY RULED OUT, from data in hand:
#    * algebraic contamination — Newton reaches ~1e-15 at EVERY level (rule 32)
#    * the 𝓝 pressure blocks   — :none vs :native agree to 4 digits
#    * bed-slope code          — flat vs variable bed differ by ~4e-4
#    * the linear core         — linear is textbook on the same basis and ladder
#  So it lives in the nonlinear (advection) block.
#
#  THE HYPOTHESIS. The MMS path integrates at degree 2*max(p_u,p_eta)+2 = 8 for
#  Q3/Q2, which in Gridap is a 5-point Gauss rule EXACT TO DEGREE 9 (measured).
#  The linear core H*u*grad(phi) (H = d, constant) is degree <=8 and exact; the
#  nonlinear advection H*(u.grad u)*phi with H = d+eta (Q2) is degree 10 and is
#  NOT — under-integrated by exactly one degree. An under-integrated term is a variational crime contributing a FIXED
#  lower-order consistency error — invisible while the discretisation error
#  dominates, emerging as it decays. That is exactly the measured shape.
#
#  THE DESIGN, and why it is two studies and not one:
#    * CONTROL (quad_extra=0) RUNS IN THIS SAME BATCH, not against the recorded
#      Campaign C numbers (rule 38c). It must reproduce 7.84/7.58/6.76/5.46. If it
#      does, the quad_extra edit to mms_driver.jl is clean and the treatment is
#      interpretable; if it does not, the edit changed something else and the
#      comparison is void.
#    * A LIVENESS CHECK RUNS FIRST (rule 38d). A dead knob yields two identical
#      curves — a clean, confident, entirely wrong negative. We integrate a
#      polynomial the default degree CANNOT do exactly and show the two degrees
#      disagree, before spending an hour on the studies.
#
#  READING THE RESULT:
#    ratios recover to ~8.00  => quadrature crime. A TEST-HARNESS artefact; the
#                                model is fine, and the fix is a higher default
#                                degree for nonlinear runs.
#    ratios unchanged         => NOT quadrature. The advection discretisation
#                                carries a genuine 2nd-order consistency error,
#                                which is a real defect in the operator.
#
#  USAGE: julia --project=. examples/local_mms/run_quadrature_test.jl
# ==============================================================================
using GridapBALFEM, Gridap, Printf, Dates

const OUT = joinpath(@__DIR__, "..", "..", "output", "local", "mms", "quadrature_test")
mkpath(OUT)

# ---- STEP 1: is the knob live? (seconds, and it buys the whole experiment) ----
println("="^78)
println("STEP 1 — quadrature liveness check (rule 38d)")
println("="^78)
let
    model, trian = GridapBALFEM.build_horizontal_model((0.0,1.0,0.0,1.0), (1,1))
    #  ⚠ PICK THE PROBE POLYNOMIAL FROM THE RULE'S ACTUAL EXACTNESS, NOT ITS
    #  NAME. Gridap's `Measure(trian, 8)` is a 5-point Gauss rule, exact to
    #  degree 9 — so an x^9 probe (my first attempt) passes at BOTH degrees and
    #  proves nothing. Measured here: degree 8 -> exact to x^9, degree 10 -> x^11,
    #  degree 12 -> x^13. x^11 therefore separates quad_extra=0 from 4.
    #  This is rule 38d biting the check itself: a liveness test that cannot fail
    #  is worth no more than no test at all.
    f(x) = x[1]^11
    exact = 1/12
    for q in (0, 4)
        deg = 2*max(3,2) + 2 + q
        dΩ  = Measure(trian, deg)
        val = sum(∫(f)*dΩ)
        @printf("  quad_extra=%d  degree=%-3d  ∫x^9 = %.16f   err = %.3e\n",
                q, deg, val, abs(val-exact))
    end
    println("  => the two degrees MUST differ, or the knob is dead and the")
    println("     study below would produce a meaningless negative.\n")
end

# ---- STEP 2: control + treatment, identical but for quad_extra ---------------
#  Parameters copied EXACTLY from the Campaign C model-3 job in
#  run_phaseB_shard.jl (C1/C2 use p_u/levels/nx0/a_eta = 3/5/4/0.8).
common = (; p_u = 3, domain = :d1, mode = :static, levels = 5, nx0 = 4,
            M = 2, p_vert = 1, a_eta = 0.8,
            regime = :nonlinear, flat_bed = true, nl_pressure = :none, a_b = 0.0,
            diag_every = 5, diag_csv = true,
            monitor_factory = () -> SolverMonitor(), verbose = true)

results = Dict{Int,Any}()
for q in (0, 4)
    println("="^78)
    @printf("STEP 2 — nonlinear/flat/:none  Q3/Q2  1-D  quad_extra=%d   (%s)\n",
            q, Dates.format(now(), "HH:MM:SS"))
    println("="^78); flush(stdout)
    t0 = time()
    r = run_conv_study(; common..., quad_extra = q,
                         output_dir = joinpath(OUT, "quad_extra_$q"))
    results[q] = r
    @printf("  elapsed %.1f min\n\n", (time()-t0)/60); flush(stdout)
end

# ---- STEP 3: the comparison --------------------------------------------------
println("="^78)
println("RESULT — e_eta ratios between levels (8.00 = optimal 3rd order)")
println("="^78)
@printf("%-28s %s\n", "Campaign C baseline (q=0)", "7.84  7.58  6.76  5.46")
for q in (0, 4)
    r = results[q]
    rat = [r.e_eta[i]/r.e_eta[i+1] for i in 1:length(r.e_eta)-1]
    @printf("%-28s %s\n", "this run, quad_extra=$q",
            join([@sprintf("%.2f", x) for x in rat], "  "))
end
println()
for q in (0, 4)
    r = results[q]
    @printf("  quad_extra=%d   p_eta pairwise %s  (fit %.3f, optimal 3)\n", q,
            join([@sprintf("%.3f",x) for x in r.pw_eta], ";"), r.fit_eta)
    @printf("               p_u   pairwise %s  (fit %.3f, optimal 4)\n",
            join([@sprintf("%.3f",x) for x in r.pw_u], ";"), r.fit_u)
    @printf("               e_eta %s\n", join([@sprintf("%.4e",x) for x in r.e_eta], " "))
end
println()
d = maximum(abs.(results[0].e_eta .- results[4].e_eta) ./ results[0].e_eta)
if d < 1e-10
    println("⚠ CONTROL AND TREATMENT ARE IDENTICAL -> the knob is DEAD, not a finding.")
elseif results[4].pw_eta[end] > 2.9
    println("✅ QUADRATURE CRIME CONFIRMED: p_eta recovers to ",
            @sprintf("%.3f", results[4].pw_eta[end]),
            " with exact integration.\n   The order loss is a TEST-HARNESS artefact; the operator is fine.")
else
    println("⛔ NOT QUADRATURE: p_eta stays at ",
            @sprintf("%.3f", results[4].pw_eta[end]),
            " with the integrand integrated exactly.\n   The advection discretisation carries a genuine 2nd-order consistency error.")
end
println("\nmax relative change in e_eta between q=0 and q=4: ", @sprintf("%.3e", d))
