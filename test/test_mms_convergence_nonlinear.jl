# ==============================================================
#  test_mms_convergence_nonlinear.jl — ORDER OF ACCURACY for the nonlinear
#  analytic-MMS models (the verification proper)
#
#  Specification: ValidationTests.tex §subsec: mms model3 / §subsec: mms model4
#                 and §subsubsec: mms measure.
#  Plan:          markdown_files/VERIFIED_SCOPE.md (the nonlinear plan itself is gone) §2.3.
#
#  WHAT THIS CERTIFIES, AND WHY IT IS DIFFERENT FROM test_selfconsistency.jl.
#  The forcing is derived from the GOVERNING EQUATIONS (src/mms.jl), never from
#  the residual, so a wrong term in problem.jl does NOT cancel — it shows up as a
#  reduced order of accuracy. Reaching the theoretical order is therefore code
#  verification in the sense of Roache, not self-consistency.
#
#  PAIRING. Q3/Q2, i.e. velocity one order above the surface. Equal order costs
#  one order in u (measured); on Q_p/Q_{p−1} both fields reach their optimum:
#      p_η → p_e+1 = 3      p_u → p_u+1 = 4
#
#  TOLERANCE. nl_tol = 1e-9, NOT 1e-12. The nonlinear Jacobians are deliberately
#  quasi-Newton (the pressure blocks' η-dependence is frozen — see the coverage
#  note in problem.jl), so Newton converges linearly and stalls around 1e-10.
#  1e-9 sits two to three orders below the smallest discretisation error measured
#  here, so it cannot affect the rate.
#
#  NEVER loosen nl_tol to make a study pass: the algebraic error would then sit
#  INSIDE the discretisation error being measured and the rate would be
#  meaningless.
#
#  ⚠ AND DO NOT REACH FOR ITERATION BUDGET EITHER, WITHOUT FIRST ASKING WHICH
#  KIND OF NON-CONVERGENCE YOU HAVE. A Newton solve can fail two ways that look
#  identical in the log: converging SLOWLY (a higher-order term missing from the
#  Jacobian — more iterations fix it) or converging to the WRONG FIXED POINT (an
#  O(1) term missing — no budget ever fixes it). This file briefly carried
#  nl_iter=400 for Model 4 on the first reading; the truth was the second, and
#  the real fix was to assemble the 𝓐/𝓚 package in ∂R/∂u̇ (2026-08-17).
#  test_jacobians_ad.jl distinguishes the two directly, by measuring how the
#  hand↔AD gap scales with state amplitude: vanishing ⇒ slow, flat ⇒ wrong.
#
#  RUNTIME. ~25–30 min per nl_pressure=false study (3 levels). The nl_pressure=true
#  studies are SUBSTANTIALLY slower — the 𝓝 forcing costs several × the 𝓝-free one
#  (Ψ carries an Nσ²×8 component sum that the outer gradient then differentiates)
#  and the solver assembles the broken Class-III skeleton. Set MMS_NL_LEVELS to shorten.
#
#  THE VERTICAL BASIS IS A PARAMETER: MMS_M (elements) and MMS_PVERT (order),
#  defaulting to the P1LFE-2 every study in this repository has ever run. The
#  rates asserted below are properties of the HORIZONTAL discretisation, so they
#  must hold for ANY vertical basis — overriding these is how the basis-agnosticism
#  the model is named for gets tested (markdown_files/COMPLETED_VBASIS_STUDY.md §1).
#  ⚠ The forcing cost scales as Nσ²; P1LFE-4 or P2LFE-3 is several times the
#  default. Prune with MMS_NL_LEVELS before sweeping.
#
#  ENV
#    MMS_NL_LEVELS  refinement levels                       3
#    MMS_NL_NX0     coarsest nx                             8
#    MMS_M          vertical elements M                     2
#    MMS_PVERT      vertical FE order p  (Nσ = M·p+1)       1
#    MMS_NL_P       1 ⇒ also run the two nl_pressure=true models  0  (opt-in, see below)
#
#  RUN:  julia --project=. test/test_mms_convergence_nonlinear.jl
#        MMS_M=3 julia --project=. test/test_mms_convergence_nonlinear.jl
# ==============================================================

using GridapBALFEM
using Printf

println("=" ^ 76)
println("  test_mms_convergence_nonlinear.jl — order of accuracy, nonlinear models")
println("=" ^ 76)

const LEVELS = parse(Int, get(ENV, "MMS_NL_LEVELS", "3"))
const NX0    = parse(Int, get(ENV, "MMS_NL_NX0",    "8"))
const TOLP   = 0.3          # |p_obs − p_opt| gate, as in §subsubsec: mms measure
#  VERTICAL basis — see the header. c_bdy is left to resolve_cbdy, which is what
#  makes M ≠ 2 legal at all: run_conv_study used to hard-wire the M=2 node set and
#  threw the length(c_bdy)==M+1 assertion for anything else (fixed 2026-08-21).
const M_VERT = parse(Int, get(ENV, "MMS_M",     "2"))
const P_VERT = parse(Int, get(ENV, "MMS_PVERT", "1"))
const RUN_NLP  = get(ENV, "MMS_NL_P", "0") != "0"
@printf("  vertical basis: P%dLFE-%d  (Nσ = %d)   levels=%d  nx0=%d\n",
        P_VERT, M_VERT, M_VERT*P_VERT + 1, LEVELS, NX0)

n_pass = 0; n_fail = 0
function check(name, cond, extra = "")
    global n_pass, n_fail
    cond ? (println("  PASS  $name $extra"); n_pass += 1) :
           (println("  FAIL  $name $extra"); n_fail += 1)
end

#  Studies, simplest first, so a broken rate is attributable to the block just added.
#  Six-model numbering (V2_SOLVER_PLAN.md §1): models 3–4 are nl_pressure=false, models 5–6
#  nl_pressure=true (all eight 𝓝 components; Class III by the broken formulation, which is
#  exact, so both fields are rate-gated). Read the pairwise sequence, not only the fit (rule 33).
studies = [
    (name = "Model 3  nonlinear / flat bed      / nlp0",
     regime = :nonlinear, flat_bed = true,  nl_pressure = false, nl_iter = 50, a_eta = 0.8,
     rate_gate_u = true),
    #  Model 4 runs at the DEFAULT budget. It briefly carried nl_iter=400, added
    #  when the stall at ‖r‖=4.8e-8 was read as "quasi-Newton convergence is just
    #  slow, give it more iterations". That diagnosis was WRONG and the extra
    #  budget would never have helped: ∂R/∂u̇ was missing the 𝓐/𝓚 slope package,
    #  whose prefactor H·∇h does NOT scale with the solution, so Newton was
    #  converging to a fixed point of the wrong map — an O(1) error, not a slow
    #  one. With that block assembled (2026-08-17) Model 4 converges in the
    #  ordinary number of iterations. Lesson: distinguish "converging slowly"
    #  from "converging to the wrong thing" BEFORE spending iterations on it;
    #  test_jacobians_ad.jl tells them apart by amplitude scaling.
    (name = "Model 4  nonlinear / variable bed  / nlp0",
     regime = :nonlinear, flat_bed = false, nl_pressure = false, nl_iter = 50, a_eta = 0.8,
     rate_gate_u = true),
]
#  ⚠ Models 5–6 are OPT-IN (MMS_NL_P=1) and NOT part of the verified scope: the full nonlinear
#  model is still under development and not ready for MMS evaluation (decision 2026-10-05). The
#  one exploratory run measured p_η 2.997/2.996 but p_u pairwise 3.94 → 1.97 (Model 5, nx 8–32).
#  Models 5–6 at a_eta = 0.4: the {3,6,7,8} and ∇h-IBP 𝓝 blocks are quasi-Newton, and at the
#  default a_eta = 0.8 (H_min = 0.2·d) Newton stalls above nl_tol. Lower the amplitude, never
#  loosen nl_tol.
RUN_NLP && append!(studies, [
    (name = "Model 5  nonlinear / flat bed      / nlp1", regime=:nonlinear,
     flat_bed=true,  nl_pressure=true, nl_iter=50, a_eta = 0.4, rate_gate_u = true),
    (name = "Model 6  nonlinear / variable bed  / nlp1", regime=:nonlinear,
     flat_bed=false, nl_pressure=true, nl_iter=50, a_eta = 0.4, rate_gate_u = true),
])
RUN_NLP || println("\n  [skip] models 5–6 (nl_pressure=true) — opt-in, MMS_NL_P=1; not yet in the verified scope.")

for s in studies
    println("\n" * "-"^76)
    println("  $(s.name)")
    println("-"^76); flush(stdout)
    #  A study that cannot COMPLETE must be reported as a failed gate, not allowed
    #  to abort the run. `run_conv_study` throws when Newton exhausts its budget,
    #  and Model 4 is exactly the study expected to do that — letting it propagate
    #  would destroy the report of every study that DID work, including Model 3's
    #  verification result. Report it, keep going, fail at the end.
    r = try
        run_conv_study(; p_u = 3, domain = :d1, mode = :static,
                         levels = LEVELS, nx0 = NX0, ny_1d = 3,
                         Lx = 1.7, Ly = 1.1, d = 2.5,
                         M = M_VERT, p_vert = P_VERT,   # c_bdy ⇒ resolve_cbdy(M)
                         dt = 1e-5, nsteps = 100, nl_tol = 1e-9,
                         nl_iter = s.nl_iter,
                         regime = s.regime, nl_pressure = s.nl_pressure, a_eta = s.a_eta,
                         flat_bed = s.flat_bed, a_b = s.flat_bed ? 0.0 : 0.2,
                         kbx = 1.3, kby = 0.0, verbose = true)
    catch e
        msg = first(split(sprint(showerror, e), '\n'))
        check("$(s.name): study completed", false, "\n        $msg")
        if occursin("did not converge", msg)
            println("        ⇒ the SOLVE is under-converged, which is NOT evidence about the")
            println("          forcing or the operator: a stalled Newton reads exactly like a")
            println("          wrong residual. The nonlinear Jacobians are quasi-Newton by")
            println("          design, so convergence is linear. Raise nl_iter (currently")
            println("          $(s.nl_iter)) and re-run. If the residual PLATEAUS rather than")
            println("          decreasing slowly, budget will not help and the real fix is to")
            println("          complete the nonlinear Jacobians — see MMS_NONLINEAR_PLAN.md §5b.")
            println("        ⇒ do NOT loosen nl_tol to make this pass: the algebraic error would")
            println("          then sit inside the discretisation error being measured.")
        end
        flush(stdout)
        continue
    end
    check("$(s.name): p_eta → $(Int(r.opt_eta))", abs(r.fit_eta - r.opt_eta) < TOLP,
          @sprintf("(got %.3f)", r.fit_eta))
    if s.rate_gate_u
        check("$(s.name): p_u   → $(Int(r.opt_u))",   abs(r.fit_u   - r.opt_u)   < TOLP,
              @sprintf("(got %.3f)", r.fit_u))
    end
    flush(stdout)
end

println()
println("=" ^ 76)
@printf("  Results: %d PASS,  %d FAIL\n", n_pass, n_fail)
println("=" ^ 76)
n_fail > 0 ? error("test_mms_convergence_nonlinear: $n_fail failed!") :
             println("  Nonlinear MMS models reach their theoretical order of accuracy.")
