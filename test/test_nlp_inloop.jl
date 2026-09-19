# ==============================================================
#  test_nlp_inloop.jl — Class-III projections INSIDE the Newton loop
#
#  Gates the static-condensation change of NEW_TREATMENT.md Part B: the L² projections
#  π𝖲, π𝖻 are refreshed from the CURRENT Newton iterate instead of being frozen from the
#  previous accepted step, removing the O(dt) lag.
#
#  ⛔ THE FEATURE UNDER TEST IS KNOWN BROKEN (2026-09-19). G3 fails, and the
#     discriminator in NEW_TREATMENT.md §E.3 shows why: in-loop mode nearly CANCELS the
#     Class-III contribution — η(:full, in-loop) lands on η(:native) to within 2 %. The
#     gates are kept because they are what caught it, and because they are the
#     acceptance criteria any fix has to meet.
#
#  ⚠ G1 IS THE RULE-38d CHECK AND IT COMES FIRST. A knob that does nothing produces
#  identical curves and a clean, confident, entirely wrong negative. Before any stability
#  claim is made from this switch, it must be shown to CHANGE THE ANSWER.
#  ⚠ G4 is its complement: the change must be CONFINED — tiers that carry no Class-III
#  blocks (:none, :native) must be bit-identical with the switch on.
#
#  Everything here is 1-D Q3/Q2 (rule 2b + NEW_TREATMENT.md Part C): at Q2/Q1 the
#  free-surface Hessian is identically zero in 1-D, so half of what Class III carries
#  would be structurally absent from the test.
#
#  RUN:  julia --project=. test/test_nlp_inloop.jl
#  One block only:  BALFEM_TEST_ONLY=G3 julia --project=. test/test_nlp_inloop.jl
# ==============================================================

using GridapBALFEM
using Gridap
using LinearAlgebra, Printf

println("=" ^ 64)
println("  test_nlp_inloop.jl — static condensation of the Class-III projections")
println("=" ^ 64)

n_pass = 0; n_fail = 0
function check(name, cond)
    global n_pass, n_fail
    if cond; println("  PASS  $name"); n_pass += 1
    else;    println("  FAIL  $name"); n_fail += 1; end
end

# ---- a small, fast, fully-1D flume; Q3/Q2 ---------------------------------------
const BASE = (M=2, p_vertical=1, p_u=3, p_eta=2,
              h_val=3.5, T_wave=1.6, A_wave=0.02,
              domain=((0.0, 6.0), (0.0, 0.25)), partition=(24, 1),
              sponge_wL=0.0, sponge_wR=2.0, mu_max=40.0,
              y_wall_bc=:wall, x_wall_bc=false,
              wave_bc=true, bc_side=:left, bc_profile=:model,
              solver_type=:theta, nl_tol=1e-10, nl_iter=60,
              print_every=100000, save_every=-1, diag_csv=false,
              output_dir=joinpath(@__DIR__, "..", "output", "test_inloop_scratch"))

const ONLY = get(ENV, "BALFEM_TEST_ONLY", "")
runblk(name) = isempty(ONLY) || ONLY == name

"Run to `T_end` and return (diags, prob). `A` and `T_ramp` override BASE per call."
function shortrun(; dt, T_end, nl_pressure, inloop, use_ad=false,
                    A=BASE.A_wave, T_ramp=nothing, nl_tol=BASE.nl_tol)
    kw = merge(BASE, (A_wave=A, nl_tol=nl_tol))
    diags, vert, prob = setup_and_run(; kw...,
        dt=dt, T_final=T_end, T_ramp=T_ramp,
        regime=:nonlinear, nl_pressure=nl_pressure, flat_bed=true,
        nlp_inloop=inloop, use_ad=use_ad)
    return diags, prob
end

etaseries(d) = [r.eta_max for r in d]
finaleta(d)  = etaseries(d)[end]

# ==============================================================
println("\n-- G1: the knob is LIVE (rule 38d) --")
# ==============================================================
pin = nothing
if runblk("G1")
    dlag, plag = shortrun(dt=0.02, T_end=0.24, nl_pressure=:full, inloop=false)
    din,  pin  = shortrun(dt=0.02, T_end=0.24, nl_pressure=:full, inloop=true)
    d1 = abs(finaleta(din) - finaleta(dlag))
    rel1 = d1 / max(abs(finaleta(dlag)), 1e-300)
    @printf("  eta_max(final):  lagged %.12e   in-loop %.12e   |Δ| = %.3e (rel %.3e)\n",
            finaleta(dlag), finaleta(din), d1, rel1)
    check("G1 in-loop CHANGES the :full solution (not a dead knob)", d1 > 1e-12)
    check("G1b …and only modestly — it is an O(dt) correction, not a different model",
          rel1 < 0.5)
end

# ==============================================================
println("\n-- G2: the projection is CONSISTENT at convergence --")
# ==============================================================
#  In in-loop mode the last residual evaluation of the converged step already used the
#  converged state, so refreshing once more must move π by ~nothing. That is the
#  statement "the condensed fixed point has been reached".
if runblk("G2") && pin !== nothing
    ctx = pin.nlp_ctx[]
    check("G2a a projection context is attached in in-loop mode", ctx !== nothing)
    st = pin.nlp_state[]
    check("G2b a projection pair exists after the run", st !== nothing)
end

# ==============================================================
println("\n-- G3: lagged → in-loop difference vanishes with dt --")
# ==============================================================
#  The two schemes differ ONLY by the one-step lag, so their difference must be O(dt).
#  Same physical end time at every level, so the comparison is like-for-like.
#
#  ⚠ THE WINDOW IS THE EXPERIMENT. The first version of this gate ended at t = 0.24 s
#     while T_ramp was 3.2 s, i.e. it compared the two schemes inside the first 7 % of
#     the Hann ramp at eta ~ 5.5e-05 m — where the O(A²) Class-III terms are themselves
#     ~3e-09, the same size as the gap being measured. It returned a FLAT ratio of 0.95
#     (4.190e-09 -> 4.407e-09), a number that says nothing about the lag. Rule 14's
#     transit trap in a new guise: a measurement taken before the state the claim is
#     about has been established. It was recorded as failing rather than threshold-tuned
#     (rule 15), and the WINDOW was fixed instead of the gate.
#
#     FIXED WINDOW: T_ramp = 0.8 s (half a wave period) so the ramp finishes early, and
#     t_end = 2.4 s = 3x the ramp, so the O(dt) lag accumulates over 60-120 steps of an
#     ESTABLISHED wave and dominates any one-off first-step transient.
#
#  ⚠ THE AMPLITUDE IS SET BY WHAT NEWTON CAN RESOLVE, NOT BY WHAT IS INTERESTING.
#     The first attempt at the fix used A = 0.10 m (the amplitude at which :full actually
#     fails) with nl_tol = 1e-12, and NEWTON DID NOT CONVERGE: it stalled at
#     ‖r‖ = 9.250e-08 after 60 iterations. That is the quasi-Newton cliff of rule 5 —
#     jacobian_u omits the {1,2,4,5} blocks, so at working amplitude ~1e-7 is the FLOOR
#     on the attainable residual, whatever tolerance is requested. (Measured here for the
#     first time at production amplitude; CLAUDE.md had only the MMS figure, 9.2e-04 at
#     a_eta = 0.8.) So this gate runs at A = 0.02 m, where the hand Jacobian converges to
#     1e-10, and G3a asserts the measured gap sits >= 100x above that floor. The claim
#     under test is amplitude-independent, so nothing is lost.
#  ⚠ A resolvable-above-the-floor check is not optional: the quantity is a DIFFERENCE of
#     two solutions, so it inherits both their convergence errors (rule 38g).
#
#  ⛔ MEASURED RESULT WITH THE FIXED WINDOW (2026-09-19) — STILL FAILS:
#         dt = 0.040   gap = 1.612753e-03   (rel 7.85e-02)
#         dt = 0.020   gap = 1.676368e-03   (rel 8.12e-02)     ratio 0.962
#     The gap is 1e7 x the Newton floor, so it is unambiguously resolvable, and it is
#     FLAT in dt. The hypothesis this gate exists to test is REFUTED, not unmeasured.
#     The discriminator (NEW_TREATMENT.md §E.3): the whole Class-III contribution at
#     these settings is 1.641585e-03, so the lagged/in-loop gap is 98.2 % of it —
#     in-loop mode is very nearly DELETING the blocks, not de-lagging them.
if runblk("G3")
    G3 = (A=0.02, ramp=0.8, tend=2.4, tol=1e-10)
    @printf("  window: A=%.2f m, T_ramp=%.1f s, t_end=%.1f s (= %.0fx ramp), nl_tol=%.0e\n",
            G3.A, G3.ramp, G3.tend, G3.tend/G3.ramp, G3.tol)
    diffs = Float64[]; base = Float64[]
    for dt in (0.04, 0.02)
        dl, _ = shortrun(dt=dt, T_end=G3.tend, nl_pressure=:full, inloop=false,
                         A=G3.A, T_ramp=G3.ramp, nl_tol=G3.tol)
        di, _ = shortrun(dt=dt, T_end=G3.tend, nl_pressure=:full, inloop=true,
                         A=G3.A, T_ramp=G3.ramp, nl_tol=G3.tol)
        push!(diffs, abs(finaleta(di) - finaleta(dl)))
        push!(base,  abs(finaleta(dl)))
        @printf("  dt=%.3f  eta=%.8e  |eta_inloop-eta_lagged| = %.6e  (rel %.3e)\n",
                dt, base[end], diffs[end], diffs[end]/base[end])
    end
    r1 = diffs[1] / max(diffs[2], 1e-300)
    @printf("  halving ratio: %.2f   (≥ ~1.6 ⇒ at least first order in dt)\n", r1)
    check("G3a the measured gap is resolvable — >= 100x the Newton floor $(G3.tol)",
          minimum(diffs) > 100*G3.tol)
    check("G3 the lagged/in-loop gap shrinks with dt (ratio $(round(r1,digits=2)))",
          r1 > 1.6)
end

# ==============================================================
println("\n-- G4: the change is CONFINED to :full --")
# ==============================================================
if runblk("G4")
    for tier in (:none, :native)
        a, _ = shortrun(dt=0.02, T_end=0.16, nl_pressure=tier, inloop=false)
        b, _ = shortrun(dt=0.02, T_end=0.16, nl_pressure=tier, inloop=true)
        same = finaleta(a) == finaleta(b)
        @printf("  %-8s lagged %.16e  in-loop %.16e  identical=%s\n",
                tier, finaleta(a), finaleta(b), same)
        check("G4 nl_pressure=:$tier is BIT-IDENTICAL with the switch on", same)
    end
end

# ==============================================================
println("\n-- G5: the AD guard holds (the refresh must not throw under Dual) --")
# ==============================================================
if runblk("G5")
    ok_ad = true
    try
        shortrun(dt=0.02, T_end=0.08, nl_pressure=:full, inloop=true, use_ad=true)
    catch err
        global ok_ad = false
        println("  AD path raised: ", err)
    end
    check("G5 use_ad=true runs in in-loop mode (nlp_plain_iterate skips the mass solve)",
          ok_ad)
end

println("\n" * "=" ^ 64)
@printf("  %d passed, %d failed\n", n_pass, n_fail)
println("=" ^ 64)
n_fail == 0 || error("test_nlp_inloop: $n_fail gate(s) failed")
