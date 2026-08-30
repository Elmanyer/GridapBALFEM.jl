# ==============================================================
#  test_mms_distributed_parity.jl — the SEQUENTIAL↔DISTRIBUTED MMS gate
#                                       (4 ranks, 2×2)
#
#  THE DEFECT THIS EXISTS TO CATCH (found 2026-08-19, fixed 2026-08-21).
#  `run_mms_case_distributed` hard-coded Model 1 — `regime=:linear`,
#  `nl_pressure=:none`, `flat_bed=true` as LITERALS plus the Model-1 closed form
#  `mms_forcing_stage1` — and accepted no model switches at all, while
#  `run_conv_study` still built its printed tag from the *requested* switches.
#  A distributed campaign over eight models therefore returned EIGHT IDENTICAL
#  MODEL-1 STUDIES UNDER EIGHT DIFFERENT LABELS, and every one of them PASSED,
#  because Model 1 converges optimally. It failed silently, confidently, and
#  inside the one tool whose entire purpose is to catch results that are wrong
#  self-consistently.
#
#  WHAT MAKES THIS GATE REAL, AND NOT ANOTHER GREEN LIGHT ON NOTHING.
#  A bare "distributed agrees with sequential" check would pass WITH THE DEFECT
#  PRESENT for any configuration whose answer happens to coincide with Model 1's
#  — and, run on Model 1 itself, it coincides by definition. So this file gates
#  in two directions and the first one is not optional:
#
#    G0  SEPARATION (the negative control) — Model 1 and the model under test
#        must produce MATERIALLY DIFFERENT errors sequentially. Without this the
#        parity gates below are vacuous: they would confirm agreement between two
#        computations of the same thing. `run_conv_study`'s tag was built from the
#        requested switches and looked perfect throughout; a label is not evidence.
#    G1  PARITY, flat_bed — Model 2 (:linear / VARIABLE bed / :none). Cheapest
#        non-Model-1 configuration (linear ⇒ one Newton step), and it moves the
#        `flat_bed` switch AND forces the general `mms_forcing` path instead of
#        the Model-1 closed form.
#    G2  PARITY, regime × nl_pressure — Model 5 (:nonlinear / flat / :native).
#        The other two switches, which G1 leaves at their Model-1 values.
#
#  Together G1 and G2 move all three switches off their hard-coded values, so the
#  defect cannot recur in any one of them without a gate turning red.
#
#  ⚠ THE REFERENCE IS COMPUTED IN-PROCESS, NEVER PINNED AS A CONSTANT. Three
#  stale reference constants in this suite were found the first time their tests
#  were re-run in a session (OPEN_ITEMS.md §3); a reference that recomputes cannot
#  go stale. Every rank computes it redundantly — cheap at this size, and it keeps
#  the comparison exact rather than to printed precision.
#
#  The vertical basis is a PARAMETER here as everywhere else (MMS_M / MMS_PVERT):
#  parity between the two execution paths is a property of the code, not of the
#  basis, so it must hold for any (M, p).
#
#  RUN:
#    ~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. \
#        test/test_mms_distributed_parity.jl
# ==============================================================

using GridapBALFEM
using Printf

const M_VERT = parse(Int, get(ENV, "MMS_M",     "2"))
const P_VERT = parse(Int, get(ENV, "MMS_PVERT", "1"))
const NX     = parse(Int, get(ENV, "MMS_PAR_NX", "12"))
const NY     = parse(Int, get(ENV, "MMS_PAR_NY", "8"))

const LX, LY = 1.7, 1.1
const DEPTH  = 2.5          # d ≠ 1: multiplication by h is then not the identity
const A_BED  = 0.2
const P_U, P_ETA = 3, 2     # Q3/Q2, the pairing used everywhere else

#  Tolerances, and a MEASURED WARNING ABOUT THIS TEST'S SENSITIVITY.
#
#  The two paths differ ONLY in the linear solve — direct LU against GMRES(Jacobi)
#  at rtol 1e-13 — so their L² errors agree to near round-off. But the SEPARATION
#  between models is far smaller than intuition suggests, and it was measured
#  before this gate was trusted (2026-08-21, nx=12 Q3/Q2, d=2.5, a_b=0.2):
#
#      Model 1 (lin/flat/:none)      e_eta=1.9844455e-04   e_u=9.859119e-06
#      Model 2 (lin/varbed/:none)    e_eta=1.9844463e-04   e_u=9.857747e-06
#      Model 5 (nl/flat/:native)     e_eta=1.9844477e-04   e_u=9.862121e-06
#                                    Δ ≈ 4e-07 rel        Δ ≈ 1.4e-04 / 3.1e-04 rel
#
#  ⚠ THE MMS L² ERROR IS A WEAK DISCRIMINATOR BETWEEN MODEL TIERS, and the reason
#  is structural: the MMS FORCES the same manufactured field for every model, so
#  the exact solution is identical and the models differ only through the
#  discretisation error of their different operators — a small perturbation on a
#  resolved mesh. **e_eta separates by only ~4e-7 and must NOT be used for G0.**
#  e_u separates by ~1e-4 and is the usable channel.
#
#  Hence: PARITY_RTOL sits well below the measured separation, and G0 does not
#  merely assert "the models differ" — it asserts the SEPARATION EXCEEDS THE
#  PARITY TOLERANCE BY A STATED MARGIN, i.e. that this test can actually resolve
#  the defect it exists to catch. A gate that cannot resolve its target is decor.
const PARITY_RTOL   = 1e-7    # LU vs GMRES(1e-13) on the SAME problem
const MIN_MARGIN    = 30.0    # required separation / PARITY_RTOL

#  If a future change makes the separation shrink below the margin, the honest fix
#  is a MORE DISCRIMINATING CASE, not a looser gate. The strongest one available is
#  `nl_pressure=:full`, whose e_u sits on a frozen-projection floor ~600x above the
#  :none value — but its LU-vs-CG projection difference then needs a parity
#  tolerance nearer 5e-3 (see test_nlpressure_distributed.jl).

is_rank0 = get(ENV, "OMPI_COMM_WORLD_RANK", get(ENV, "PMI_RANK", "0")) == "0"
n_pass = 0; n_fail = 0
function check(name, cond, extra = "")
    global n_pass, n_fail
    cond ? (n_pass += 1) : (n_fail += 1)
    is_rank0 && println(cond ? "  PASS  $name $extra" : "  FAIL  $name $extra")
    is_rank0 && flush(stdout)
end

#  ω = 0 ⇒ ∂ₜu* ≡ 0 and the temporal discretisation error is IDENTICALLY zero, so
#  what the two paths are compared on is the spatial operator alone. A parity
#  difference then cannot be blamed on the integrator.
vert  = assemble_vertical_tensors(M_VERT, P_VERT, resolve_cbdy(M_VERT))
field = MMSField(vert.N_dof; Lx=LX, Ly=LY, omega=0.0)
bed   = bathymetry_field(; d0=DEPTH, a_b=A_BED, kbx=1.3, kby=0.0)

#  ONE parameter set for both paths. Anything that differs between them other than
#  the linear solver would make a parity difference uninterpretable — two tests
#  sharing a bathymetry and a physics tier can still be different problems.
common = (; nx=NX, ny=NY, dt=1e-5, T_final=1e-4, Lx=LX, Ly=LY, d=DEPTH,
            vert_override=vert, p_horizontal=P_U, p_eta=P_ETA, field=field,
            nl_tol=1e-12, nl_iter=400, verbose=false)

models = (
    m1     = (; regime=:linear,    flat_bed=true,  nl_pressure=:none,   hfun=nothing),
    model2 = (; regime=:linear,    flat_bed=false, nl_pressure=:none,   hfun=bed),
    model5 = (; regime=:nonlinear, flat_bed=true,  nl_pressure=:native, hfun=nothing),
)

if is_rank0
    println("=" ^ 76)
    println("  test_mms_distributed_parity.jl — sequential vs distributed MMS (4 ranks)")
    println("=" ^ 76)
    @printf("  vertical basis P%dLFE-%d (Nσ=%d) | %dx%d | Q%d/Q%d | d=%.1f | static (ω=0)\n",
            P_VERT, M_VERT, vert.N_dof, NX, NY, P_U, P_ETA, DEPTH)
    flush(stdout)
end

# ---------------------------------------------------------------------------
#  Sequential references — every model, computed here and now.
# ---------------------------------------------------------------------------
is_rank0 && println("\n--- sequential references ---")
seq = map(m -> run_mms_case(; common..., m...), models)
if is_rank0
    for k in keys(seq)
        @printf("  %-7s e_eta=%.10e  e_u=%.10e\n", k, seq[k].e_eta, seq[k].e_u)
    end
    flush(stdout)
end

# ---------------------------------------------------------------------------
#  G0 — SEPARATION. Read this one first: if it fails, G1/G2 mean nothing.
# ---------------------------------------------------------------------------
is_rank0 && println("\n--- G0: separation — can this test RESOLVE a silent fallback to Model 1? ---")
for k in (:model2, :model5)
    de = abs(seq[k].e_eta - seq.m1.e_eta) / seq.m1.e_eta
    du = abs(seq[k].e_u   - seq.m1.e_u)   / seq.m1.e_u
    sep = max(de, du)                      # e_u is the usable channel; see the note
    check("G0 $k separation exceeds $(MIN_MARGIN)x the parity tolerance",
          sep > MIN_MARGIN * PARITY_RTOL,
          @sprintf("(sep=%.3e = %.0fx parity; e_eta %.2e / e_u %.2e)",
                   sep, sep/PARITY_RTOL, de, du))
end
if n_fail > 0 && is_rank0
    println("  ⇒ STOP. The parity gates below cannot resolve the hard-coded-Model-1")
    println("    defect on this case: a distributed run silently computing Model 1")
    println("    would land INSIDE the parity tolerance and be reported as a PASS.")
    println("    Fix by making the case MORE DISCRIMINATING (coarser mesh, larger")
    println("    a_b, or nl_pressure=:full) — never by loosening PARITY_RTOL.")
    flush(stdout)
end

# ---------------------------------------------------------------------------
#  G1 / G2 — PARITY. Same model, both execution paths.
# ---------------------------------------------------------------------------
for (gate, k) in (("G1", :model2), ("G2", :model5))
    m = models[k]
    is_rank0 && println("\n--- $gate: $k  $(m.regime) / $(m.flat_bed ? "flat" : "varbed") / $(m.nl_pressure) ---")
    is_rank0 && flush(stdout)
    dst = run_mms_case_distributed(; common..., m..., cpu_grid=(2,2),
                                     ls_rtol=1e-13, ls_maxiter=5000)
    re = abs(dst.e_eta - seq[k].e_eta) / seq[k].e_eta
    ru = abs(dst.e_u   - seq[k].e_u)   / seq[k].e_u
    check("$gate e_eta parity (rel < $PARITY_RTOL)", re < PARITY_RTOL,
          @sprintf("(seq %.10e  dist %.10e  rel %.2e)", seq[k].e_eta, dst.e_eta, re))
    check("$gate e_u   parity (rel < $PARITY_RTOL)", ru < PARITY_RTOL,
          @sprintf("(seq %.10e  dist %.10e  rel %.2e)", seq[k].e_u, dst.e_u, ru))
    #  A distributed run that silently computed Model 1 would land on the Model-1
    #  reference instead. Name that explicitly, so the failure mode is readable
    #  from the log rather than only from a number that is merely "off".
    if re ≥ PARITY_RTOL || ru ≥ PARITY_RTOL
        d1e = abs(dst.e_eta - seq.m1.e_eta) / seq.m1.e_eta
        d1u = abs(dst.e_u   - seq.m1.e_u)   / seq.m1.e_u
        if is_rank0 && d1e < PARITY_RTOL && d1u < PARITY_RTOL
            println("  ⇒ THE DISTRIBUTED RUN MATCHED MODEL 1, not the requested model.")
            println("    That is the 2026-08-19 defect exactly: check that")
            println("    run_mms_case_distributed forwards regime/flat_bed/nl_pressure/hfun")
            println("    and calls mms_forcing, not mms_forcing_stage1.")
        end
    end
end

if is_rank0
    println()
    println("=" ^ 76)
    @printf("  Results: %d PASS,  %d FAIL\n", n_pass, n_fail)
    println("=" ^ 76)
    flush(stdout)
end
n_fail > 0 && error("test_mms_distributed_parity: $n_fail failed!")
is_rank0 && println("  Distributed MMS reproduces the sequential result for the REQUESTED model.")
