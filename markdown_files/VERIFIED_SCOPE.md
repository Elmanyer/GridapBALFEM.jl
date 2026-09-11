# VERIFIED_SCOPE.md — what is actually proven about the residual

> **Scope — what belongs here.** The *record of what is proven*: the analytic MMS campaign, the
> Jacobian-vs-AD oracle, the verified model scope and its boundaries. Answers **"is the solver
> correct, and how far does that claim reach?"**
>
> **What does NOT belong here:** the inventory of test files → [`TEST_SUITE.md`](TEST_SUITE.md);
> known defects → [`OPEN_ISSUES.md`](OPEN_ISSUES.md); work not yet run →
> [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md).


**What this file is.** The verification record: the analytic MMS campaign, the Jacobian-vs-AD
oracle, the verified scope, and the scope boundaries. It is the answer to "is the solver correct?",
stated precisely enough to be quotable.

Related: [`MODEL.md`](MODEL.md) (the operator being verified), [`TEST_SUITE.md`](TEST_SUITE.md)
(the full gate inventory).

---

## 0. "Validated" vs "verified" — read this before claiming correctness

Most of the suite is **self-consistency**: each test compares the solver against an analytical or
physical expectation computed *using the solver's own residual*. Writing `R = R_true + E`, the error
`E` appears on both sides and cancels identically — **such a test passes for any residual, right or
wrong.**

That is still strong evidence: a wrong residual would have to be wrong in a way that simultaneously
preserves the dispersion relation, the shallow-water limit, the Airy vertical profile, the rest
state, mass conservation and the `O(A²)` ordering of every physics tier. But it is not verification.

**The analytic MMS is the only instrument here without that property.** Its forcing is derived from
the governing equations in closed form and **never touches `problem.jl`** (enforced by a grep gate
in `test_mms_forcing.jl`), so the error does not cancel. Reaching the **theoretical order of
accuracy** certifies that the discretised operator *is* the intended operator — verification in
Roache's sense.

**The supporting logic, because it is the crux:** the operator code is byte-identical across runs
that differ only in FE spaces or mesh, and a wrong coefficient cannot be repaired by changing
function spaces. So recovering the theoretical rate under a stable pairing proves the operator was
right; a collapsed rate on an *isolated* study proves it was not.

---

## 1. Verified scope — SIX of the eight models

Measured `Q3/Q2`, 1-D static unless noted.

| model | `regime` / `flat_bed` / `nl_pressure` | `p_η` (opt 3) | `p_u` (opt 4) | |
|---|---|---|---|---|
| 1 | `:linear` / flat / `:none` | `p_e+1` exactly | `p_u+1` exactly | ✅ |
| 2 | `:linear` / **variable** / `:none` | **3.000** | **4.000** | ✅ |
| 3 | `:nonlinear` / flat / `:none` | **2.996** | **3.995** | ✅ |
| 4 | `:nonlinear` / **variable** / `:none` | **2.996** | **3.997** | ✅ |
| 5 | `:nonlinear` / flat / **`:native`** | **2.996** | **3.997** | ✅ |
| 6 | `:nonlinear` / **variable** / **`:native`** | **2.996** | **3.998** | ✅ |
| 7–8 | `:nonlinear` / any / **`:full`** | 2.99 → 2.59 | **−0.00** | ⛔ not MMS-verifiable **by construction** — §4 |

Model 2 additionally confirmed **transient** (`2.999`/`3.998`) and in **2-D** (`3.000`/`3.963`).

**So the verified scope is the complete `:none` model over arbitrary bathymetry, PLUS the whole
`:native` nonlinear-pressure tier (`𝓝` components `{3,6,7,8}`) on both flat and variable beds** —
`H`-weighting, advection, the full three-component leading pressure, the `O(ε²)` surface-slope
package and the bed-slope terms.

> **Say "the `:none` and `:native` models are verified".** Never bare *"the residual is verified"*
> (that would wrongly include `:full`), and never *"the `𝓝` tiers are verified"* (same error).

Model 4's final levels, `nl_iter=50`:

| `nx` | `e_eta` | rate | `e_u` | rate |
|---|---|---|---|---|
| 8 | 2.688582e-04 | — | 1.285103e-05 | — |
| 16 | 3.368105e-05 | 2.997 | 8.058666e-07 | 3.995 |
| 32 | 4.221570e-06 | 2.996 | 5.040873e-08 | 3.999 |

---

## 2. How the MMS is built (`src/mms.jl`)

**`strong_residual_model` is a SINGLE parent evaluator for all models**, with
`regime`/`flat_bed`/`nl_pressure` applied at **one control point each**, mirroring `resolve_physics`
so a forcing/model mismatch is unrepresentable. `strong_residual_linear` is a thin wrapper over it
(agrees with the retained legacy implementation to 2.1e-17). `mms_forcing` is the single entry
point, memoised one point deep (Gridap evaluates `Seta`/`Sx`/`Sy` at the same point consecutively,
so 3 evaluations collapse to 1) and `H>0` guarded.

Model 1 keeps a hand-written closed form because it is ~7× faster than the AD parent **and** because
their agreement (1.8e-15) is itself a check on both.

**Measured forcing cost per point** — this is why the `𝓝` rate studies take ~1 h each:

| tier | cost/point |
|---|---|
| closed form (Model 1) | ~3.5 µs |
| AD `:none` | ~25 µs |
| `:native` | ~100–150 µs |
| `:full` | ~230–300 µs |

The blow-up is **not** AD nesting per se: `Ψ` carries an `Nσ²×8` component sum that the outer
gradient then differentiates.

**Run the forcing-level gates first** — they cost seconds and bracket each new model between two
already-verified ones, so a wrong term must break one of them:

| gate | result |
|---|---|
| Model 4 at constant `h` ≡ Model 3 | **exactly 0.0** |
| Model 2 at constant `h` ≡ Model 1 | **exactly 0.0** |
| nonlinear−linear gap scales `ε²` | ratio **3.913** (expect 4) |
| eigenmode `𝓛(u*)` | 3.6e-15 |
| closed form ≡ ForwardDiff | 1.8e-15 |
| `flat_bed` control point holds under `𝓝` | rel **0.000e+00**, both tiers |
| `𝓝` excess is `O(A²)` | 1.927 / 1.998 / 1.968 / 1.967 |
| grep gate: `src/mms.jl` never references residual code | PASS |

The last `𝓝` result independently corrects a former blanket `O(A³)` claim for the whole package —
the `𝓟` block is `O(A²)` and dominates.

---

## 3. The FE pairing costs one convergence order

**Equal order `Q_p/Q_p` converges at `p`, not `p+1`.** Mixed order `Q_p/Q_{p-1}` fixes the surface
universally and the velocity at `Q3/Q2`. Established by a 12-study campaign (3 pairings × 1-D/2-D ×
static/transient, 4 mesh levels):

| pairing | domain | `p_η` (opt) | `p_u` (opt) | verdict |
|---|---|---|---|---|
| Q2/Q1 | 1-D | **2.000** (2) | 2.237 (3) | η ✓ · u short |
| Q2/Q1 | 2-D | **1.980** (2) | 2.602 (3) | η ✓ · u short |
| **Q3/Q2** | **1-D** | **3.000** (3) | **3.991** (4) | **both optimal** |
| **Q3/Q2** | **2-D** | **2.997** (3) | **3.930** (4) | **both optimal** |
| Q4/Q3 | 1-D | **4.000** (4) | — (5) | η ✓ · u at the **round-off floor** |
| Q4/Q3 | 2-D | **3.994** (4) | 4.735 (5) | η ✓ · u short |

**Established:** `η` reaches its optimal `p_e+1` in **12/12** studies to 3–4 significant figures;
`Q3/Q2` is optimal in **both** fields across four independent studies. Static (`ω=0`) and transient
agree to four significant figures in all 12, so the rates are clean *spatial* rates.

**NOT established — do not claim `p+1` for velocity as a general rule.** `u` falls short at `Q2/Q1`
(→2.4 vs 3) and `Q4/Q3` (→4.65 vs 5), with pairwise rates still *falling* at the finest level.
Genuinely suboptimal vs merely pre-asymptotic is undetermined.

**Cause.** `η` enters momentum undifferentiated via `∇·v` after integration by parts, so it plays
the pressure role and equal-order continuous spaces are inf-sup deficient (the Stokes analogue);
⚠ **as of 2026-09-06 Taylor-Hood (`p_u = p_eta + 1`) is REQUIRED and enforced by
`check_taylor_hood`** — equal order is what produced the year-long "nonlinear instability"
(`CLAUDE.md` rules 2b, 12b). Every result in this document was measured on Taylor-Hood;
velocity one order above the surface is the Taylor–Hood pairing.

> **Prefer `Q3/Q2` — but do NOT switch production on the rate alone.** At `nx=24`, `Q3/Q3` gave
> `e_eta = 5.97e-7` against `Q3/Q2`'s `2.40e-5` — **40× more accurate at that mesh** despite the
> worse rate, because `η` sits in a richer space. An error-vs-DOF study at production resolution is
> **still not run** and must precede any default change. Equal order remains the default;
> `p_eta` / `BALFEM_P_ETA` opts in.

Time integration and `R_P` were each eliminated as explanations by experiment first: a steady `ω=0`
field reproduces the same rates to 3 decimals, and `B=0` makes `η` *worse*, so `R_P` is stabilising.

---

## 4. ⛔ `:full` can never pass an MMS rate test — and that is a design fact

> 🔴 **THE EXPLANATION BELOW IS WRONG IN ONE LOAD-BEARING DETAIL — found 2026-08-21, not yet
> re-derived.** The numbers stand; their attribution does not.
>
> This section says the floor comes from projections **lagged one step**. In the MMS driver they are
> **not lagged — they are ABSENT**. `run_time_loop` defaults `nlp = nothing` (`src/timeloop.jl:130`)
> and `run_mms_case` never passes one (`src/mms_driver.jl:99–105`; only `setup_and_run` builds a
> context, `src/utilities.jl:863`), so `src/problem.jl:396`'s `st !== nothing` gate means the frozen
> `{1,2,4,5}` contribution is never assembled at all. The `nl_pressure_full` branch then adds
> **nothing whatsoever on a flat bed**, and only the `𝓐` ∇h IBP half on a sloping one.
>
> The studies remain valid and the floors are real — the *forcing* is still selected by
> `nl_pressure` and still computes `{1,2,4,5}` exactly, so forcing and solver do encode different
> operators. But the difference measured is **omission**, not **lag**. Two claims below therefore do
> not follow from this experiment:
> * that the lag is "negligible … now **measured** rather than argued" — the lag was never exercised;
> * that this floor characterises a **production** `:full` run, where `setup_and_run` *does* build
>   the context and the projections genuinely are frozen-and-lagged.
>
> **Decisive next step:** give `run_mms_case` an `nlp` kwarg (or build the context internally when
> `nl_pressure=:full`), re-run the two `:full` studies, and compare. If the floor drops sharply, it
> was omission; if it barely moves, the original wording was accidentally right. Until then, quote
> the floor as *"the magnitude of the `{1,2,4,5}` terms the MMS path omits"*.

The solver's `:full` does **not** implement the exact `𝓝` operator: components `{1,2,4,5}` carry
irreducible `∂²η`, so their `∇H` half and `𝓟` part are evaluated from **frozen L² projections**
lagged one step *(in a production run; see the box above for what the MMS path actually does)*. The
analytic MMS forcing computes those same components **exactly**. Forcing and solver therefore encode
*different operators*, and the study measures the surrogate.

| model / tier | `e_u` @ nx=8 | @ nx=16 | @ nx=32 | `p_u` |
|---|---|---|---|---|
| M3 `:native` | 1.320e-05 | 8.295e-07 | 5.194e-08 | **3.99** |
| M4 `:native` | 1.291e-05 | 8.097e-07 | 5.067e-08 | **4.00** |
| M3 `:full` | 5.983e-03 | 5.988e-03 | 5.988e-03 | **−0.00** |
| M4 `:full` | 5.444e-03 | 5.448e-03 | 5.448e-03 | **−0.00** |

The two `:full` floors agree within **10 %** despite one bed being flat and the other sloping —
itself the evidence that the floor is a property of the **approximation**, not of the bathymetry or
of any `∇h` term. (`e_η` degrades 2.994 → 2.673 as the floor is approached, for the same reason.)

**The useful result is the number, not the rate.** `e_u` settles at **5.988e-03**, moving 0.081 %
across a 16× mesh refinement — a **mesh-independent accuracy floor**. At `nx=32` that is **~10⁵×
larger than `:native`'s discretisation error on the same mesh**: a `:full` run's velocity error is
entirely the frozen-projection approximation, and refining the mesh will not reduce it. The standing
claim that the lag is "`O(dt)` on an already `O(A²)` term, so negligible" is now **measured rather
than argued** — and it is negligible only relative to `O(A²)`, not relative to the discretisation
error at production resolution.

> **Consequence for production: `nl_pressure=:native` is the tier to refine with.** `:full` buys
> extra physics at a fixed velocity-error floor.

> ⚠ **This is the IDENTICAL fingerprint to a genuinely wrong operator** — `e_u` pinned at a constant
> while `e_η` holds its optimal rate. **The MMS cannot distinguish "wrong operator" from
> "deliberately approximated operator"; only knowing what the solver implements can.** Anyone
> re-running `:full` and reading `p_u = 0.00` as a bug will be chasing a defect that does not exist.
> **Do not enable `:full` as a rate gate, and do not "fix" it.**

---

## 5. The Jacobian oracle, and the three defects it and the MMS found

**`test/test_jacobians_ad.jl` (17/17 gates over 8 models)** assembles the hand `∂R/∂u` and `∂R/∂u̇`
and compares them **entry by entry** against AD of the *same* residual, on a **sloping** bed. It
gates the linear branch on **equality** (`0.000e+00`, bit-exact) and the nonlinear branch on how the
gap **scales with state amplitude** — the distinction that matters, because a vanishing gap is the
deliberate quasi-Newton choice while an `O(1)` one is a defect. Measured: `∂R/∂u̇` exact in all 8;
`∂R/∂u` vanishing at order 1.11–1.16. Gate A0 additionally pins the Gridap fork by `hasmethod`.

Three defects were found against the assembly invariant (`MODEL.md` §7) and the MMS. Each is fixed;
what is kept here is the **rule** each one established.

| defect | what it was | rule it established |
|---|---|---|
| `𝓐/𝓚` package assembled **twice** under `:linear, flat_bed=false` | two consumers for one classification row, each guarded by only *part* of its activation condition | **the assembly invariant** — one consumer per row, guarded by the conjunction |
| `𝓐/𝓚` package **missing from `∂R/∂u̇`** | prefactor `H·∇h` does not scale with the solution ⇒ an **`O(1)`** error in the effective mass matrix. Newton converged to the fixed point of the *wrong map* and stalled at `‖r‖=4.8e-8` | **"converging slowly" and "converging to the wrong thing" look identical in a solver log.** Tell them apart by amplitude scaling, not by raising the iteration budget — no budget rescues an `O(1)` error |
| nonlinear gravity branch **missing the `−η∇h` half** of its own IBP | the `H∇η` identity has two pieces; the linear branch had both, the nonlinear branch one | the shared half is now assembled **once, outside the branch** |

**Two inferences that are invalid, and were both made:**

1. **"AD converges where the hand Jacobian fails, therefore the residual is correct."** AD
   differentiates the *same assembled residual*, so its converging proves residual↔Jacobian
   **consistency**, never residual **correctness**. **AD is an oracle for the JACOBIAN, never for
   the RESIDUAL.** Only a forcing derived independently of `problem.jl` verifies the residual.
2. **"The suite is green, therefore the model is verified."** The gravity defect survived 21 green
   files. The term was absent from the residual **and** from its own Jacobian, so AD agreed with the
   hand Jacobian throughout and every self-consistency check cancelled the error identically.

**How the gravity defect was localised — 3-case bisection, not code reading:**

| case | config | `rate_u` |
|---|---|---|
| A | sloping bed, `flat_bed=false` | −0.001 |
| B | *same code path*, `a_b=0` so `∇h ≡ 0` numerically | 3.993, 3.998 |
| C | `flat_bed=true` control | **bit-identical to B** |

One variable differed between A and B ⇒ the defect was in the `∇h` **values**, not the code path.
(C ≡ B also independently confirms the `flat_bed` switch is sound.)

---

## 6. Method rules earned by this campaign

1. **A refinement study measures the rate of whichever error DOMINATES — verify isolation in BOTH
   directions before interpreting any slope.** A saturated slope and a genuinely wrong coefficient
   produce the *same* observable. **Guards come in pairs**: "halving `dt` must not move the spatial
   study" needs the mirror "refining the mesh must not move the temporal study".
2. **Read the pairwise rate SEQUENCE, not the fitted slope.** A least-squares fit over a drifting
   sequence looks like a rate and is not one; and coarse-end over-convergence can cancel fine-end
   saturation to give a plausible fitted number from a contaminated window. Four levels make drift
   visible; three hide it.
3. **Check error MAGNITUDE before trusting a fine-level high-order rate.** `Q4` in 1-D exhausts
   double precision (`e_u` 1.79e-11 → 1.14e-11 over the last refinement).
4. **`d = 1.0` makes `h`-weighting invisible** — multiplying by `h ≡ 1` is the identity. Use `d ≠ 1`.
5. **A flat-bed regression cannot test `∇h` code.** New bed-slope terms need a sloping-bed gate from
   the start.
6. **Gates that test one side prove one side.** Passing forcing gates say nothing about
   residual↔Jacobian agreement. For a linear problem, "Newton converges in one iteration" is a free,
   sharp Jacobian check.
7. **Single-point gates are fragile** — one gave a false FAIL because the sample sat exactly at
   `k_y y = π/2` where `cos(k_y y)` vanishes. Sample a grid; compare against typical magnitude.
8. **Two tests sharing a bathymetry function and a physics tier can still be different problems.** A
   reference-vs-reference comparison is meaningful only when *every* discretisation parameter
   matches; mismatched domains once manufactured a phantom 2.5 % sequential/distributed discrepancy.
9. **When an error message names a library type, that is where the bug SURFACED, not where it
   lives.** Before searching that library, print the type of *every* input to the failing
   expression.
10. **Test a diagnosis against a case it cannot explain, rather than looking harder where it
    points.** A one-line bisection refuted a two-day-old diagnosis in one run.

---

## 7. The nested-closure hazard (structural in this codebase)

> **In Julia, a nested function that assigns a name already local to an enclosing function ASSIGNS
> THE ENCLOSING VARIABLE.** It does not create a new one. Declare `local` in nested helpers.

This codebase is full of long functions with nested helper closures (`strong_residual_model`'s
`Nvec`/`Ψ`/`Lvec`, the residual's CellField helpers), so the hazard is **structural, not
incidental**. One real instance cost two days: nested closures assigned `H`, `ukx`, `uky`, `L`, `Nc`
— all locals of the enclosing function — so a `Dual` leaked into the outer `H` under
`ForwardDiff.gradient` and surfaced 40 lines away in an unrelated array store. The `local`
declarations that fix it are commented as load-bearing.

**Audit mechanically, not by eye.** Parse every `src/*.jl` for nested functions assigning a name
also assigned by their enclosing function. That audit found the one real instance and four benign
ones (loop-carried state in the time loops, a `Ref` mutation, and a parser artefact).
**Re-run it after adding any nested helper.**
