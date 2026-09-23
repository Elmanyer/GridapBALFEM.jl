# OPEN_ISSUES.md — open work, each with its decisive next step

> **Scope — what belongs here.** *Known defects and gaps in work already done*, each with its
> decisive next step. Answers **"what is known to be wrong, missing, or unexplained?"**
>
> **What does NOT belong here:** work that is simply not started → [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md)
> (a run nobody has launched is a plan, not a gap); what is proven → [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md);
> test coverage → [`TEST_SUITE.md`](TEST_SUITE.md).


Ordered by whether the next step is *decided* or still needs a judgement call. Nothing in the solver
is half-built; every item below is a gap in verification, performance, or follow-through.

Related: [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md), [`TEST_SUITE.md`](TEST_SUITE.md),
[`CONFIGURATION.md`](CONFIGURATION.md).

---

## 0. ✅ Grid-scale instability in the fully nonlinear model — RESOLVED 2026-09-06

**It was the equal-order `Q2/Q2` horizontal pairing, not the model.** `η` enters momentum
undifferentiated via `∇·v`, so it plays the pressure role of a Stokes system and equal-order
continuous spaces are inf-sup deficient, admitting a spurious checkerboard at `λ ≈ 2·dx` — exactly
where the measured growth spectrum peaked. On Taylor-Hood the mode does not exist and the refinement
signature **inverts** (settled `η` 0.11055 → 0.10593 → 0.10481 as `dx` halves, converging on the
delivered 0.102, where refinement used to move the blow-up from t≈24 s to t≈10 s).

The fix cost **no residual change and no re-verification** — `Q2/Q1` and `Q3/Q2` are the pairings the
MMS campaign already certifies. It is now enforced by `check_taylor_hood` (`CLAUDE.md` rule 2b), so
the configuration cannot be selected again. The complete account — ten refuted hypotheses, the
campaign results across three meshes and three integrators, and the method lessons — is **`CLAUDE.md`
rule 12b**, which absorbed the former `NONLINEAR_INSTABILITY.md`.

⚠ **The quantitative phenomenology is void, not superseded**: the `σ ∝ A⁴` rates, the "no amplitude
threshold" claim, the tier ordering and the `dt`-masking numbers all describe a discretisation no
longer in use.

**What this leaves open, and it is real work:**
* **The 17 physics/smoke tests need re-baselining onto Taylor-Hood.** They call `setup_and_run`
  without an explicit `p_eta`, so they moved from `Q2/Q2` to `Q2/Q1` when the default changed. Their
  pinned constants were measured on equal order and must be **re-measured, not re-thresholded**.
  ✅ The verification tier is unaffected — it always passed `p_eta` explicitly.
* **A long-duration nonlinear regression is still missing.** No test runs long enough to have caught
  this (the mode needed 50–80 s to emerge); that gap is why it reached production. The Taylor-Hood
  80 s flume run is the natural basis for one.
* **Production runs made on equal order should be re-examined**, since they were made on a
  discretisation now known to be deficient.

---

## 0b. 🔴 THE NONLINEAR MODELS LOSE AN ORDER IN `p_η` — narrowed to the advection block

> Found 2026-09-11 by Campaign C, cause narrowed 2026-09-12. Letter-suffixed so the ten
> sections below keep their numbers (the §4 renumbering hazard, in miniature).

**THE OBSERVATION.** Campaign C, 1-D, Q3/Q2, `a_eta = 0.8`, five-level ladder `nx = 4…64`,
`nl_tol = 1e-14`. Same basis, same ladder, same mesh; only `regime` differs:

| | pairwise `p_η` | `e_η` ratios (8.00 = optimal 3rd order) |
|---|---|---|
| **linear** (models 1, 2) | 2.990, 2.997, 2.999, **3.000** | 7.95  7.99  8.00  8.00 |
| **nonlinear** (models 3, 4, 5, 6) | 2.971, 2.922, 2.757, **2.450** | 7.84  7.58  6.76  5.46 |

⚠ **This is an ORDER REDUCTION, not a floor.** `e_η` keeps falling — it is the *rate* that decays,
from 3 toward 2, and it decays MONOTONICALLY WITH REFINEMENT. That is the signature of a fixed
lower-order error component overtaking the third-order one (rule 33), and it means the problem gets
worse on better meshes, not better.

⚠ **`p_u` moves the OTHER WAY on the same runs** — 2.77 → 3.19 → 3.58 → 3.69, rising toward its
optimal 4 while `p_η` falls. Any explanation has to account for both signs at once; several
otherwise-plausible ones do not.

**WHERE IT DOES AND DOES NOT SHOW, and why that is consistent.**

| pairing | optimal `p_η` | measured | reading |
|---|---|---|---|
| Q2/Q1 | 2 | 2.000 exactly | a lower-order component at order 2 CANNOT be seen — it is degenerate with the optimal rate |
| Q3/Q2 | 3 | 2.450 and falling | visible, because the third-order part decays past it |
| Q4/Q3 | 4 | 4.000–4.002 | ladder stops at `nx=32`; discretisation error still dominates |

So the absence at Q2/Q1 and Q4/Q3 is expected under the same explanation, not evidence against it.

**WHAT IT IS NOT — five eliminations, each with its own measurement.**

| ruled out | evidence |
|---|---|
| algebraic contamination (rule 32) | Newton reaches **~1e-15 at every level of every model** (`nl_tol=1e-14`); from the per-level `diagnostics.csv` the §4b work added |
| the `𝓝` pressure blocks | `:none` vs `:native` agree to **4 digits** — m3 2.4498 vs m5 2.4499, m4 2.4494 vs m6 2.4496 |
| bed-slope (`∇h`) code | flat vs variable bed agree to 4 digits; reproduced on an independent process (m4 control, 2.4494) |
| the linear core | linear models are textbook 3.000 on the identical basis, ladder and mesh |
| **quadrature under-integration** | refuted quantitatively — see below |
| the test harness generally | the `quad_extra` control reproduced Campaign C to **5 significant figures**, so the driver edit perturbed nothing |
| **a 1-D posing artefact** (`ny=1`, `:wall`, `ky=0`) | **reproduces in full 2-D**: `T9_tier2_2d` model 3 measured `p_η = 2.5565` on the `nx = 8…64` ladder, a separate campaign months earlier, against 2.450 in 1-D. `HORIZONTAL_CONVERGENCE.md` §2c |
| a 2-D-only or dimension-dependent effect | the **absence** reproduces too: at Q4/Q3 the nonlinear models are optimal in BOTH dimensions (2-D 3.9997–4.0004, 1-D 4.000–4.002, 6/6 studies). The defect is present at Q3/Q2 and absent at Q4/Q3 **in 1-D and 2-D alike** |

**THE QUADRATURE HYPOTHESIS, AND WHY IT FAILED — worth keeping, because the shape of the refutation
is the lesson.** The default MMS quadrature is `2·max(p_u,p_η)+2` = degree 8 at Q3/Q2, which in
Gridap is a 5-point Gauss rule **exact to degree 9** (measured; do not read exactness off the
`Measure` argument). Counting from `problem.jl`, per direction: the linear core is ≤ 8 and exact,
while the nonlinear advection `H·(𝓜₃:(U⊗∇U))·W` (line 149) and the nonlinear pressure `H²·(sP·DW)`
(line 110) are **degree 10** — genuinely under-integrated. `src/mms_driver.jl` hard-coded the degree,
so no MMS study had ever been able to vary it; `quad_extra` is now threaded through
`run_conv_study → run_mms_case → Measure`.

Assembling the residual directly at degrees 8/10/12/16 (`output/local/mms/quadrature_test/
knob_liveness.log`):

```
LINEAR    core  ‖R‖ changes 4e-14 across ALL degrees        -> exact at the default
NONLINEAR core  ‖R‖ changes 3.49e-08 from degree 8 -> 10,
                then 4e-14 for 10 -> 12 -> 16               -> under-integrated, exact from 10
```

So the knob is **live** and the crime is **real**. It is also **six orders of magnitude too small**:

| | |
|---|---|
| quadrature error in the residual | 3.5e-08 relative |
| its effect on converged `e_η` | ≤ 4.7e-07 relative ≈ **4.6e-13 absolute** at `nx=64` |
| change needed to restore 3rd order at `nx=64` | `e_η` 9.77e-07 → 6.67e-07 = **3.1e-07 absolute** |

✅ **COMPLETED 2026-09-13 — the full dose-response, both bed types, all identical:**

| case | `quad_extra` | exact to | `p_η` | `e_η` ratios |
|---|---|---|---|---|
| m3 flat, control | 0 | x⁹ | 2.450 | 7.84 7.58 6.76 5.46 |
| m3 flat | **2** | x¹¹ | 2.4498 | 7.840 7.578 6.758 5.463 |
| m3 flat | **4** | x¹³ | **2.450** | 7.84 7.58 6.76 5.46 |
| m4 varbed, control | 0 | x⁹ | 2.4494 | 7.840 7.577 6.756 5.462 |
| m4 varbed | **4** | x¹³ | **2.4494** | 7.840 7.577 6.756 5.462 |

The control reproduced the Campaign C baseline to 5 significant figures, so the driver edit
perturbed nothing. Max relative change in `e_η` between `q=0` and `q=4`: **1.422e-08** — the knob is
live and moves the answer, by **seven orders less** than the 3.1e-07 the rate deficit requires.
`quad_extra=2` already covers the degree-10 integrand and changes nothing, which is the shape a
genuine quadrature cause could not produce.

⚠ **Two method traps fired during this probe and both are worth remembering.** (a) The first
liveness check integrated `x⁹`, which a degree-8 rule already does exactly — a check that *could not
fail*, and it passed a knob whose reach was still unproven. (b) The `e_η` columns came back BITWISE
identical, which reads as "dead knob, result void" — I called it that — when it was really a live
knob whose effect fell below the CSV's 7-significant-digit print precision. **Neither a liveness
check nor a null result means anything until the effect size is measured against the effect being
explained** (rule 12b's standing lesson, and the same shape as the skew-advection refutation).

**WHAT REMAINS, AND THE DECISIVE NEXT STEP.** By elimination the error sits in the **nonlinear
advection block** (`problem.jl` §"nonlinear advection (𝓜/𝓖 block)", lines 139–149) — a second-order
consistency error in `𝓜₃/𝓖` as discretised. The next step is **derivation, not another run**: check
that block term-by-term against `BALFEM_models/` (`GlobalResidual`, and the `𝓖` identity recorded in
`CLAUDE.md` under the reverted skew-advection work), since the LaTeX is the single source of
mathematical truth and a residual term it does not contain is one nobody can check. Two cheap
discriminators exist if a run is wanted later:

* **amplitude scaling** — a consistency error in the advection block scales with the manufactured
  amplitude. Re-run model 3 at `a_eta = 0.4` and `0.2`: if the crossover moves to finer mesh, the
  offending term is quadratic in the solution; if the rate degrades identically, it is not.
* **extend Q4/Q3 to `nx = 64`** — it should then show the same decay once its fourth-order part has
  decayed past the second-order component. ⚠ Check `e_u` against the ~1e-10 algebraic floor first
  (`PLANNED_CAMPAIGNS.md` §0); at Q4/Q3 that level may be unmeasurable.

⚠ **IT WAS ALREADY IN THE 2-D DATA AND WENT UNSEEN.** The 2-D measurement above predates Campaign C
by months. It was not noticed because the old CSV schema persisted only the fitted slope and the
finest error — and a *fitted* 2.78 reads as "slightly low", where the pairwise sequence
`2.971 → 2.922 → 2.756 → 2.450` reads as a rate collapsing under refinement. The per-level sequence
was added in `fdc4774`. **The instrument, not the physics, is why this is a 2026-09 finding**, which
is the concrete cost of rules 32/33 having been unenforceable in the data format.

**Scope of the damage.** This affects *measured convergence rates of the nonlinear models*, not the
verified scope as stated: `VERIFIED_SCOPE.md` records 30/30 spatial studies at optimal order, and
those were run at Q3/Q2 with **four levels to `nx=32`**, where `p_η` still reads 2.76–3.00. The
degradation is a **fine-mesh** effect that the shorter ladder never reached. The claim "the `:none`
and `:native` models are verified" therefore stands *as far as it was measured*, but it should now
carry the qualifier that nonlinear `p_η` degrades beyond `nx = 32` at Q3/Q2.

---

## 0c. ⛔ THE CLASS-III (`nl_pressure=:full`) NONLINEAR INSTABILITY — complete account

> The authoritative record for this topic (`INDEX.md`). Written 2026-09-17, after the Yang & Liu
> collapse verification closed off two of the three candidate causes.
> ⚙ **AMENDED 2026-09-23 after the mixed-formulation campaign** (`NEW_TREATMENT.md` Part H,
> `CLAUDE.md` §5.2e): §8's option table is **re-scored** and §9 is **rewritten**. Read those two
> before acting on anything in §1–§7, which pre-date the campaign.

### 1. The short version

`nl_pressure=:full` diverges on a **flat bed** in 8 wave periods at `A = 0.10` m, while `:native` —
identical in every other knob — completes 100 periods. Three things could have caused it: a wrong
**operator**, an incomplete **Jacobian**, or the **assembly** of the Class-III blocks. The first two
are now eliminated by measurement. **What remains is the assembly**, specifically the frozen `L²`
projections used to carry the `{1,2,4,5}` blocks.

⚙ **AMENDED 2026-09-23 — THE ASSEMBLY IS CONFIRMED AS *A* CAUSE BUT NOT *THE* CAUSE, AND THE
THREE-WAY SPLIT ABOVE WAS NOT EXHAUSTIVE.** Removing the projections entirely (the mixed unknown
`𝖦 ≈ ∇𝖲`) fixes the Q2/Q1 flat case — **≥61 periods against 12.6 s** — and beats the projected
control in all eight matched pairs, so the assembly was real and large. **But `:full` still dies at
Q3/Q2, and still degrades with `dx`, in a formulation that has no projection and forms no second
derivative.** So a fourth possibility, absent from the original list, is now live: the **coupled
discrete system** itself — in particular the new `𝖦↔u` pairing, which sits at **equal order**
(`src/horizontal.jl:59`) and has never been analysed. ⚠ And the Class-III attribution at Q3/Q2 is
**not established at all**, because no `:native` Q3/Q2 control has been run (§9, step 1).

### 2. Why the operator is not at fault — the literature check

Our `p = 1` reduction was verified term by term against Yang & Liu (2024), whose non-linear LFE-*M*
is stable and reproduces laboratory physics (sideband instability, harmonic generation over a bar).
Full account: `VERIFIED_SCOPE.md` §0b, `CLAUDE.md` §5.2d, and chapter 4 of
`latex_docs/BALFEM_models/`.

| object | source | worst relative difference |
|---|---|---|
| continuity weights `Φ` | their ⌊2.31⌋ | **0** |
| linear coefficients `A`, `B`, `D` | their Appendix A | **7.8e-15** |
| vertical velocity `w` (non-linear) | supplementary §A | **1.0e-14** |
| non-hydrostatic pressure `p_nh` | supplementary §B | **1.7e-15** |
| **weighted momentum residual** | **supplementary §C + ⌊2.28⌋–⌊2.30⌋** | **2.1e-12** |

⚠ **Two of these test physics, not agreement**: the pressure stage against the vertical momentum
equation, the residual stage against the projection of ⌊2.12⌋ — so an error *shared* by both
derivations would still have been caught. The last stage uses the **solver's own assembled tensors**,
so `src/vertical.jl` is verified too, not just the LaTeX.

**Conclusion: the equations are right, and so is the code that builds the vertical tensors.**

### 3. Why the Jacobian is not at fault

See rule 17b. Every cell of the `dx`×`dt` factorial was run with both the hand Jacobian (which omits
`{1,2,4,5}`) and the exact AD Jacobian: **onset identical or one output interval apart, `u_max` at
failure agreeing to three significant figures, only the iteration count differing.** Completing
`jacobian_u` would buy cost, not stability.

### 4. What "Class III" means, and why those terms exist

The classification is in `GlobalResidual.tex` §`sec: pressure operator implementation`, and it follows
from one ground rule: **only first derivatives of the unknowns are admissible in a `C⁰` Lagrangian
integrand.** Sorting the eight `𝓝ₖⱼ` components by what they demand:

| class | components | demand | status |
|---|---|---|---|
| **I** | 6, 7, 8 | products of first-order fields only | directly implementable |
| **II** | 3 | one second derivative, but it lands on the **analytic** bed Hessian `∂²h` | admissible (ground rule 2) |
| **III** | **1, 2, 4, 5** | second derivatives of the **unknowns** | ⛔ not directly implementable |

Expanding Class III, the demand reduces to two objects:
* `∇sₖ` where `sₖ = ∇·(Huₖ)` — contains the **velocity Hessian** `∂²u` and, through `∇H`, the
  **free-surface Hessian** `∂²η`. Appears in components 1, 2, 5.
* `∇bⱼ` where `bⱼ = uⱼ·∇H` — its only inadmissible part is `∂²η`. Appears in component 4.

⚠ **The count is smaller than four terms, and smaller than the document currently states.** Two
identities collapse it:
* `𝓝₂ = −𝓝₁ + sₖ∇·uⱼ`, so regrouping `𝓝₁Θ₁ + 𝓝₂Θ₂ = 𝓝₁(Θ₁−Θ₂) + (sₖ∇·uⱼ)Θ₂` moves part of the pair
  into admissible territory and leaves the rest on a **reduced** weight;
* `𝓝₅` is the `(k,j)` transpose of `𝓝₁`: `Σₖⱼ(−uₖ·∇sⱼ)T'ᵢₖⱼ = Σₖⱼ(−uⱼ·∇sₖ)T'ᵢⱼₖ` after relabelling,
  so it folds into the same contraction with weight `(Tᵢₖⱼ + T'ᵢⱼₖ)`.

**Net: one `∇s` contraction and one `∂²η`.** That is the whole of Class III.

⚠ **And `∂²η` is the Hessian of the field in the LOWER Taylor-Hood space.** At **Q2/Q1** — the pairing
every instability run has used — `η` is piecewise linear, so in 1-D `∂²η ≡ 0` **identically**:
component 4's free-surface half, and the `∇H` half of `∇s`, contribute nothing there at all. On each
element `H·u` is (Q1 × Q2) = cubic, so the exact `∂²(Hu)` is piecewise **linear** and discontinuous
across every face; at Q3/Q2 it is piecewise cubic and `∂²η` is no longer identically zero.
⚠ **CORRECTED 2026-09-18: the object Class III demands is a SECOND derivative of the unknowns, not a
third.** `s = ∇·(Hu)` is first order and `∇s` is second — `GlobalResidual.tex` classifies Class III as
"second derivatives of the unknowns" and that is right.
⚠ **MEASURED 2026-09-23, AND THE EARLIER CORRECTION WAS ITSELF TOO STRONG.** Both derivative counts
are right, of different objects. Probing each quantity with a perturbation `ε(x−x₀)^k` — which
nulls every derivative below order `k` at `x₀` — gives, at `x = 3.7`:

| object | `∂²u` | `∂³u` | `∂⁴u` |
|---|---|---|---|
| `w` | independent | independent | independent |
| `p_nh` | **depends** | independent | independent |
| their §C residual ⌊2.30⌋ | depends | **depends** | independent |
| our residual, tensor form | depends | **depends** | independent |

So **the model does carry `∂³u` in its strong-form momentum equation** — `p_nh` carries the Hessian,
momentum carries `∇p_nh` — and Yang & Liu's structural claim is correct. Our `∂²` statement is a
claim about the **weak form**: `R_P` is assembled as `∫H²(𝓝∴𝓟)·(∇·v)`, so the integration by parts
that produced it has already moved one derivative onto the test function, leaving the trial field
carrying at most `∂²`. **Class III is a Hessian problem only because of that IBP**, not because the
model is lower-order than the literature says.

### 5. Why projection was chosen in the first place

The admissibility problem is not a property of the components alone — it depends on the **prefactor**
of the block each one sits in. `𝓝` enters the residual three times:

| block | prefactor `Ψ` | can integration by parts repair it? |
|---|---|---|
| (i) non-linear pressure, **bed slope** | `H ∇h · v` — prescribed, analytic | ✅ **yes** — `∇Ψ` produces only admissible objects. Vanishes identically on a flat bed |
| (ii) non-linear pressure, **surface slope** | `H ∇H · v` — contains the unknown `η` | ⛔ **no** — `∇Ψ ∋ ∂²η`; IBP trades `∂²u` for `∂²η`, both inadmissible, no cancellation |
| (iii) **leading pressure** `R_P` | `H² (∇·v)` — the test function already carries a derivative | ⛔ IBP would put `∂²v` on the **test** function |

With IBP unavailable for (ii) and (iii), two sound options remained: a **mixed formulation** (promote
`∇s` to a genuine unknown with its own weak equation) or a **projection**. The projection was chosen
because the mixed formulation enlarges the system substantially, against the goal of a fast solver —
**and on the expectation that the projection error would be an accuracy cost, not a stability one.**
That expectation is what the measurements contradict.

### 6. What the measurements actually implicate

`CLAUDE.md` §5.2c. At fixed `dx`, onset moves **out** as `dt` falls; at fixed `dt`, it moves **in** as
`dx` falls.

| `dt` (at `dx` = 0.25) | 0.04 | 0.02 | 0.005 |
|---|---|---|---|
| onset | 12.60 s | 21.20 s | **40.40 s** |

⚠ **These are two different error sources and they should not be conflated:**

| | mechanism | evidence |
|---|---|---|
| **recovery error** | the recovered `∇s` ≠ the true `∇s`. The true `∂²` of a `C⁰` field is **not an `L²` function**: it is `{∂²u}` cell-wise **plus a Dirac layer on the skeleton** weighted by `[∂ₙu]`. Projecting `s` and then differentiating replaces that layer with a smoothed surrogate | `dx`-refinement **advances** onset (×2.5) |
| **lag error** | the recovered field is from the **previous step** — `update_nlp_state!` runs *after* a step, so the residual is evaluated against a stale state | `dt`-refinement **delays** onset (×1.7 per halving) — the stronger, cleaner trend |

**The lag is not intrinsic to projecting. It is intrinsic to *freezing*.**

⚠ **AND THE RECOVERY ERROR HAS BEEN MEASURED DIRECTLY, ON AN ANALYTIC SOLUTION — it did not have to
be inferred from a crash time.** The `T8_full_projection` studies ran models 7 and 8 at
`a_eta = 0.4` on four vertical bases (`PLANNED_CAMPAIGNS.md` §3), an amplitude at which Newton
converges to ~2e-09 and the Jacobian is therefore discharged (rule 17b), so the floors are the
**residual's own** error:

| basis | `p_η` (optimal 3) | `p_u` (optimal 4) | `e_u` floor |
|---|---|---|---|
| P1LFE-2 | 2.56 | 1.63 | 1.4e-06 |
| P1LFE-3 | 2.45 | 1.88 | 5.4e-06 |
| P1LFE-4 | 2.41 | 2.10 | 9.4e-06 |
| P2LFE-1 | 2.57 | 2.65 | 2.4e-07 |

**One to two orders short in every field on every basis, and the floor grows with `Nσ`** — with the
number of `∇s` objects that have to be recovered. Independently, the single-study P1LFE-2 model 7
figure after the `nlp`-context fix is `p_u = 1.948` (`VERIFIED_SCOPE.md` §4). **The projection is
not a small perturbation of the exact operator at production resolution: it is the dominant error in
the velocity field, and it converges too slowly to be refined away.**

### 7. ⚠ A correction to the original reasoning

`sec: pressure operator implementation` rules out IBP for block (iii) because it "would place a second
order derivative on the test functions — inadmissible on `C⁰` test spaces for exactly the same
distributional reason as for the trial fields."

**That equivalence is too strong.** The distributional objection concerns the *trial* space — whether
the discrete solution is well defined and convergent. A **test** function is a chosen, known object:
`∂²v` is a computable cellwise polynomial (non-zero from Q2 up), discontinuous only across faces.
Keeping second derivatives cellwise on both sides and repairing the face jumps with penalty terms is
the standard method for fourth-order problems (`C⁰` interior penalty).

⚠ **And one of the original grounds no longer binds**: Gridap exposes `∇∇` — we already use it for the
analytic bed Hessian at `src/nlpressure.jl:35` — and full skeleton machinery
(`SkeletonCellFieldPair`, jump/mean). **This was never a missing-API problem.** `∇∇` returns the
**broken (cell-wise) Hessian**, which is a perfectly good object — but it is *not* the distributional
`∂²` of a `C⁰` function, which carries a Dirac layer on the skeleton that cell quadrature never sees.
⚠ **AND RAISING THE POLYNOMIAL ORDER DOES NOT CHANGE THAT.** `FESpace(model, ReferenceFE(lagrangian,
…, p); conformity=:H1)` is exactly `C⁰` for **every** `p` — nodal DOFs match function values across a
face, never normal derivatives. `Q2`, `Q3`, `Q4` are all `C⁰`; only a different element family
(Hermite, Argyris, splines/IGA) is `C¹`. **The task is to choose a formulation that is consistent with
a broken `∂²` — not to find a library or an order that makes `∂²` classical.**

### 8. Options, viability and complexity

> ⚙ **RE-SCORED 2026-09-23 AFTER THE MIXED CAMPAIGN** (`NEW_TREATMENT.md` Part H, `CLAUDE.md`
> §5.2e). **Five of the six verdicts below have changed and the table as written is superseded** —
> it is kept because the *reasoning* in each row is still the right way to think about the choice.
> The outcomes:
>
> | # | option | original verdict | **what actually happened** |
> |---|---|---|---|
> | 0 | run at Q3/Q2 | ✅ do first | ⛔ **DONE, REFUTED, SIGN REVERSED** — Q3/Q2 is *uniformly worse*, all four cells, and worse even at matched DOF count |
> | 1 | de-lag the recovery | ✅ cheapest test | ✅ **DISCHARGED** by the mixed construction (it solves `𝖦` at the current iterate). Both `dt` arms improved — the lag was real |
> | 2 | algebraic reduction | ✅ do regardless | ✅ **DONE, EXACT** (4.4e-16; `test_class3_reduction` 19/19, parity 4/4) |
> | 3 | `C⁰` interior penalty | the principled route | ⛔ **RATIONALE REMOVED** — it fixes a broken `∂²`, and the `dx` signature it was conditioned on survived a formulation with **no `∂²` to break** (H4). Do not launch on the old justification |
> | 4 | **mixed formulation** | ⛔ **rejected (system size)** | ⚙ **BUILT ANYWAY, AND IT IS THE SOURCE OF EVERY RESULT ABOVE.** The rejection was wrong to treat as final: it bought the first long-running `:full` configuration in the project (Q2/Q1 flat, ≥61 periods vs 12.6 s). Its cost objection stands and then some — Q3/Q2 quoted **186–242 h** ETAs — so it is a **diagnostic instrument, not a production path**, unless the auxiliary rows get hand Jacobians |
> | 5 | `C¹` / IGA | high, no element | untouched |
>
> ⚠ **The generalisable lesson (rule 39b): option 4 was rejected on cost and became the only
> experiment that could refute option 3's premise.** A construction that *cannot commit* the
> suspected error is worth building as a discriminator even when it is unaffordable as a solver.


| # | option | complexity | grows the system? | addresses | viability |
|---|---|---|---|---|---|
| **0** | **Run `:full` at Q3/Q2** | **trivial** (one env var) | no | nothing — a diagnostic | ✅ do first: every run so far was Q2/Q1, where `∂²η ≡ 0` identically and `∂²(Hu)` is only piecewise linear |
| **1** | **De-lag the recovery**: evaluate it inside the Newton loop from the current iterate instead of freezing it from the previous step | **low** — a change to *when* `update_nlp_state!` is called | no | the **lag** half | ✅ cheapest test of the mechanism the `dt` ladder implicates. Costs one mass-matrix solve per residual evaluation; that matrix is constant, so factor once and back-substitute |
| **2** | **Algebraic reduction** via the two identities of §4 | **low** | no | shrinks what any later treatment must carry | ✅ worth doing regardless; does not remove the problem |
| **3** | **`C⁰` interior penalty** on blocks (ii) and (iii) | **moderate–high** | no | the **recovery** half, principled | the only genuinely projection-free route that keeps the system size. Needs a penalty parameter and its own stability argument |
| **4** | **Mixed formulation** (`G ≈ ∇s` as an unknown) | moderate | ⛔ **yes, substantially** | both halves | ⛔ **rejected** — against the goal of a fast solver |
| **5** | **`C¹` space** (Hermite / splines / IGA) | high | no | both halves | no standard `C¹` Lagrange element in Gridap; a new element family |

### 9. Current state of the solver, and the next steps

*Rewritten 2026-09-23 after the mixed-formulation campaign — `NEW_TREATMENT.md` Part H.*

**State.** `:none` and `:native` are sound: `:native` flat-bed completes 100 wave periods at
`A = 0.10` m, Newton flat at 5.12. **`:full` is no longer unconditionally unstable** — with the
Class-III projections replaced by the mixed unknown `𝖦 ≈ ∇𝖲`, the Q2/Q1 flat case runs **≥61 wave
periods flat on the `:native` trace** where it used to die at 12.6 s, and beats its projected control
in **all eight** matched pairs. ⚠ But the cause is **only half identified**: every **Q3/Q2** arm
still dies, earlier than its Q2/Q1 twin and earlier even at matched DOF count, and the **`dx`
signature survives a formulation that never differentiates a `C⁰` field twice**. ⚠ Separately, **any
variable bed** grows a lee-shoulder mode at a rate set by `|∇h|` (§0d) — a different defect, not
addressed by any of this.

⚙ **FINAL CAPTURE 2026-09-23 15:29 — THE Q2/Q1 RESULT SPLITS ON AMPLITUDE.** `A` = 0.10 reached
**t = 100.6 s (62.9 periods)** flat in the fourth decimal with Newton constant at 6 — a genuine
post-fill stability result. **Its `A` = 0.15 twin diverged at ~87 s** with the same velocity-led
signature (η 0.207 → 0.278, u 0.851 → 1.144, Newton 8 → 15 in one window). **So the mixed
formulation DELAYS the instability rather than removing it**, and amplitude is a **third independent
axis** beside pairing and `dx`. ⚠ Quoting the `A` = 0.10 trace alone yields a confident and wrong
"fixed" verdict — it took the same-batch amplitude twin to see the mode was merely late.

**What the campaign settled.**
* ✅ The frozen projection was a **large real contributor** (8/8 matched pairs, 1.9×–>16×).
* ⛔ It was **not the whole cause** (Q3/Q2, and the surviving `dx` sign).
* ⛔ The **broken-Hessian recovery explanation of the `dx` signature is REFUTED as stated** — it was
  measured in a path that cannot commit that error. No replacement explanation yet; three candidates
  remain unseparated (H4: it was never the recovery error / `𝖦`'s own approximation error / a
  grid-scale problem in the new `𝖦↔u` pairing).
* ⛔ Raising the polynomial order is **not** "more physics" and **not** a one-variable change
  (rule 2d) — it moves representable content *and* effective resolution together. A low-order pass
  can be a pass on a **partially masked operator**, which is what Q2/Q1 turns out to be here.

**Next steps, in order.**
1. ⛔ **`:native` at Q3/Q2 — DO THIS FIRST, IT GATES EVERYTHING.** Every Q3/Q2 arm run is `:full`
   and the only `:native` runs are Q2/Q1, so nothing in this batch can distinguish "Class III is the
   carrier" from "Q3/Q2 is unstable here for an unrelated reason" (rule 14c). ⚠ This is the **second
   time on this branch** the control was the thing not run — see `NEW_TREATMENT.md` F.0, where the
   omission voided an entire harness campaign.
2. **`𝖦` one order below `u`.** `Vaux` is built from `reffe_U` (`src/horizontal.jl:59`), so `𝖦` is
   at **equal order with the velocity**. Its own block is a Gram matrix and coercive — not the
   classic rule-2b failure — but the *coupled* pairing has never been analysed, and this is the only
   probe of H4 reading (3). Open since F.3.2, still never run.
3. **Q3/Q2 at nx = 120**, to finish H3's isolation by matching `dx` and not only DOF count.
4. **Hand Jacobians for the auxiliary rows** — at 186–242 h ETAs the Q3/Q2 mixed tier cannot
   produce a stability claim at all. A prerequisite, not an optimisation.
5. **Let `c3v_base_mixed` reach 100 periods**; it is the only candidate reference trace for a
   `:full` long-duration regression gate (`CLAUDE.md` §5.7 item 4).
6. ⛔ **NOT `C⁰`-IP**, until someone re-derives a justification for it (see the re-scored §8).

⚠ **What is still NOT established.** *(The loose end below is unchanged by this campaign and is now
joined by a second.)* The mode pins at the **inflow** in runs that die mid-fill and near the
**front/sponge** in the one that died after fill completion — a fill-state explanation is plausible
but remains a hypothesis, and both VTK series are on disk. **And now:** with the recovery
explanation refuted, there is **no standing mechanism** for the `dx` signature at all. `:full`
having a stable configuration does not mean it is understood.

---

## 0d. ⛔ ANY VARIABLE BED GROWS A LEE-SHOULDER MODE; `|∇h|` SETS THE RATE

> Found 2026-09-15, same campaign. Full context: `CLAUDE.md` §5.2b.

Five runs identical but for the bar's shoulder length — height 2.0 m on `d` = 3.5 m, span 26–34 m,
crest depth 1.5 m, `max|∇h| = hbar/(2·sramp)` — all `:native`, all `A = 0.10` m:

| shoulder | `max\|∇h\|` | face | onset | periods |
|---|---|---|---|---|
| 0.5 m | 2.0 | 63° | t = 37.2 s | 23 |
| 1.0 m | 1.0 | 45° | t = 57.8 s | 36 |
| 1.5 m | 0.67 | 34° | t = 93.8 s | 59 |
| 2.0 m | 0.5 | 27° | t = 137.2 s | 86 |
| flat | 0 | — | **none** | ✅ 100 |

One mode, four growth rates. Identical fingerprint in all four: growth **pinned at the downwave
shoulder** (x ≈ 32–34.5 m — rule 40, never read growth without `x_at_max`; a pinned location is a
stationary mode, not a travelling wave), `u_max` 2.2–3.3 against the flat control's steady 0.42, η to
0.8–3.0 m from a 0.10 background, Newton degrading 12 → 15 → 20 → cap.

**Onset falls monotonically with slope, and the flat bed is the limit point of the same family** —
a *rate*, not a threshold. A badly-posed bathymetry would give a threshold; a rate is the signature
of an instability in the formulation. The `t = 93.8` point was a **prediction** (it had to land
between 58 and 137) and was met.

⚠ **THE FIRST READING OF THIS WAS WRONG AND THE ERROR IS INSTRUCTIVE.** The 1:2 bar was called
stable at t = 80 s — it died at t = 137 s. *A "threshold" may just be a run that ended too early*
(rule 12b's standing lesson) applies to the shoulder ladder exactly as it applied to the equal-order
mode. Any bar result quoted before 100 periods is provisional.

**DECISIVE NEXT STEP — NOT YET RUN.** Halve `dx` on a bar case. Rule 38b separates the two whole
classes of explanation in one run: refinement **delaying** onset ⇒ under-resolution of a steep bed;
refinement **advancing** it ⇒ a grid-scale problem in the ∇h terms, which the pinned location and
velocity involvement already suggest.

---

## 1. 🔴 Cluster memory attribution

`--mem-per-cpu=4G` is in every launcher and is **required**: a `rome` job at the node-default
2 GB/core was OOM-killed with the same error as the earlier crashes. So the compile-spike
explanation is not sufficient on its own. ***Why* is the open question, not *whether*.**

Four hypotheses, one decisive measurement each:

| # | hypothesis | decisive measurement | verdict if true |
|---|---|---|---|
| H1 | per-rank JIT (image absent / stale / incomplete) | **when** the RSS peak occurs — a compile peak is **before step 1**; confirm `-J` was on the command line and the freshness check did not warn | rebuild the image, re-test 2 GB/core |
| H2 | GMRES cache over-allocation | was the job's commit after the `krylov_m` fix? predicted ≈5.6 MB/rank at `krylov_m=100` | if pre-fix, re-run post-fix |
| H3 | per-step leak (Gridap caches, `:full` frozen projections, VTK buffers) | RSS rising with step index | escalate; bisect by disabling `write_w`/`write_pressure`, then `nl_pressure` |
| **H4** | **baseline footprint simply exceeds 2 GB/core** | RSS flat and already > 2 GB after step 1 | **4 GB/core is correct and permanent** |

**H4 leads on local evidence.** On a *tiny* sequential case (3366–4669 DOFs) the first diagnostics
sample reads **1390 MB**, and a 400-step run went 2046 → 2089 MB (**+2 %, flat**). A Julia + Gridap +
GridapBALFEM process costs ≈1.4 GB before solving anything, so 2 GB/core leaves ~0.6 GB/rank for the
computation. Flat memory also means **H3 is not supported at this scale** — though a 32-rank
`:full` run exercises paths this case does not, so H3 is not eliminated.

⚠ Those are **local** numbers and must not be reported as the cluster answer: the cluster runs
against a ~1 GB sysimage (code mapped and shared) and MPI adds per-rank buffers.

> **NEW EVIDENCE FOR H3, 2026-08-22 — the flat-memory result does NOT hold at length.** The
> vertical-basis campaign ran 8-10 independent sequential processes for ~16 h. Every one grew
> steadily and none plateaued:
>
> | elapsed | typical process RSS |
> |---|---|
> | start | 1.5 GB |
> | +4 h | 2.0-2.5 GB |
> | +10 h | 2.5-3.1 GB |
> | +14 h | up to **3.9 GB** |
>
> Machine swap went 0 → 630 → 996 → 1820 → 2171 MB over the same period and **three processes had to
> be killed** to keep the run off swap. This is a *sequential* Julia+Gridap process — no MPI, no
> distributed assembly — so it isolates the growth to the per-process solve path rather than to
> anything rank-related.
>
> That does not confirm H3 by itself (a long-lived process accumulating Gridap caches and compiled
> specialisations is not the same as a per-step leak), but it **removes the "flat memory ⇒ H3 not
> supported" argument above**, which was measured over 400 steps on one tiny case. The decisive
> measurement is unchanged — RSS against step index within a single run — and it can now be taken
> locally, cheaply, from any of these shard logs, without cluster access.

**Evidence still to collect (needs cluster access):**

```bash
sacct -j <jobid> --format=JobID,JobName,State,ExitCode,MaxRSS,MaxRSSTask,AveRSS,Elapsed
seff <jobid>
grep -c '^\[dist\] step' <job>.out        # how far it got before the kill
grep -i 'sysimage\|stale\|-J ' <job>.out  # was the image used, was it fresh
```

Exit code **137** = OOM-kill; an `oom-kill event` line in `.err` confirms it.

**Decision rule.** Drop to 2 GB/core only if global-max RSS stays below ~1.6 GB/rank for a
full-length run of the largest small-domain case **and** the peak is not a pre-step-1 compile spike.
Any monotonic growth (H3) escalates instead — it means run length is bounded by memory, which no
per-core request fixes. **Cheapest decisive job:** `run/dist_small/run_lin_periodic_plane_small.sh`.

---

## 2. ✅ The MMS path never assembles the `:full` frozen projections — FIXED 2026-09-01

> ⚠ **THIS SECTION WAS STALE UNTIL 2026-09-15.** The decisive next step below WAS taken: both MMS
> drivers now build the `nlp` context and both time loops prime it from the initial condition.
> Measured on P1LFE-2 model 7, `e_u` fell **1.19e-03 → 1.37e-06 (~870×)** and `p_u` went from a flat
> −0.00 to **1.948** — so the recorded floor measured **omission**, as the last paragraph predicted,
> and `:full` is *not* "MMS-unverifiable by construction". `CLAUDE.md` carries the full account.
> ⚠ Consequence recorded there: with the blocks in the residual but still absent from `jacobian_u`,
> the quasi-Newton gap has a **cliff in amplitude** — tier-3 studies must drop `a_eta` to ≤ 0.4.
> ⚠ **And see §0c: completing `jacobian_u` does NOT make `:full` stable.** The two are separate.
> The historical text follows.

Found 2026-08-21 while designing the vertical-basis campaign. **Reported, deliberately not patched**
— it is a design-level question about an existing interface, not a usage error.

`run_time_loop` defaults `nlp = nothing` (`src/timeloop.jl:130`) and **`run_mms_case` never passes
one** (`src/mms_driver.jl:99–105`); only `setup_and_run` builds a context (`src/utilities.jl:863`).
`src/problem.jl:396` gates the frozen `{1,2,4,5}` contribution on `st !== nothing`, so through the
MMS driver the `nl_pressure_full` branch adds **nothing at all on a flat bed**, and only the `𝓐` ∇h
IBP half on a sloping one.

The studies stay valid — `mms_forcing` is selected by the same symbol and *does* compute `{1,2,4,5}`
exactly, so forcing and solver really do encode different operators — but **what the floor measures
is OMISSION, not LAG**. That invalidates two statements in `VERIFIED_SCOPE.md` §4: that the lag is
"negligible, now *measured* rather than argued" (it was never exercised), and that this floor
characterises a *production* `:full` run (where the context exists and the projections genuinely are
frozen-and-lagged).

**Decisive next step:** add an `nlp` kwarg to `run_mms_case` — or build the context internally when
`nl_pressure=:full` — re-run the two `:full` studies and compare. **Floor drops sharply ⇒ it was
omission; floor barely moves ⇒ the original wording was accidentally right.** Either way the
`:full`-is-not-MMS-verifiable conclusion probably survives; only its reason is in question.

---

## 3. 🔴 `test_mms_convergence` G7 — needs a specification decision, not a fix

Diagnosed in full in `TEST_SUITE.md` §6: the temporal window is contaminated at both ends, `η`'s
asymptotic window is one refinement wide, and the cause is structural (the `Q3/Q2` pairing puts `η`
in the lower-order space). **Deliberately not re-specified, because the fix is a runtime-cost
decision:** re-specify at ~3× runtime, gate the two fields differently, or refine to `nx ≥ 48`.

**Do not widen the ±0.3 tolerance.** Also fix **G10**, which already computes the per-field
contamination (3.4 % for `u` vs 0.5 % for `η`) but collapses both into one 10 % threshold.

---

## 4. 🟠 No MPI tier in `runtests.jl`

The distributed set — now **four** files, `test_mms_distributed_parity.jl` added 2026-08-21 — needs
`mpiexecjl` and is only ever run by hand. That cost **three stale reference constants**, all found
the first time they were re-run in a session.

**Fix:** an MPI tier that runs them when `mpiexecjl` resolves and **skips LOUDLY** when it does not.
A silent skip recreates the same blind spot. `MPI_TESTS` in `runtests.jl` already lists all four with
their rank counts, so the tier has its inventory; what is missing is execution and verdict capture
(from **gate output**, never exit codes).

---

## 5. 🟠 Preconditioner replacement — the single biggest performance item

GMRES needs ~300–770 iterations where a well-preconditioned solve of this size should need tens.
Both cheap drop-ins are measured dead (`CONFIGURATION.md` §5). Remaining, in order of effort:
**field-split/Schur** → **geometric multigrid** → **restricted Schwarz** (needs library support
GridapSolvers 0.7.1 does not expose).

**Cheap reproducible test bed:** the ring-vs-plane gap on an identical mesh (654–695 vs 451–517).
A preconditioner that closes it likely helps everywhere.

---

## 6. 🟠 Four run-output gaps

In priority order:

1. **The discrete-equation check does not run under the default integrator.** It is `is_theta`-gated,
   so `res_theta = NaN` under SDIRK. Highest value: it is the only *in-run* verification that the
   accepted state satisfies the discretised governing equations.
2. **No per-Newton-iteration linear counts.** `SolverMonitor` samples the `ConvergenceLog` once per
   `solve!`, i.e. per RK stage, so `gmres 451/517` is min/max over *stages*. A Newton step whose
   first linear solve is cheap and second expensive is invisible.
3. **No wall-time breakdown.** `t_solve` lumps Jacobian assembly and the linear solve together, so a
   55 s/step 2-D cost cannot be attributed from the log alone.
4. **The linear solve's achieved tolerance is not recorded** — only its iteration count. A solve that
   stopped on `atol` rather than `rtol` is indistinguishable from one that converged well.

---

## 7. Verification gaps

* **`Q2/Q1` and `Q4/Q3` velocity shortfall** — ✅ **answered for the vertical-basis case
  (2026-08-29): PRE-ASYMPTOTIC.** P1LFE-4's `pw_u` ran 3.33 → 3.45 → 3.78 over `nx = 8…64` and
  reaches **3.94 at `nx = 128`**, error still falling ~15× per refinement
  (`MMS_VBASIS_CAMPAIGN.md` §2.5). The method is the transferable part: **extend the ladder and read
  the pairwise SEQUENCE**, never a slope fitted over a window that never reached the asymptotic
  regime — the fit read 3.515/3.729 against a true 3.94.
  ⚠ Still open for the **horizontal FE pairings** (`Q2/Q1`, `Q4/Q3`), which were not re-run on an
  extended ladder. The same diagnostic applies and is now cheap to run.
* **Error-vs-DOF study at production resolution** — **must precede any change to the default FE
  pairing.** `Q3/Q3` was 40× more accurate than `Q3/Q2` at `nx=24` despite the worse rate.
* **`:full` dynamic sensitivity has never been measured.** The block's relative effect scales as
  `A²`: 0.013 % at `A=1e-3` but **~1.3 % at `A=1e-2`**. The 2-D bar case already ran at `A=0.01`
  with `:full`; **its `:none` counterpart was never run**, and that single pair would both resolve
  the difference and confirm the `A³` scaling dynamically.
* **Only one vertical resolution and one `kd`** exercised locally (P1LFE-2, `Nσ=3`, `kd=5.5`).
  The convergence half of this is designed AND now **unblocked** — see
  [`COMPLETED_VBASIS_STUDY.md`](COMPLETED_VBASIS_STUDY.md) §1: both `src/mms_driver.jl` defects were fixed 2026-08-21,
  every MMS driver/test/example takes `(M, p_vert)`, and the sweep driver is
  `examples/local_mms/run_vertical_basis_study.jl`. **Decisive next step: run tier 1.** Two spot
  checks already reach theoretical order off the default basis (P1LFE-3 and P2LFE-2, Q2/Q1 1-D
  static Model 1: `p_η` 1.994, `p_u` 2.997/2.996).
  A `T` sweep measuring sponge reflection at `kd = 1, 3, 5.5` would directly test the requirement
  that the sponge cover the longest component — that half is still untouched.
* A run-and-reconstruct **pressure** profile test; a **distributed-gauge** utility.

---

## 8. Runs not yet made → moved

The list of runs that are designed and unblocked but not yet executed now lives in
[`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md) §7, which is the single place forward work is
tracked. **This file is for gaps in work already done**; a run that has simply not been launched is
not a gap, it is a plan.

---

## 9. Naming follow-through (outside this checkout)

The `GridapLFEM → GridapBALFEM` rename is complete, committed and pushed. Remaining:

* the **GitHub repository** is still named `GridapLFEM.jl` — rename it, then
  `git remote set-url origin …GridapBALFEM.jl.git`;
* the **cluster checkout** must be renamed to match (`run/balfem_env.sh` defaults
  `BALFEM_PROJ=$HOME/GridapBALFEM.jl`);
* the **sysimage must be rebuilt** as `GridapBALFEM_sysimage.so`.

**Deliberately NOT renamed:** `../LFE-M_2D_solver/` and its `LFEModel2D` module (a real legacy
directory), and `CFC2027_LFEMultilayer_abstract.pdf` (an existing artefact).

---

## 10. Rules for anyone doing a repository-wide rename here

Three earned the hard way, and all three are cheap to respect:

1. **Use a single regex alternation, longest-first, with protected tokens mapping to themselves.**
   An *ordered cascade* of find/replace rules re-matches its own output
   (`GridapLFEM → GridapBALFEM → GridapBABALFEM`).
2. **Renaming file CONTENTS and file NAMES are two jobs.** 39 launchers ended up sourcing a helper
   under its new name while the file still had the old one. Verify by **resolving every sourced
   path**, not by grepping.
3. **A blanket rename corrupts semantics exactly where the distinction matters most** — it rewrote
   the block *defining* the three model names. Guard such passages with an `assert count == 1`.
