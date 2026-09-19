# NEW_TREATMENT.md — replacing the frozen Class-III projection

> **Branch `new-classIII-treatment`, opened 2026-09-18.**
> Two changes to how the Class-III components `{1,2,4,5}` of `𝓝ₖⱼ` reach the residual:
> **(A)** an exact algebraic reduction collapsing `{1,2,5}` onto a single contraction, and
> **(B)** moving the `L²` projection *inside* the nonlinear solve, turning it from a lagged
> approximation into static condensation of the mixed formulation.
> Neither changes the continuous model. Background and evidence: `OPEN_ISSUES.md` §0c.

---

## 0. What problem this addresses, in one paragraph

`nl_pressure=:full` diverges on a flat bed in 8 wave periods at `A = 0.10` m while `:native`
completes 100 periods. The operator is verified against Yang & Liu to 2.1e-12 (`CLAUDE.md` §5.2d) and
the Jacobian is exonerated (rule 17b), so what remains is the **assembly** of `{1,2,4,5}` — carried
by `L²` projections that are **frozen from the previous time step**. The refinement signature splits
into two errors: a **recovery** error with a `dx` signature (halving `dx` advances onset ×2.5) and a
**lag** error with a `dt` signature (halving `dt` delays it ×1.7, the stronger and cleaner trend).
**This branch removes the lag entirely and shrinks what the recovery has to carry.** It does not
attempt the `C⁰`-interior-penalty treatment of the recovery error; that remains open, and this work
is a prerequisite for judging whether it is needed.

---

# PART A — the algebraic reduction

## A.1 Which components, and why those

From `GlobalResidual.tex` §`sec: pressure operator implementation`, the eight components sort by what
they demand of the unknowns:

| class | components | demand |
|---|---|---|
| I | 6, 7, 8 | products of first-order fields only |
| II | 3 | one second derivative, landing on the **analytic** bed Hessian `∂²h` |
| **III** | **1, 2, 4, 5** | **second derivatives of the unknowns** |

⚠ **The demand is a SECOND derivative, not a third.** `s_k = ∇·(Hu_k)` is first order; `∇s_k` is
second order in `H·u`. Any statement in the project history saying `∂³` is wrong — corrected
2026-09-18 across `OPEN_ISSUES.md` §0c, `CLAUDE.md` §5.2c and `PLANNED_CAMPAIGNS.md` §6b.

Writing the four Class-III components out (`GlobalResidual.tex`, eq. *N1, N2, N4, N5 class III
expressions*):

```
𝓝¹ₖⱼ = −uⱼ·∇s_k
𝓝²ₖⱼ = ∇·(s_k uⱼ) = uⱼ·∇s_k + s_k(∇·uⱼ) = −𝓝¹ₖⱼ + s_k(∇·uⱼ)
𝓝⁴ₖⱼ = u_k·∇bⱼ
𝓝⁵ₖⱼ = −u_k·∇sⱼ
```

**Components 1, 2 and 5 are the three that are built from `∇s`.** Component 4 is built from `∇b`,
`b_j = uⱼ·∇H` — a *different* field, requiring its own projection. So `{1,2,5}` is exactly the right
set: it is the maximal subset sharing one inadmissible object.

## A.2 The reduction

Two facts, both exact:

1. **`𝓝²` is `−𝓝¹` plus an admissible remainder.** `s_k(∇·uⱼ)` is a product of first-order fields —
   Class I in everything but name.
2. **`𝓝⁵` is the `(k,j)` transpose of `𝓝¹`.** `𝓝⁵ₖⱼ = −u_k·∇sⱼ` is `𝓝¹` with the roles of the two
   layer indices exchanged.

Fact 2 only helps once the **tensor weights** are accounted for, and that is where it is easy to get
wrong. The residual contracts each component against its own constant third-order tensor:

```
contribution_i  =  Σ_{ℓ} Σ_{k,j}  T^{(ℓ)}[i,k,j] · 𝓝^{(ℓ)}[k,j]
```

with `T ∈ {𝓐₃, 𝓚₃, 𝓟₃}` depending on the residual block. Write `G_a := ∂_a π(s)` (a layer-vector
field, `a ∈ {x,y}`) so that

```
𝓝¹ = −Σ_a G_a ⊗ U_a            𝓝² = +Σ_a G_a ⊗ U_a + S ⊗ DU            𝓝⁵ = −Σ_a U_a ⊗ G_a
```

using the code's own convention `(a⊗b)[k,j] = a[k]·b[j]` (`src/tensors.jl:61`). Substituting and
relabelling the dummy pair `k ↔ j` in the `𝓝⁵` term:

```
Σ_{k,j} T⁵[i,k,j] U_a[k] G_a[j]  =  Σ_{k,j} T⁵[i,j,k] G_a[k] U_a[j]
```

which gives the result:

> ### The reduction
> ```
> Σ_{ℓ∈{1,2,5}} T^{(ℓ)} ⊙ 𝓝^{(ℓ)}   =   W ⊙ ( Σ_a G_a ⊗ U_a )   +   T⁽²⁾ ⊙ ( S ⊗ DU )
>
> with        W[i,k,j]  =  −T¹[i,k,j] + T²[i,k,j] − T⁵[i,j,k]
> ```
> ⚠ **`T⁵` enters TRANSPOSED IN ITS LAST TWO INDICES.** This is the whole content of fact 2 and the
> one place an implementation can silently go wrong — `𝓐₃`, `𝓚₃`, `𝓟₃` are **not** symmetric in
> `(k,j)`, so `T⁵[i,k,j]` and `T⁵[i,j,k]` are different tensors.

**`W` is constant**, built once from the vertical tensors at problem-construction time, exactly like
`𝓐₃`/`𝓚₃`/`𝓟₃` themselves.

## A.3 Verification — done before implementing

Checked against **Gridap's own** `double_contraction` and `outer`, not against an assumed index
convention (`Nσ = 3`, random tensors and fields):

| check | result |
|---|---|
| `double_contraction(T,S)[i] == Σ_{kj} T[i,k,j]S[k,j]` | **0.0** — convention confirmed |
| `outer(a,b)[k,j] == a[k]b[j]` | **0.0** — convention confirmed |
| **direct `{1,2,5}` sum vs the reduced form** | **4.4e-16** on a scale of 1.64 |
| **control: same, with `T⁵` NOT transposed** | **0.283** — must be large, and is |

The control is the important row: it proves the test can *detect* the error the reduction is most
likely to make. This is rule 34 (the resolution principle) applied before writing any code.

## A.4 What the reduction does and does not buy

**Does:**
* Reduces the Class-III contractions per frozen block from **four to three** (`W`, `T²` with the
  admissible `S⊗DU`, and `T⁴` with `𝓝⁴`).
* **Isolates exactly one inadmissible object per block**: `∇π(s)` with weight `W`. The `s_k(∇·uⱼ)`
  half of `𝓝²` is moved into manifestly admissible territory and never needs a projection at all.
* Makes the target of any future treatment unambiguous. If `C⁰`-IP is ever built, it has **one**
  term to act on, not three.

**Does not:**
* It is **not** a significant cost saving. The current code already forms `∇π(s)` once and reuses it;
  the saving is one tensor contraction per block per quadrature point.
* It does **not** remove the Class-III problem. One `∇s` contraction and one `∇b` contraction remain.
* ⚠ It changes **nothing** about the answer. The identity is exact, so the reduced residual must
  equal the current one to round-off — which is precisely the regression gate in §C.2.

## A.5 The bed-slope (`𝓐`) block is deliberately NOT touched

`nlp_gradh_contrib` handles `{1,2,4,5}` in the `∇h` block by **exact integration by parts** onto the
test function — no projection, no lag, no approximation. The same algebraic identity would apply
there, but:

* that path is already exact and gated by `test_nlpressure`;
* it **vanishes identically on a flat bed**, which is where every instability on record occurs, so
  changing it cannot affect the phenomenon under investigation;
* rewriting a working exact path to save one contraction is risk without reward.

**Decision: `nlp_gradh_contrib` is left byte-for-byte unchanged.** The reduction is applied only to
the frozen `𝓚` and `𝓟` blocks, i.e. only where a projection is actually used.

## A.6 An alternative that was checked and rejected

`s_k = H(∇·u_k) + b_k` identically (it is how `S` is built in the code). So
`∇b_j = ∇s_j − (∇H)(∇·u_j) − H∇(∇·u_j)`, which would let component 4 reuse `π(s)`.

**Rejected**: it trades the projection of `b` for a projection of `∇·u` — still two projections, and
it introduces `H∇(∇·u)`, a *new* second-derivative object with no better standing than the one it
replaces. No reduction in count, no reduction in inadmissibility.

⚠ A second variant **is** worth recording for later, though it is out of scope here: expanding
`∂_a b_j = (∂_a u_{j,b})(∂_b H) + u_{j,b}(∂_a∂_b h + ∂_a∂_b η)` shows that the *only* inadmissible
part of `𝓝⁴` is `∂²η`. Projecting `∇η` alone — instead of the whole of `b` — would leave the rest of
component 4 exact. That is a genuine accuracy improvement, deferred so this branch changes one thing
at a time.

---

# PART B — static condensation: the projection inside the Newton loop

## B.1 The projection and the mixed formulation are the same construction

`GlobalResidual.tex` presents two "sound options" for Class III as alternatives:

* **mixed formulation** — promote `G ≈ ∇s` to a genuine unknown with its own weak equation, solved
  monolithically. Rejected on system size.
* **`L²` projection** — project `s`, differentiate the projection, solve a small mass system.

⚠ **These are not different methods. They are the same auxiliary-field construction with different
elimination strategies.** The mixed route keeps `G` in the monolithic unknown vector; the projection
route eliminates it against a constant, SPD, once-factorised mass matrix. Eliminating an auxiliary
field whose equation is cheap and decoupled is **static condensation** — standard, exact, and it does
not enlarge the system.

**The lag is not part of either construction.** It was added purely to decouple the mass solve from
the Newton iterate:

> *"…defined as we just did, (projection) still couples `π(g)` to the current Newton iterate through
> `g = g(H,U)`. In order to decouple this projection from the global system solve, projections are
> evaluated from the previous accepted time step…"*

That decoupling is what introduces the `O(Δt)` consistency error the `dt` ladder measures. **Removing
it costs no new unknowns and no change of function space** — only a mass solve per residual
evaluation, on a matrix that is already assembled and factorised once at start-up.

## B.2 What "inside the Newton loop" means here

The residual is evaluated at every Newton iteration with the current iterate `u^(m)`. Refresh the
projection **there**, from `u^(m)`, rather than once per accepted step from `u^n`:

```
lagged  (now):   r( u^(m) ; π(s(u^{n-1})) ) = 0        ← π from the PREVIOUS STEP
in-loop (new):   r( u^(m) ; π(s(u^(m)))  ) = 0        ← π from the CURRENT ITERATE
```

At convergence the second is the **fully consistent** condensed system: `π` is evaluated at the same
state the residual is. The `O(Δt)` lag error is gone — not reduced, *gone*.

**The Jacobian stays quasi-Newton.** We do not differentiate through `M⁻¹`. By rule 17b that costs
iteration count and nothing else: the converged root is defined by `r(u) = 0`, and `r` is now the
consistent residual. This is the same quasi-Newton choice already made for these blocks, so nothing
about the Jacobian story changes.

## B.3 The one real hazard: the AD path

The default operator is `TransientFEOperator(r, j, jt, U, V)` — hand Jacobians, and **the residual is
never differentiated**. But `use_ad=true` builds `TransientFEOperator(r, U, V)`, and Gridap then
autodiffs `r`, calling it with **Dual-valued** cell data. A mass solve inside the residual would meet
`ForwardDiff.Dual` numbers where it expects `Float64`.

**Guard:** the refresh runs only when the incoming iterate carries plain `Float64` free values. On
the AD path it is skipped, and the Jacobian is taken with `π` held at the value the preceding
residual call computed. That is exactly the intended quasi-Newton treatment, and it keeps the
hand-vs-AD comparison of `test_jacobians_ad` apples-to-apples.

⚠ **This must be verified, not assumed** (rule 38d): the gate in §C.3 asserts that the refresh
actually fires on the hand path — a dead knob would give a clean, confident, entirely wrong negative.

## B.4 Cost

| item | cost |
|---|---|
| mass matrix assembly + factorisation | unchanged — once, at start-up |
| mass solves | 2 per **residual evaluation** instead of 2 per **step** — so ≈ `n_newton ×` more |
| new unknowns | **none** |
| new Jacobian blocks | **none** |
| conditioning of the Newton system | **unchanged** — the mass solve is outside it |
| stencil / sparsity | **unchanged** |

At the measured ~4–5 Newton iterations/step this is ~4–5× the projection work, against a mass solve
that is a back-substitution on a pre-factorised SPD matrix — small next to the LU solve of the full
system each iteration. **Measure it** (§C.4); do not assume it is negligible.

⚠ The post-step `update_nlp_state!` becomes redundant in this mode and is **skipped**, recovering two
of those solves.

---

# PART C — the testing programme

**All stability testing is at 1-D Q3/Q2**, for two reasons: it is the production pairing and the one
the entire MMS scope was measured on; and at Q2/Q1 the free-surface Hessian `∂²η ≡ 0` identically in
1-D, so half of what Class III is supposed to carry is structurally absent and a Q2/Q1 result cannot
speak for the treatment.

## C.1 `test/test_class3_reduction.jl` — the algebra gate (fast, `:fast` tier)

Pure tensor algebra on the **real** vertical tensors from `assemble_vertical`, not random ones.

| gate | assertion |
|---|---|
| G1 | Gridap's `double_contraction` contracts the trailing two indices — pins the convention the derivation assumes |
| G2 | reduced ≡ direct for `𝓚₃`, on P1LFE-2 · P1LFE-3 · P2LFE-2 — `< 1e-14` relative |
| G3 | same for `𝓟₃` |
| G4 | **control**: the same comparison with `T⁵` *not* transposed must FAIL by a margin ≥ 1e-3 |
| G5 | `W` is built once and is independent of the state (constant tensor) |

G4 is the gate that makes G2/G3 meaningful.

## C.2 `test/test_class3_residual_parity.jl` — the regression gate

The reduction is **exact**, so with an identical frozen state the new residual must reproduce the old
one to round-off. Assemble both forms on the same mesh, same state, same projections, and compare
the assembled residual **vectors**.

| gate | assertion |
|---|---|
| G1 | `‖r_reduced − r_direct‖_∞ / ‖r_direct‖_∞ < 1e-12`, flat bed, Q3/Q2 |
| G2 | same on a **sloping** bed — exercises the `𝓐` IBP path alongside, confirming it was not disturbed |
| G3 | same at `Nσ = 4` (P1LFE-3) — the identity must not depend on the layer count |

⚠ This is the gate that would catch a wrong `W`. It is a *self-consistency* test (rule 28) and that is
appropriate here: the claim being tested is an algebraic identity between two implementations, not a
physical claim.

## C.3 `test/test_nlp_inloop.jl` — the condensation gate

| gate | assertion |
|---|---|
| G1 | **the knob is LIVE** — with in-loop on, `π` changes between the first and last Newton iteration of a step by more than `1e-12` (rule 38d) |
| G2 | **consistency at convergence** — after the step converges, recomputing `π` from the accepted state changes it by `< nl_tol`; i.e. the fixed point is reached |
| G3 | **lagged and in-loop agree as `dt → 0`** — the difference in `η` falls at ≥ first order in `dt`, confirming the two differ only by the lag |
| G4 | in-loop mode leaves `:none` and `:native` **bit-identical** — the change must not touch tiers that carry no Class-III blocks |
| G5 | AD path still runs (the guard does not throw) and hand-vs-AD Jacobians still agree |

## C.4 The stability factorial — the point of the exercise

The case that dies: `P1LFE-2`, 1-D flume, `A = 0.10` m, `kd = 5.5`, `dx = 0.25`, `dt = 0.04`,
SDIRK_2_2, `nl_pressure=:full`, flat bed. **Known onset at Q2/Q1 lagged: t = 12.60 s.**

| # | pairing | projection | reduction | expectation |
|---|---|---|---|---|
| 0 | Q2/Q1 | lagged | no | **control — must reproduce t = 12.60 s** |
| 1 | Q3/Q2 | lagged | no | isolates the pairing |
| 2 | Q3/Q2 | **in-loop** | no | isolates the lag — **the measurement that matters** |
| 3 | Q3/Q2 | in-loop | **yes** | must match #2 to round-off (the reduction is exact) |

⚠ **Run #0 first and confirm 12.60 s.** A factorial whose control does not reproduce the known number
is measuring something else (rule 38c, and rule 12b's root-cause-of-the-root-cause).
⚠ **Budget ≥ 100 wave periods** for any run that does *not* crash. The flume fills in 36 s (rule 14),
the 1:2 bar looked stable at 80 s and died at 137 s, and a run that ends early proves nothing
(rule 12c). Nothing before t ≈ 36 s may be read as growth.

**Reading the outcome:**
* #2 pushes onset out substantially or removes it ⇒ **the lag was the dominant mechanism**, and the
  `C⁰`-IP work may not be needed at all.
* #2 ≈ #1 ⇒ the lag was not the mechanism; the **recovery** error is, and `C⁰`-IP (or a `C¹` space)
  becomes the next step — now with one clean target thanks to Part A.
* #1 ≫ #0 ⇒ the pairing mattered and every earlier Q2/Q1 conclusion needs re-reading.

---

# PART D — the implementation plan

Ordered so that every step is verifiable before the next one depends on it.

### D.1 Tensors — build `W`
`src/tensors.jl`: add `alg_class3_weight(T1, T2, T5)` returning
`W[i,k,j] = −T1[i,k,j] + T2[i,k,j] − T5[i,j,k]`. Pure, constant, no CellFields.

### D.2 Problem — carry `W` and the projection context
`src/problem.jl`:
* `BALFEMProblem` gains `WK3`/`WP3` (the reduced weights for the `𝓚` and `𝓟` families) and
  `nlp_ctx :: Base.RefValue{Any}`. **Presence of a ctx means in-loop mode**; `nothing` keeps the
  legacy lagged behaviour, so the default is unchanged.
* `build_problem_raw` computes `WK3`, `WP3` from `vert.Kcal`, `vert.Pcal` and passes
  `Ref{Any}(nothing)`.
* ⚠ Adding struct fields ⇒ **rule 27**: restart the session, do not trust Revise.

### D.3 Nonlinear pressure — reduced blocks + in-loop refresh
`src/nlpressure.jl`:
* `nlp_class3_fields(Ux,Uy,S,DU,piS,pib)` → `(GU, SD, N4)` where `GU = Σ_a ∂_a(π𝖲) ⊗ U_a`.
* `nlp_gradH_reduced_contrib` / `nlp_P_reduced_contrib` using `W ⊙ GU + T² ⊙ SD + T⁴ ⊙ N4`.
* Keep `nlp_frozen_N` and the two existing contributors — the parity gate compares against them, so
  they are the reference implementation, not dead code.
* `refresh_nlp_state!(prob, u)` — the guarded in-loop projection.

### D.4 Residual — the hook
`src/problem.jl` `global_residual`: before assembling the Class-III blocks, if
`prob.nlp_ctx[] !== nothing`, call `refresh_nlp_state!`. Then use the reduced contributors.

### D.5 Time loops — skip the now-redundant post-step update
`src/timeloop.jl` and `src/timeloop_dist.jl`: keep the priming call (it seeds the first residual),
skip the post-step `update_nlp_state!` when in-loop mode is active.

### D.6 Drivers — expose the switch
`setup_and_run` (+ distributed, + MMS driver): `nlp_inloop::Bool = false`, wired to
`BALFEM_NLP_INLOOP` in `examples/local_1d/run_flume_1d.jl`. **Default off** until the factorial says
otherwise, so no existing result changes silently (rule 2c's lesson).

### D.7 Tests, then runs
`test_class3_reduction.jl` → `test_class3_residual_parity.jl` → `test_nlp_inloop.jl`, then the
four-run factorial of §C.4.

---

## Risks, and what would falsify this work

| risk | detection |
|---|---|
| `W` built with the wrong transpose | §C.1 G4 control + §C.2 parity |
| the refresh never fires (dead knob) | §C.3 G1 — the rule-38d check |
| the refresh throws on the AD path | §C.3 G5 |
| in-loop changes `:none`/`:native` | §C.3 G4 |
| the outer coupling fails to converge | §C.3 G2; symptom would be Newton counts rising without bound |
| **it simply does not fix the instability** | §C.4 #2. **This is a real possible outcome** — the projection is the leading hypothesis, not a proven cause, and the `dx` signature would survive this branch untouched |

---

# PART E — RESULTS AS MEASURED (2026-09-18)

## E.1 Implementation — complete

| step | file | state |
|---|---|---|
| D.1 `alg_class3_weight` | `src/tensors.jl` | ✅ |
| D.2 `WK3`/`WP3` + `nlp_ctx` on the problem | `src/problem.jl` | ✅ |
| D.3 reduced contributors + `refresh_nlp_state!` + AD guard | `src/nlpressure.jl` | ✅ |
| D.4 residual hook | `src/problem.jl` `global_residual` | ✅ |
| D.5 post-step update skipped in in-loop mode | `src/timeloop{,_dist}.jl` | ✅ |
| D.6 `nlp_inloop` kwarg + `BALFEM_NLP_INLOOP` | `src/utilities.jl`, `examples/local_1d/run_flume_1d.jl` | ✅ |
| D.7 three gates registered | `test/runtests.jl` | ✅ |

## E.2 Gates

| gate file | result |
|---|---|
| `test_class3_reduction.jl` | ✅ **19/19** |
| `test_class3_residual_parity.jl` | ✅ **4/4** |
| `test_nlp_inloop.jl` | ⚠ **G1, G1b, G2a, G2b, G4(×2) pass; G3 FAILS; G5 not reached** |

**The reduction is exact, measured, and controlled.** Reduced ≡ direct as assembled residual
vectors at **3.2e-16** (flat bed), **3.9e-16** (sloping), **5.6e-16** (`Nσ=4`), while the
no-transpose control lands at **0.283** and a 1 % mutation of `W` is detected at **1.6e-3**. The
`𝓐` IBP path is undisturbed (G2 runs with it live).

**The in-loop switch is live and confined.** `:full` moves by 4.4e-09 (rel 8.0e-05) — not a dead
knob (rule 38d) — while `:none` and `:native` are **bit-identical** to the last digit. The banner
prints which mode is active.

⚠ **Measured cost of in-loop mode: `4.85 s/step` against `4.85 s/step` lagged on the gate's mesh**
— i.e. not resolvable there, because the mass solve is a back-substitution on a pre-factorised SPD
matrix. This must be re-measured on the production mesh before being quoted (rule 38g).

## E.3 ⛔ G3 FAILED TWICE — AND THE SECOND FAILURE IS AN IMPLEMENTATION DEFECT

### First measurement (bad window)
```
dt = 0.040 → 4.189748e-09 ;  dt = 0.020 → 4.406931e-09 ;  ratio 0.95
```
Diagnosed as a defect in the **test**: the window ended at t = 0.24 s inside a 3.2 s ramp, at
`η ≈ 5.5e-05`, where the gap was the same order as the `O(A²)` terms themselves.

### Window fixed, re-measured
`T_ramp = 0.8 s`, `t_end = 2.4 s` (3× ramp, established wave), `A = 0.02 m`.

⚠ The first attempt at the fix used `A = 0.10` m and **Newton did not converge** — it stalled at
`‖r‖ = 9.250e-08` after 60 iterations. **That is the quasi-Newton cliff (rule 5) measured at
production amplitude for the first time**: with `{1,2,4,5}` absent from `jacobian_u`, ~1e-7 is the
floor on the attainable residual whatever tolerance is asked. `CLAUDE.md` had only the MMS figure
(9.2e-04 at `a_eta = 0.8`). The gate therefore runs at `A = 0.02`, where the hand Jacobian reaches
1e-10.

```
dt = 0.040   η = 2.05528672e-02   gap = 1.612753e-03   (rel 7.85e-02)
dt = 0.020   η = 2.06435371e-02   gap = 1.676368e-03   (rel 8.12e-02)
ratio 0.962        resolvable: gap is 1e7 × the Newton floor
```

**The window fix worked — the measurement is now unambiguous — and the gate STILL FAILS.** The gap
is flat in `dt` (it grows slightly). So the hypothesis G3 exists to test is **refuted, not
unmeasured**.

### ⛔ The discriminator: the in-loop path nearly CANCELS the Class-III blocks

`:native` omits `{1,2,4,5}` entirely, so `|η(:full lagged) − η(:native)|` is the size of the whole
Class-III contribution. Measured at the same settings:

| quantity | value |
|---|---|
| `η(:native)` | 1.89112821e-02 |
| `η(:full, lagged)` | 2.05528672e-02 |
| **Class-III contribution** | **1.641585e-03** |
| **lagged → in-loop gap** | **1.612753e-03** |
| **ratio** | **0.982** |

`η(:full, in-loop) ≈ η(:native)` to better than 2 %. **In-loop mode is not "the same operator without
the lag" — it is very nearly switching the Class-III blocks OFF.**

⚠ **THIS IS A DEFECT IN THE BRANCH'S IN-LOOP PATH, NOT A SUBTLETY ABOUT CONSISTENCY ORDER**, and the
project has seen this exact failure mode before: `OPEN_ISSUES.md` §2, where `run_mms_case` never
built an `nlp` context and `:full` silently degenerated to `:native` — blocks **absent**, not lagged.
The signature is identical.

⚠ **CONSEQUENCE FOR §E.4: run #2's survival past the control's 12.60 s onset is NOT evidence that
de-lagging cures the instability.** It is equally consistent with `:full` degenerating to `:native`,
which is already known to survive 100 periods. **No stability conclusion may be drawn from any
in-loop run until this is fixed.**

**Not yet diagnosed.** The refresh demonstrably fires (G1: the answer changes) and is confined (G4:
`:none`/`:native` bit-identical), so the fault is in *what* it computes or *when*, not in whether it
runs. Candidates not yet separated: the projections being read before they are written within an
assembly; the stage state at `t_{n+θ}` versus the accepted state; lazy evaluation of the
`DomainContribution` capturing a different `nlp_state` than the one the refresh stored.

### What still stands
The **algebraic reduction is unaffected and verified**: 19/19 + 4/4 at 3e-16, it is active in every
tier including the control, and the control reproduced the known onset exactly (§E.4). Part A is
sound; Part B is not.

## E.4 The stability factorial — LAUNCHED, NOT FINISHED

| # | pairing | projection | status at write-up |
|---|---|---|---|
| 0 | Q2/Q1 | lagged | ✅ **died at t = 12.600000, `u_max` = 2.5046, Newton at the 51 cap** — reproduces the documented onset to five digits (rule 17b records 2.5046). **The branch is bit-faithful to the pre-branch solver with in-loop off.** |
| 1 | Q3/Q2 | lagged | not started |
| 2 | Q3/Q2 | in-loop | ran past t = 13.4 s with `u_max` = 0.451 and Newton 5 — ⛔ **UNINTERPRETABLE**, see §E.3: in-loop nearly disables the blocks, so this may simply be `:native` |
| 3 | Q2/Q1 | in-loop | not started |

⚠ **NO STABILITY CONCLUSION MAY BE DRAWN YET.** Run #0 has not reached its onset, so the control is
unconfirmed; and a 100-period arm is a **multi-hour** proposition on this workstation (the earlier
Q2/Q1 100-period run cost 13 h 22 m, and Q3/Q2 with in-loop projections is measurably slower per
step). Whether removing the lag fixes the `:full` instability **is not yet known**, and §C.4's
reading rules apply only once #0 reproduces 12.60 s.

## E.5 Also corrected on this branch

The project record said Class III needs a **third** derivative. It needs a **second** (`s = ∇·(Hu)`
is first order, `∇s` second) — `GlobalResidual.tex` had it right. Corrected in `OPEN_ISSUES.md` §0c,
`CLAUDE.md` §5.2c/§5.7/rule 1b and `PLANNED_CAMPAIGNS.md` §6b, together with the Q2/Q1 degree counts
that followed from it (`∂²(Hu)` is piecewise **linear** there, not a piecewise constant).
