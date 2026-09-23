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
second order in `H·u`, so the **weak form's trial requirement** is a Hessian.
⚠ This does NOT contradict the literature's third-derivative claim: that is about the strong-form
momentum equation, which carries `∇p_nh` and therefore `∂³u`. Measured both ways 2026-09-23 —
see `OPEN_ISSUES.md` §0c. The reduction to `∂²` is bought by `R_P`'s integration by parts.

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

The project record said Class III needs a **third** derivative. In the WEAK FORM it needs a **second** (`s = ∇·(Hu)`
is first order, `∇s` second) — `GlobalResidual.tex` had it right. Corrected in `OPEN_ISSUES.md` §0c,
`CLAUDE.md` §5.2c/§5.7/rule 1b and `PLANNED_CAMPAIGNS.md` §6b, together with the Q2/Q1 degree counts
that followed from it (`∂²(Hu)` is piecewise **linear** there, not a piecewise constant).

---

# PART F — THE MIXED (projection-free) FORMULATION
*Branch `mixed-formulation-solver`, 2026-09-19/20.*

## F.0 ⛔ READ THIS FIRST: the harness these were run in is UNFIT, and most of Part F's
## stability numbers are VOID

> ✅ **SUPERSEDED BY PART H (2026-09-23) — the proper campaign has now been run.** Part F's step-1
> demand (re-run on the real 60 m flume with in-batch controls) was carried out as a 2 × 4 × 2
> factorial. **Read Part H for every stability number.** Part F is kept for the static results in
> F.1, the layout rationale in F.2, the still-open gaps in F.3, and the Class-III isolation in F.4,
> which was always trustworthy because it ran in the production launcher against a reproduced
> control.

The mixed formulation was exercised in a small purpose-built harness (6 m domain, 24 cells,
interior Gaussian wavemaker at x = 1.0, sponge from x = 4.0, `A` = 0.02). **A control was
not run in the same batch until the very end, and when it finally was, it failed too:**

| arm (identical settings) | n = 10 | n = 20 |
|---|---|---|
| **`:native` / projected** — the production tier that runs **100 periods** on the real flume | OK, η = 1.588e-01 | ⛔ **DIED** |
| `:full` (gs) / projected | OK, η = 1.580e-01 | ⛔ DIED |
| `:full` (gs) / **mixed** | ⛔ DIED | — |

⚠ **`:native` DYING IS THE FINDING.** It carries no Class-III terms at all, so this harness
destabilises for reasons that have nothing to do with anything Part F is about. The tell was
visible earlier and I did not act on it: `η` reaches **0.1588 from a requested `A` = 0.02** —
eightfold — and goes 0.0269 → 0.1588 in five steps. That is not a wave, it is a blow-up
already under way. Likely cause: `make_wavemaker_line` defaults to `sigma_wm = 1.5` m, so in
a **6 m** box the Gaussian source spans the domain and reaches into the sponge (rule 14b: the
interior source does not deliver `A_wave`, and the factor is geometry-dependent).

**Therefore: the mixed formulation has NOT been shown to be unstable.** Every "mixed dies"
statement in this session — the 7-field NaN, the 5-field NaN, the identical 38-iteration
failures — is a statement about the harness. Rule 14c, in the form that hurts: *a crash is
not evidence until its control runs*, and here the control was the last thing tried.

## F.1 What IS solid

| result | evidence |
|---|---|
| ✅ **𝖦 recovers ∇𝖲 essentially exactly** | vs an analytic `∂ₓ𝖲`: **1.2e-11** interior, **6.1e-12** at the boundary (test_mixed_formulation.jl G4) |
| ✅ **The boundary term is load-bearing** | dropping `∮𝖲Ψn` leaves the boundary cells wrong by a **relative 22.4**; assembling it fixes them to 6e-12. The `:full` mode is velocity-led AT THE INFLOW, so a diagnostic wrong by 22× there could not have discriminated anything |
| ✅ **The index-1 DAE solves** | zero rows in `∂R/∂u̇` are not a problem; `:none` with the auxiliary rows present is **bit-identical** to `:none` without them |
| ✅ **The Class-III blocks are live in the mixed path** | mixed `:full` − mixed `:native` = 6.05e-06 at 5 steps, the same order as the projected treatment's contribution — the G2 criterion Part B failed |
| ✅ **The two treatments differ by ~40 % of the terms' own size** | mixed `:full` − projected `:full` = 2.42e-06 against a 6.05e-06 contribution. Consistent with one being a poor approximation of `∇𝖲` — which is the hypothesis — and *not* the ~100 % signature of a deleted block |

⚠ All five are 5-step, early-time, tiny-amplitude results taken before the harness defect was
known. The first two are static and unaffected; the last three should be re-taken.

## F.2 The `𝖦`-only (5-field) layout

`c3_mask` now sizes the mixed layout:

```
c3_mask = (true, false)  →  5 fields  [η, 𝖴x, 𝖴y, 𝖦x, 𝖦y]      ← the diagnostic
c3_mask = (true, true)   →  7 fields  (adds 𝖥x, 𝖥y ≈ ∇𝖻)
```

guarded so a mask/layout mismatch errors instead of silently mis-indexing `u[4]`.

**Justified by measurement, not preference** (§F.4): component 4 is 1140× smaller than the
`∇𝖲` family and runs clean on its own. ⚠ **Consequence for comparisons: the 5-field path is
the projection-free version of the `gs` ARM**, so its like-for-like control is the *projected
gs* run, not plain `:full`.

⚠ **Dropping `𝖥` did NOT change the failure**: 7-field and 5-field both NaN'd at the
identical iteration count (38). Under F.0 that is explained — both were failing on the
harness — but it also means the cost saving was real and the conditioning argument untested.

## F.3 Known gaps in the mixed path

1. **THE DIAGNOSTICS DO NOT SUPPORT IT.** `build_run_diagnostics` interpolates onto a
   **3-field** space, so `x_at_max`, `diagnostics.csv` and the divergence guard are all
   unavailable for mixed runs. ⚠ This is why no mixed failure has a *location*, and location
   is the single most diagnostic quantity here — the projected failure pins at the inflow.
   **Fix this before the next mixed run.**
2. **Inf-sup untested.** `𝖦` sits in the velocity space (Q3). One order lower is the standard
   first alternative and has not been tried.
3. **Cost.** 7-field AD was ~180 s/step at `nx` = 24; 5-field is cheaper but unquantified.
   Hand Jacobians for the auxiliary rows are the obvious remedy if this path continues.
4. **Non-rest starts are inconsistent.** `make_initial_conditions_mixed` zeroes `𝖦`; from a
   non-rest state the constraint would have to be solved first.

## F.4 ✅ THE ONE UNAMBIGUOUS RESULT OF THIS SESSION: Class III is a SINGLE term

Run on the **real** 60 m flume — the configuration with a known, reproduced control
(`c3_0_Q2Q1_lagged` died at **t = 12.600000**, `u_max` = 2.5046, Newton 51):

| arm | outcome |
|---|---|
| **`gs`** — ∇𝖲 family only ({1,2,5} collapsed) | ⛔ **reproduces the failure**: dying at t = 12.00, `u_max` 0.61 → **3.17** in one output interval, `x_at_max` pinned at **0.25 m (the inflow)**, Newton → cap |
| **`gb`** — component 4 only (∇𝖻) | ✅ **ran the full 25 periods to t = 40.0**, η = 0.1065, `u_max` = 0.435, **Newton flat at 4** |

Statically, on the assembled residual (test_class3_split.jl): `|∇𝖲 arm|` = 1.2837e-02 against
`|∇𝖻 arm|` = 1.1271e-05 — **a ratio of 1140** — with the partition exact to 8.65e-15.

**So `∂²η` is exonerated twice over: by magnitude and dynamically.** The target is one object,
`W ⊙ (Σₐ ∂ₐπ(𝖲) ⊗ Uₐ)`.

⚠ `gb` survived 25 periods (40 s) and the flume fills in ~36 s, so this is "did not reproduce
the failure", **not** a stability claim (rules 12c/14). 100 periods would be needed for that.

⚠ **THIS RESULT IS TRUSTWORTHY WHERE PART F's ARE NOT** because it ran in the production
launcher configuration against a control that reproduced a known number to five digits.

⚠ **And it nearly did not survive its own first attempt**: the first pair of isolation
launchers put `export BALFEM_C3_MASK` **after** `balfem_local_run`, where it is a silent
no-op (rule 38h). Both arms then reproduced the CONTROL, bit-identical to nine significant
figures — which is the tell. The knob was moved before the run line and liveness re-checked
(the arms now differ at t = 1.8 s) before any of the above was read.

## F.5 Next steps, ordered

> ✅ **ALL THREE ARE NOW DONE OR RE-SCOPED — see Part H.** (1) the real-flume re-run is Part H;
> (2) the diagnostics were made mixed-aware before launch (H.0); (3) `𝖦` at one order lower is
> **still not run** and has been promoted — H4 reading (3) makes it the probe of the one mechanism
> that could explain the surviving `dx` signature. It is now step 2 of H.4.

1. **Re-run the mixed comparison on the REAL flume configuration**, not the small harness —
   the 60 m Q2/Q1 case with the reproduced 12.60 s control. Nothing else about the mixed
   formulation can be believed until this exists.
2. **Make the diagnostics mixed-aware** (F.3.1) first, or the run yields a time but no place.
3. Then, and only then: `𝖦` at one order lower, to close the inf-sup ambiguity.

---

# PART G — COMPLETING THE MIXED JACOBIAN

*Plan written 2026-09-22, before implementation.*

## G.1 Why

The mixed operator is solved with a **block-diagonal** Jacobian (`mixed_aux_jacobian`): the hand
Jacobian on the physics rows, the mass matrix on the auxiliary rows, and **both coupling blocks
omitted**. Writing the true stage Jacobian in block form, with `u = [η,𝖴x,𝖴y]` and `𝖦 = [𝖦x,𝖦y]`:

```
        ┌ A   B ┐            A = ∂R_phys/∂u   (hand Jacobian — present)
   J =  │       │            M = ∂R_𝖦/∂𝖦     = the auxiliary mass matrix (present)
        └ C   M ┘            B = ∂R_phys/∂𝖦   ⛔ OMITTED
                             C = ∂R_𝖦/∂u      ⛔ OMITTED
```

Newton with `[A 0; 0 M]` is a block Jacobi iteration on the coupling, and its convergence rate is
governed by `‖A⁻¹B M⁻¹C‖`. **Measured cost of that omission:**

| pairing | Newton iterations/step |
|---|---|
| projected (3-field, hand Jacobian) | 6–8 |
| **mixed Q2/Q1** | **26–30** |
| **mixed Q3/Q2** | **64–88, oscillating** — and the run stalled at t = 8.0 s |

⚠ **At Q3/Q2 this is not merely slow, it is disqualifying.** A step that terminates at the
iteration cap is **not converged**, and rule 17b's licence to ignore an inexact Jacobian applies
only to *converged* steps. A capped mixed run therefore has a Jacobian-dependent trajectory and
cannot be compared against anything. Completing `C` (and `B`) is what makes the Q3/Q2 arm
interpretable at all.

## G.2 The derivation — both blocks are first-order and admissible

The auxiliary row for `𝖦ₐ` (test function `Ψ`) is

```
R_𝖦ₐ = ∫ 𝖦ₐ·Ψ  +  ∫ 𝖲·(∂ₐΨ)  −  ∮ 𝖲 Ψ nₐ ,        𝖲 = ∇·(H𝗎)
```

so `C` is entirely the `𝖲`-dependence. Differentiating `𝖲` in the direction `(dη, d𝗎)` and using
`H = h + η` (so `δH = dη`, `δ∇H = ∇dη`, the bed being data):

> ```
> δ𝖲 = ∇·(dη·𝗎) + ∇·(H·d𝗎)
>     = dη(∇·𝗎) + ∇dη·𝗎  +  H(∇·d𝗎) + ∇H·d𝗎
> ```

⚠ **Every term is a product of FIRST derivatives.** `C` introduces no new admissibility problem —
which is the point: `𝖲` was always first-order, and it is `∇𝖲` that was inadmissible, which is
exactly why `𝖦` exists. In the code's stacked notation:

```
δ𝖲 = dη*DU + (∂ₓdη)*𝖴x + (∂ᵧdη)*𝖴y  +  H*(∂ₓd𝖴x + ∂ᵧd𝖴y) + (∂ₓH)*d𝖴x + (∂ᵧH)*d𝖴y
```

and then

```
C·(dη,d𝗎) :   ∫ δ𝖲·(∂ₓΨgx) − ∮ δ𝖲 Ψgx nₓ      (the 𝖦x row)
              ∫ δ𝖲·(∂ᵧΨgy) − ∮ δ𝖲 Ψgy nᵧ      (the 𝖦y row)
```

`B` is simpler still. The Class-III blocks enter the physics rows through
`GU = Σₐ 𝖦ₐ ⊗ 𝖴ₐ`, which is **linear in 𝖦**, so

```
δGU = d𝖦x ⊗ 𝖴x + d𝖦y ⊗ 𝖴y
```

and `B` is that object pushed through the same two reduced contributors the residual uses
(`nlp_gradH_reduced_contrib`, `nlp_P_reduced_contrib`) with `SD` and `N4` held fixed — since
neither depends on `𝖦`.

⚠ **With `𝖥` active** (`c3_mask[2]`, the 7-field layout) the same applies with
`δ𝖻 = ∇dη·𝗎 + ∇H·d𝗎` and `δN4 = Σₐ 𝖴ₐ ⊗ d𝖥ₐ`.

## G.3 What is still omitted afterwards, deliberately

`∂(\text{Class-III})/∂(η,𝗎)` — the dependence of `GU`, `SD`, `N4` on the *physics* unknowns through
`𝖴ₐ`, `H` and `DU`. ⚠ **This is not a new approximation:** it is precisely the omission
`jacobian_u` already makes on the projected path, where the `{1,2,4,5}` blocks are absent from the
Jacobian entirely. Keeping it keeps the two paths' Jacobians comparable, which matters because the
whole experiment is a projected-vs-mixed comparison. After G.2 the mixed Jacobian is **no less
complete than the projected one**, which is the right target — not exactness for its own sake.

## G.4 Implementation steps

| # | change | file |
|---|---|---|
| G4.1 | `mixed_aux_jacobian` → `mixed_coupling_jacobian(prob, u, du, v, trian, dΩh, dΓ, nΓ)`: mass block **+ C + B** | `src/mixed.jl` |
| G4.2 | `build_ode_operator_mixed`: pass the state `u` and the boundary measure into the Jacobian closure (currently it gets neither) | `src/mixed.jl` |
| G4.3 | keep the old block-diagonal form behind `coupling=false`, so the two can be compared | `src/mixed.jl` |
| G4.4 | expose `BALFEM_MIXED_COUPLING` (default **on**) | `run_flume_1d.jl` |

⚠ **G4.2 is the structural change.** The present Jacobian closure is `jac(t,u,du,v)` and never sees
`dΓ`/`nΓ`; `C` needs the boundary term, because the residual has one and an inconsistent Jacobian
would reintroduce exactly the boundary defect §F measured at a relative 22.4.

## G.5 How it will be tested, before any run

1. **AD oracle** (the decisive one). Assemble the mixed Jacobian both ways — hand-coded with
   coupling, and by AD of `global_residual_mixed` — and compare **block by block** on a prescribed
   state, as `test_jacobians_ad.jl` does for the 3-field path. ⚠ Compare the `C` and `B` blocks
   *individually*, not just the total: a total agreeing while two blocks are wrong with cancelling
   errors is exactly what a summed comparison cannot see (rule 36).
2. **Amplitude scaling.** The remaining gap (G.3) must vanish with amplitude at the documented
   order 1.11–1.16; a gap that does not scale means something is wrong that is *not* the deliberate
   omission (rule 5).
3. **Newton count, the acceptance criterion.** On the Q2/Q1 baseline the count must fall from
   26–30 toward the projected path's 6–8. ⚠ **If it does not fall, the implementation is wrong** —
   that is the whole purpose of the change, and no stability result should be read from it until it
   does.
4. **Residual unchanged.** The converged solution must match the block-diagonal one to solver
   tolerance: the Jacobian sets the path, not the root (rule 17b). This is the regression that
   catches a `C` or `B` that is subtly wrong rather than merely inefficient.

## G.6 Then, and only then

Re-launch the Q3/Q2 arms of the §G.1 configurations — the same four set-ups the Q2/Q1 batch uses —
at 100 periods.

---

# PART H — THE MIXED CAMPAIGN ON THE REAL FLUME (2026-09-22/23)

*This is the run F.5 step 1 demanded: the mixed formulation on the production 60 m configuration,
with a projected control for every arm in the same batch. It supersedes every stability number in
Part F, which F.0 had already voided for harness defects.*

## H.0 Design

15 arms, `P1LFE-2`, `nl_pressure=:full`, **flat bed**, SDIRK_2_2, boundary-generated regular wave,
`kd` = 5.5, 60 m flume, `ny` = 1 + walls, 100 periods requested, diagnostics every 0.2 s.
A 2 × 4 factorial, each cell run **twice** — once projected, once mixed:

* **pairing**: Q2/Q1 (`c3v_*`) and Q3/Q2 (`c3q_*`);
* **cell**: `base` (nx=240, dt=0.04, A=0.10) · `amp` (A=0.15) · `dx` (nx=480) · `dt` (dt=0.02).

Launchers `run/local/run_1dc3{v,q}_{base,amp,dx,dt}_{proj,mixed}.sh`; output `output/local_1d/c3*`.

**Three prerequisites were discharged before launch, and each had previously invalidated a batch:**
1. **diagnostics made mixed-aware** (F.3.1) — `build_run_diagnostics` now uses `_n_multifields(U)`
   and `fill(zf, nf-3)` instead of a hard-coded 3-field list. Without this a mixed failure has a
   time but **no place**, and location is the single most diagnostic quantity here;
2. **the coupled Jacobian** `C = ∂R_𝖦/∂(η,u)` and `B = ∂R_phys/∂𝖦` (Part G), gated at 1.99e-11 and
   2.08e-09 against a finite-difference oracle, and measured to restore Newton to 3.10 it/step
   against the block-diagonal 7.17 and the projected 3.13;
3. **`BALFEM_MIXED` verified LIVE** — the two arms differ from t = 0.4 s and print different DOF
   counts (14406 for 3-field vs 31710 for 5-field at Q3/Q2). Rule 38h cost this branch a full
   isolation batch once already.

## H.1 Results — every arm

Onset = last diagnostics row before Newton fails to converge (50-iteration cap, or NaN).

| cell | pairing | free DOFs | **mixed** | **projected** | gain |
|---|---|---|---|---|---|
| `base` | Q2/Q1 | 14901 | ✅ **running > 97.4 s** | 12.60 (documented control) | **> 7.7×** |
| `amp` A=0.15 | Q2/Q1 | 14901 | ✅ running > 93.2 s ⚠ watch | 5.8 | > 16× |
| `dt` = 0.02 | Q2/Q1 | 14901 | ✅ running > 57.3 s | 21.2 | > 2.7× |
| `dx` nx=480 | Q2/Q1 | 29781 | ⛔ **36.6** | 5.0 | 7.3× |
| `base` | Q3/Q2 | 31710 | ⛔ **16.0** | 7.0 | 2.3× |
| `amp` A=0.15 | Q3/Q2 | 31710 | ⛔ 9.6 | 4.6 | 2.1× |
| `dt` = 0.02 | Q3/Q2 | 31710 | running 19.3 | 10.4 | > 1.9× |
| `dx` nx=480 | Q3/Q2 | 63390 | ⛔ 6.4 (NaN) | 3.4 | 1.9× |

**The Q2/Q1 `base` trace, windowed maxima over 10 s bins** — this is the headline object:

```
 t= 0- 10  eta 0.1129  u 0.4147      t= 50- 60  eta 0.1112  u 0.4043
 t=10- 20  eta 0.1168  u 0.4358      t= 60- 70  eta 0.1112  u 0.4076
 t=20- 30  eta 0.1168  u 0.4341      t= 70- 80  eta 0.1113  u 0.4097
 t=30- 40  eta 0.1157  u 0.4277      t= 80- 90  eta 0.1114  u 0.4112
 t=40- 50  eta 0.1112  u 0.4023      t= 90-100  eta 0.1115  u 0.4122
```

Flat to the fourth decimal from t = 40 on, Newton 5–6, i.e. **indistinguishable from the `:native`
reference trace** (η 0.10494 → 0.10300, Newton 5.12) at the tier that previously died in 12.6 s.

⚠ **All eleven failures share one signature**: `x_at_max` snapping back and **pinning at the inflow**
(x = 0.50, 0.25, 0.12, 0.56) while `η` stays bounded and `u` runs away — the velocity-led mode of
`CLAUDE.md` §5.2b/§5.6. **The mixed formulation changed *when*, never *what*.**

⚠ **Diagnosing the deaths correctly required reading past the traceback.** All four logs end inside
`paraview_collection` / `createpvd` (WriteVTK), which reads as a VTK failure; it is the enclosing
`do`-block. The real line is higher up: `ERROR: LoadError: Newton solver did not converge after 50
iterations` (or `NaN`). **A stack whose top frame is a writer does not mean the writer failed.**

## H.1b ⚙ FINAL CAPTURE 2026-09-23 15:29 — AMPLITUDE SPLITS THE Q2/Q1 RESULT

*Taken immediately before a forced suspend/restart. The four LIVE arms had not finished; all
artifacts are on persistent disk (H.3b).*

| arm | state | t | η_max | u_max | Newton |
|---|---|---|---|---|---|
| `c3v_base_mixed` | LIVE | **100.6 s (62.9 T)** | 0.1116 | 0.403 | 6 |
| `c3v_amp_mixed` | LIVE | 95.2 s | **0.2300** | **1.133** | **15** |
| `c3v_dt_mixed` | LIVE | 58.8 s | 0.1290 | 0.417 | 6 |
| `c3q_dt_mixed` | LIVE | 19.7 s | 0.1312 | 0.501 | 8 |

### ✅ `A` = 0.10 — a genuine post-fill stability result

```
 t= 40- 50  eta 0.1112  u 0.4023  NL 6      t= 80- 90  eta 0.1114  u 0.4112  NL 6
 t= 50- 60  eta 0.1112  u 0.4043  NL 6      t= 90-100  eta 0.1115  u 0.4127  NL 6
 t= 60- 70  eta 0.1112  u 0.4076  NL 6      t=100-110  eta 0.1116  u 0.4082  NL 6
 t= 70- 80  eta 0.1113  u 0.4097  NL 6
```

Flat in the **fourth decimal** over seven consecutive windows, Newton **constant at 6**, entirely
post-fill (the flume fills in ~36 s, rule 14). This satisfies rule 12c's bar in a way no earlier
`:full` run has: it is stability, not an unexpired clock.

### ⛔ `A` = 0.15 — diverges at ~87 s, same mode

```
 t= 60- 70  eta 0.1865  u 0.6949  NL  8
 t= 70- 80  eta 0.2113  u 0.7379  NL 10
 t= 80- 90  eta 0.2072  u 0.8514  NL  8
 t= 90-100  eta 0.2780  u 1.1443  NL 15    <- eta +34%, u +34%, Newton ~2x in one window
```

⚠ **THE CONCLUSION: THE MIXED FORMULATION *DELAYS* THE INSTABILITY, IT DOES NOT *REMOVE* IT.**
At `A` = 0.10 the delay exceeds the observation window and the run looks cured. At `A` = 0.15 the
same construction buys **≈16×** over its projected control (5.8 s → ~87 s) and then fails with the
identical velocity-led signature. **A pass at one amplitude is not a pass.**

⚠ **And note how nearly this was missed.** `c3v_base_mixed` alone — the single most quotable run on
the branch — would have produced a confident *"the projection was the cause, `:full` is fixed"*. It
took a **same-batch amplitude twin** to see that the mode was merely late. This is rule 14c and rule
38c's in-batch-control discipline paying out for the third time on this branch.

⚠ **This amplitude sensitivity is NOT old rule 13.** That cap (`A_wave ≤ 0.001`) was an equal-order
artefact, lifted 2026-09-06 when Taylor-Hood removed an inf-sup failure. This is measured **on**
Taylor-Hood at 100× that cap, in the Class-III path, and it is a **third independent axis** beside
pairing and `dx`: both `amp` arms die, each earlier than its `base` twin.

⚠ **Both `dt` arms were showing early warning at capture** — Newton maxima creeping 4 → 6
(`c3v_dt_mixed`) and 4 → 12 with `u` 0.434 → 0.707 (`c3q_dt_mixed`). Neither had failed; neither is
a pass.

⚠ **THESE FOUR ARMS WERE TERMINATED BY AN OPERATOR REBOOT AT 2026-09-23 15:38, NOT BY DIVERGENCE.**
A frozen GNOME session forced a restart; suspend was tried first and failed. Final states:
`c3v_base_mixed` **t = 101.4 s (63.4 T), η 0.1054, u 0.411, Newton 6 — still flat and healthy**;
`c3v_amp_mixed` t = 95.6 s (diverging, η 0.211, u 0.849, Newton 14); `c3v_dt_mixed` t = 59.2 s;
`c3q_dt_mixed` t = 19.8 s. **Do NOT read `c3v_base_mixed` stopping at 101.4 s as a failure** — its
diagnostics simply end where the machine went down. It never reached the 160 s / 100-period target,
so its result is a **lower bound**: ≥63.4 periods stable, not a completed 100-period run. Re-running
it is the cheapest way to convert that bound into the regression-gate reference trace §5.7 item 4
needs. All artifacts survived (`output/local_1d/<arm>/` + `_logs_2026-09-23/`).

## H.2 The four hypotheses and their verdicts

### H1 — "the frozen projection is what destabilises `:full`." → ✅ PARTIALLY CONFIRMED — **a delay, not a cure**

Mixed outlives projected in **8/8** matched pairs, by 1.9× to >16×, and produces the first
long-running `:full` configuration in the project's history (`A`=0.10, Q2/Q1: >100 s flat). The
projection is a large, real contributor. **It is not the whole cause** — H2–H4 survive its removal —
**and it is not sufficient even at Q2/Q1**: the `A`=0.15 twin diverged at ~87 s (H.1b).

### H2 — "Q3/Q2 will help, because Q2/Q1 under-represents Class III." → ⛔ REFUTED, SIGN REVERSED

This was option 0 of `OPEN_ISSUES.md` §0c and step (0) of `CLAUDE.md` §5.7 — *"free, do it first"*.
The reasoning was sound: at `η ∈ Q1`, `∂²η ≡ 0` **identically** in 1-D, so part of the Class-III
operator never contributes and Q2/Q1 was never an honest test.

**It is an honest test, and it fails harder.** All four Q3/Q2 arms die earlier than their Q2/Q1
twins — 16.0 vs >97.4, 9.6 vs >93.2, 6.4 vs 36.6, `dt` tracking the same way.

⚠ **The mechanism is NOT "more nonlinear effects."** The residual is byte-for-byte the same operator
at both pairings — same `:full`, same eight `𝓝` components. What changes is **representable
content**: with `η ∈ Q2`, `∇η` is piecewise **linear** rather than piecewise constant, so the
`∇η·u` half of `𝖲 = ∇·(Hu)` — the thing `𝖦` approximates the gradient of — becomes a real object for
the first time. (Note the mixed form never builds `∂²η` explicitly; the content still enters, through
`𝖲`.) **Consequence: the Q2/Q1 pass is a pass on a partially masked operator and must always be
quoted as such.** Generalised as `CLAUDE.md` rule 2d.

### H3 — "the Q3/Q2 penalty is just resolution." → ⛔ NOT PURELY, on a matched-DOF pair

Q3 on the same cells is 2.1× the DOFs, and rule 38b says more resolution advances onset here. The
batch contains an unplanned but decisive near-matched pair:

| run | pairing | nx | dx | DOFs | onset |
|---|---|---|---|---|---|
| `c3v_dx_mixed` | Q2/Q1 | 480 | 0.125 | 29781 | 36.6 s |
| `c3q_base_mixed` | Q3/Q2 | 240 | 0.25 | 31710 | **16.0 s** |

**At matched problem size, order still costs 2.3×.** ⚠ **Partial isolation only** — DOFs match but
`dx` does not (0.125 vs 0.25), so the grid scale, where a `λ ≈ 2–3·dx` mode would live, is free. The
clean run is **Q3/Q2 at nx = 120**, which matches `dx` = 0.5 against nothing yet run, or better, a
Q3/Q2 ladder read against the Q2/Q1 one at equal `dx`. Not done.

### H4 — "the `dx` signature is the broken-Hessian recovery error." → ⛔ REFUTED. **The key result.**

§5.2c attributed the `dx` sign to **recovering `∂²` of a `C⁰` field**: the exact distributional `∂²u`
is `{∂²u}` cell-wise **plus a Dirac layer on the skeleton** weighted by `[∂ₙu]`, and cell quadrature
sees only the first half — an error that does not improve as `h → 0`.

**The mixed formulation cannot commit that error.** `𝖦` is a genuine FE unknown defined by
`∫𝖦·Ψ = −∫𝖲 ∇·Ψ + ∮𝖲Ψ·n`; nothing is differentiated twice anywhere in the path.

**The `dx` signature survived it unchanged** — `c3v_dx_mixed` (nx=480) died at 36.6 s while
`c3v_base_mixed` (nx=240) passed 97.4 s: halving `dx` still advances onset by ≥2.7×.

⚠ **THIS IS A DEDUCTION FROM TWO RUNS, NOT A MEASUREMENT OF A MECHANISM.** At least one of the
following holds and they are **not yet separated**:

1. the `dx` signature was **never** the recovery error — it was something else all along;
2. `𝖦`'s own FE approximation of `∇𝖲` carries an `h`-dependent error of the same sign (it is an
   `L²`-type projection onto a finite space, so it has its own consistency error);
3. there is a genuine **grid-scale instability in the coupled `(η,u,𝖦)` system** — a pairing problem
   in the new `𝖦↔u` block, not in `η↔u`.

**Reading (3) has concrete support in the code:** `Vaux` is built from `reffe_U`
(`src/horizontal.jl:59`), so **`𝖦` sits in the velocity space at equal order with `u`** at both
pairings. `𝖦`'s own diagonal block is a Gram matrix and therefore coercive, so this is *not* the
classic rule-2b inf-sup failure — but the **coupled** pairing has never been analysed, and the
obvious probe (`𝖦` one order below `u`) has still never been run. It was already open as F.3.2.

⚠ **Immediate consequence for queued work: `C⁰`-IP loses its stated rationale.** Option 3 of
`OPEN_ISSUES.md` §0c is *"`C⁰` interior penalty if the `dx` signature survives 0–2"*. It survived —
but `C⁰`-IP exists to make a **broken `∂²` consistent**, and the surviving signature was measured
where there is no `∂²` to break. **Do not launch it on the old justification.** Generalised as
`CLAUDE.md` rule 39b.

### ⚠ H5 — "is Class III implicated at Q3/Q2 at all?" → **NOT TESTED. RULE 14c GAP.**

**Every Q3/Q2 arm in this batch is `BALFEM_NL_PRESSURE=full`, and the only `:native` runs on record
are Q2/Q1.** So the batch cannot distinguish:

* Class III is the carrier and Q3/Q2 is merely where it is finally represented → a `:native` Q3/Q2
  arm survives 100 periods; from
* Q3/Q2 is unstable here for a reason with nothing to do with Class III → `:native` dies too, and
  §0c, §5.7 item 0 and this whole branch are aimed at the wrong object.

**This is exactly the failure mode F.0 recorded three days earlier** — the harness batch whose
`:native` control was run last and died, voiding everything. *A crash is not evidence until its
control runs*, and the control is again the thing not run. **It is the single highest-value run
available and it is cheap.**

## H.3 Caveats that must travel with these numbers

* ⚠ **`c3v_amp_mixed` is UNDER WATCH, not a pass.** Its `u_max` window maxima run
  **0.695 → 0.738 → 0.851 → 0.951** over t = 60–100 with `η` bounded at 0.205–0.211 — the onset of
  the velocity-led signature. It may yet join the failures.
* ⚠ **"Running" is not "stable" (rules 12c, 14).** The flume fills in ~36 s; only post-fill behaviour
  counts, and only ≥100 periods earns the word. `c3v_base_mixed` at 97.4 s is 61 periods — good, not
  finished. ⚠ `c3v_dx_mixed`'s 36.6 s onset lands **exactly at fill completion**, so it deserves the
  rule-14 caveat even though the inflow pinning and η 0.113 → 0.153 make it a real divergence.
* ⚠ **Cost makes the Q3/Q2 tier impractical as built.** Those arms quoted ETAs of **186–242 h** for
  100 periods at ~31.7k DOFs with an AD-coupled 5-field Jacobian. None could have finished
  regardless of stability. Hand Jacobians for the auxiliary rows (F.3.3) are now a prerequisite for
  any Q3/Q2 stability *claim*, not an optimisation.
* ⚠ **Still a `p = 1` vertical basis, 1-DH, flat bed, one integrator.** Nothing here speaks to the
  variable-bed lee-shoulder mode (§0d), which is a separate defect.

## H.3b Where the outputs live, and the one way to lose them

Every arm directory under `output/local_1d/` holds **`diagnostics.csv`** (time series),
**`sol_t_*.vtu` + `solution.pvd`** (fields at 0.2 s), and **`run.log`** (config banner, per-step
Newton trace, crash reason). Index: `output/local_1d/README_c3_mixed_campaign.md`.

⚠ **`output/` IS GITIGNORED — none of it survives a fresh clone**, which is why every number that
matters is transcribed into this file and `CLAUDE.md` §5.2e rather than left in the CSVs.

⚠ **AND THE LOGS WERE ONE REBOOT FROM GONE.** The solvers write stdout into the session scratchpad
under **`/tmp`**, and this machine's `tmpfiles.d` carries **`D /tmp`** — *emptied on every boot*.
They were copied onto persistent disk on 2026-09-23 (`<arm>/run.log`, plus all 45 session logs in
`output/local_1d/_logs_2026-09-23/`, which is the only copy of the isolation arms `iso_gs`/`iso_gb`,
the `c3_0`/`c3_2` factorial and the unit gates). **Any future batch should redirect stdout to
`output/…` directly rather than to the scratchpad.**

⚠ **The `diagnostics.csv` does NOT record why a run stopped** — it simply ends. The reason is only in
`run.log`, and its traceback is misleading: every one terminates inside `paraview_collection` /
`createpvd`, which is the enclosing `do`-block, **not** a VTK failure. Use
`grep -m1 '^ERROR' run.log`.

⚠ **The Q2/Q1 `base` projected control is not in this batch** — it is the earlier reference run
`c3_0_Q2Q1_lagged` (t = 12.600000, `u_max` = 2.50462190), which is what every "12.6 s" claim cites.

## H.4 Next steps, ordered by what would change a conclusion

1. **`:native` at Q3/Q2** — closes H5, the rule-14c gap. Gates everything below it. Cheap.
2. **`𝖦` one order below `u`** — the only probe of reading (3) of H4; `Vaux` currently uses
   `reffe_U`. Open since F.3.2 and never run.
3. **Q3/Q2 at nx = 120** — completes H3's isolation by matching `dx` rather than only DOF count.
4. **Hand Jacobians for the auxiliary rows** — makes any Q3/Q2 claim affordable at all (H.3).
5. **Let `c3v_base_mixed` reach 100 periods** and re-read; it is the only candidate reference trace
   for a `:full` regression gate, and §5.7 item 4 needs one.
6. ⛔ **`C⁰`-IP is NOT next** — H4 removed its rationale. Re-derive a justification before spending
   on it.
