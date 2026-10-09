# STABILISATION_ORDER_RULE.md — consistency terms against stabilisation terms, and which jump orders a penalty must reach

*Opened 2026-10-06 (branch `v2-solver`). It records a line of reasoning and the tests that bear on
it. The v1 measurements quoted here were made with component 4 of `𝓝` omitted (`C3_MASK=gs`) and
are labelled (v1). Their detailed records are in [`BROKEN_FORMULATION_PLAN.md`](BROKEN_FORMULATION_PLAN.md)
§4–5 and [`GHOST_PENALTY_PLAN.md`](GHOST_PENALTY_PLAN.md) §5. LaTeX: chapter 6 (the broken audit),
chapter 9 (stabilisation), appendix C (the jump penalties).*

---

## 1. Two different jobs for skeleton terms

On `C⁰` Lagrange spaces an element-wise integration by parts leaves integrals on the interior facets.
These can play two quite different roles, and confusing them cost the project a year (CLAUDE.md
rule 39b).

**(a) Consistency terms.** They complete the weak form, so that the exact solution satisfies the
discrete equations (Galerkin orthogonality), or so that the discrete operator is the operator
actually applied to the discrete field. The broken audit (LaTeX §6.1.5) sorts every facet term by
which quantity jumps:
* **kind 0**: the jump of a `C⁰` field; it vanishes identically and the term does not exist;
* **kind S**: a jump of the trial solution. It vanishes at the exact solution, so including or
  omitting it gives two different but equally consistent operators;
* **kind T**: a jump of a test-function derivative. It does not vanish at the exact solution and is
  mandatory.

In this model:
* **No kind-T term arises** in the residual as written. The leading-pressure block would need one
  only in its integrated-by-parts form (B); the solver uses form (A).
* **The broken Class-III layer** (`broken_class3_residual`) is the kind-S term that makes the
  operator distributionally complete. It is always assembled with `nl_pressure=true`.
* **The leading-pressure jump `⟦𝒫ᵢ⟧`** is kind S and **excluded**: it makes the effective mass matrix
  `h⁻⁴`-conditioned and the scheme non-convergent.

None of these terms has a sign in the energy, so none removes energy. **Consistency terms are not
stabilisers.**

**(b) Stabilisation terms.** They are added after the weak form is complete, against spurious modes
that the consistent operator lets grow. They are sign-definite, vanish on functions smooth across
the facets (so they stay consistent), and are not needed for consistency at all.

**The classical "C⁰ interior penalty" does both jobs.** In the C⁰-IP method for fourth-order
problems (Engel et al. 2002; Brenner & Sung 2005) the consistency terms and their symmetric
counterpart make the method consistent and adjoint-consistent. A *separate* penalty
`σ h⁻¹ ⟦∂ₙu⟧⟦∂ₙv⟧` makes it coercive.

In this solver the two roles are split between two pieces of code:
* **Consistency:** the broken layer (`src/broken.jl`, form A), no penalty involved.
* **Stabilisation:** `stab_contrib`, either `:jumpgrad` (jumps of `∂ₙʲ`, `j ≤ cip_order ≤ 2`) or
  `:ghostvolume` (the direct ghost penalty, the complete jump family `0…p`). One of the two per run,
  never both.

So **`:jumpgrad` with order 1 is only the penalty half of a classical C⁰-IP.** The consistency half
is the broken layer, and it was there before any penalty was tried.

Moreover, the instability being fought is not a coercivity defect of an elliptic problem. It is a
transfer of carrier energy towards the grid scale by the nonlinear terms, into a band that `C⁰`
elements represent with no null and no dissipation (LaTeX chapter 8). The penalty is chosen against
that, not against a loss of coercivity.

---

## 2. What a jump penalty leaves undamped: its kernel

A penalty on the normal-derivative jumps of orders `1…m` vanishes exactly on the piecewise
polynomials that are `C^m` across the facets. That set is its **kernel**. Oscillatory modes that fit
in the kernel are not damped, and in the linear eigen-analysis they rearrange into it.

The linear eigen-analysis showed this directly for the first-order penalty (v1, Q3/Q2,
16 cells/λ). The wave modes moved into the `C¹` subspace and stayed almost undamped. Their damping
was at most 0.016 / 0.003 / 0.0003 s⁻¹ at γ_u = 0.01 / 0.1 / 1, i.e. it **falls** like `1/γ`. The
penalty acts as a `C¹` *constraint*, not as a band-limited damper.

For `Q_p` in one horizontal dimension:

| penalty reaches order `m` | kernel | wave content per cell | consequence |
|---|---|---|---|
| `m ≤ p − 2` | `C^m` splines | ≥ 2 | still carries the element-scale (optical) branches beyond the cell Nyquist |
| **`m = p − 1`** | **maximally smooth `C^{p−1}` splines** | **1** | resolves up to the cell Nyquist; no optical branch |
| `m = p` (ghost penalty) | one polynomial over the patch, i.e. globally | 0 | no undamped wave at all, so γ must stay inside a window or the carrier locks |

The `m = p − 1` row links to a known result. Maximally smooth splines have no optical branches in
their discrete spectrum (Cottrell, Reali, Bazilevs & Hughes 2006, *Isogeometric analysis of
structural vibrations*, CMAME 195). Those branches are exactly the `C⁰` content beyond the node
Nyquist that LaTeX chapter 8 identifies as where the transferred energy accumulates.

---

## 3. Two candidate rules for how far the penalty must reach

* **(W) — the order of the weak form.** The weak form carries trial second derivatives (rule 1b:
  `∂³u` in the strong form, `∂²u` after the leading-pressure integration by parts), so penalise up to
  order 2 for every field. Read as a conformity argument, (W) actually says less: a second
  derivative is an `L²` function, with no skeleton layer, as soon as the field is `C¹`. So `m = 1`
  would suffice for a conforming `H²` weak form, which is the classical C⁰-IP reasoning.
* **(P) — the element order.** Penalise up to `m = p − 1` per field, so the kernel is the maximally
  smooth spline space: the richest space with no optical branch, and no locking.

**What the evidence says (v1, closed box, Crank–Nicolson, Q3/Q2 = velocity `p = 3`, surface `p = 2`):**

| penalty | velocity `m` vs `p − 1` | surface `m` vs `p − 1` | result |
|---|---|---|---|
| none (broken, γ = 0) | — | — | Newton failure at ≈ 44 s, every band from `5k₀` up at +0.16…+0.29 s⁻¹ |
| first order, (γ_u, γ_η) = (0.3, 0.3) | 1 < 2 | 1 = 1 | grid-scale growth removed; **late mid-band growth from ≈ 110 s at +0.13 s⁻¹** (λ ≈ 4–10 nodes) |
| first order, velocity only | 1 < 2 | none | delays onset (≈ 50 → 80 → 110 s at γ_u = 0, 0.03, 0.3, Q2/Q1, 32 cells); no effect on `:native` |
| jump order ≤ 2, γ = 2e-3 (both fields) | 2 = 2 | 2 = p (full) | ✅ 160 s, every band decaying; carrier −1.5e-4 s⁻¹ |
| ghost, γ ∈ [0.0033, 0.03] (both fields) | 3 = p (full) | 2 = p (full) | ✅ 160 s, decaying; carrier loss 0.5 / 1.5 / 4.3 % |

How each rule fares against these results:
* **The first-order failure refutes the conformity reading of (W).** `m = 1` makes the weak form
  `H²`-conforming, yet the run grew.
* **(W) and (P) coincide for the velocity at Q3.** `p − 1 = 2` is also the weak-form order, so the
  order-≤2 success cannot tell them apart.
* **The surface did the decisive work in the order-≤2 success.** For the Q2 surface, order ≤ 2 is
  already the *full* family (`m = p`). In the eigen-analysis, velocity-only orders ≤ 2 left mid-band
  modes undamped; the window existed only with the surface penalised, whose kernel is then trivial.
* **Q2/Q1 has no window at any γ or order (v1 eigen-analysis).** A Q1 surface admits only `m = 1 = p`,
  i.e. locking, and its `p − 1 = 0` means no surface penalty at all.

---

## 4. The discriminating test: Q4/Q3 with orders 1 and 2

At **Q4/Q3** (velocity `p = 4`, surface `p = 3`) the rules part:

| rule | velocity, order ≤ 2 | surface, order ≤ 2 | prediction |
|---|---|---|---|
| (W) | sufficient | sufficient | every mid/high/sub-element wave mode damped |
| (P) | insufficient (needs 3) | sufficient (`2 = p − 1`) | **undamped modes left in the `C²` quartic splines** (two wave dofs per cell) |

Order 3 cannot be built with `:jumpgrad`, because Gridap evaluates basis derivatives only up to
order 2. Order ≤ 2 is exactly what the test needs, though.

**Run** `examples/local_1d/cip_eigen_order_rule.jl`, launched 2026-10-06:
* the linearised flat-bed box operator `−M_eff⁻¹(K + J_h)`, diagonalised densely, 16 cells/λ,
  P1LFE-2, `d` = 3.5 m;
* Q4/Q3 and, as the same-batch control, Q3/Q2;
* `:jumpgrad` orders 1 and 2, penalty on both fields and on the velocity only;
* γ ∈ {1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2, 3e-2};
* the ghost penalty at Q4/Q3 as the complete-family reference.

**Output, revisitable without re-running:** `output/local_1d/cip_eigen_order_rule/`
* `summary.csv`: one row per case, with per-band counts, minimum and maximum damping, and the
  frequency of the least-damped mode;
* `modes_<case>.csv`: every eigenvalue;
* `run_manifest.toml`;
* `run.log`.

**How to read it.** The linear operator has no instability, so these are dissipation budgets, not
stability verdicts. The criterion is the one used in v1:
* the weakest x-varying mid-band mode is damped above the late growth rate (≈ +0.13 s⁻¹), with margin;
* the carrier is damped at ≪ 1/T_run ≈ 0.006 s⁻¹.

(P) is supported if Q4/Q3 with order ≤ 2 on both fields leaves mid- or high-band modes at ≈ 0 damping
at every γ, while the ghost penalty does not. (W) is supported if order ≤ 2 opens a window at Q4/Q3
as it did at Q3/Q2. A closed-box nonlinear run at the chosen γ would then confirm either reading.

### Results

**First pass (2026-10-06, `output/local_1d/cip_eigen_order_rule/`): inconclusive, because the mode
classifier was blind to transverse modes.**
* **The operator is right.** The Q3/Q2 control reproduces the v1 carrier damping of the order-≤2
  penalty exactly: 7.61e-5 / 1.50e-4 / 2.25e-4 s⁻¹ at γ = 1e-3 / 2e-3 / 3e-3, against v1's 7.6e-5 /
  1.5e-4 / 2.25e-4.
* **The mid-band minimum is wrong.** The minimum damping of the "mid band" was ≈ 0 in *every*
  case, including the ghost penalty at Q4/Q3 and the Q3/Q2 order-≤2 control, where v1 recorded
  0.27–0.82 s⁻¹. The modes responsible all sit at ω ≈ 7.55–7.83 rad/s. That is the frequency of the
  **transverse** standing mode of the one-cell-wide box and of the saturated short waves.
* **Cause.** The classifier samples only the centre line `y = h/2`. A mode `∝ cos(πy/h)·cos(kx)`
  vanishes there, so the "x-varying" test (`x_frac > 0.5`) passes on what little it sees.
* **Why these modes must be excluded.** They are carried by the y-structure inside the single cell
  (`Q4` has three interior y-node rows where `𝖴y` is free, `Q3` has two). They cannot be excited by a
  y-invariant run, and they cannot be reached by any x-facet penalty, because there are no interior
  y-facets.
* **Consequence.** The first-pass mid-band minima discriminate nothing.
* **⚠ The v1 window values were partly luck.** The v1 file (`output_v1/local_1d/cip_eigen/cip_eigen.csv`,
  same classifier) gives, for Q3/Q2 order ≤ 2:
  * 0.27 / 0.55 / 0.82 at γ = 1e-3 / 2e-3 / 3e-3, but
  * **−6e-14 at γ = 0.01**, breaking its own monotone trend.

  The transverse modes are degenerate, so LAPACK returns an arbitrary mixture of them and of
  x-waves at the same frequency. Whether a mixture passes the centre-line filter depends on that
  arbitrary choice, which differs between runs (BLAS threads, version). The v1 γ window
  (`BROKEN_FORMULATION_PLAN.md` §5.8, LaTeX appendix C) must therefore be re-read with the
  y-filter before it is trusted.
* **The v1 closed-box nonlinear results do not depend on this classifier.** They used band energies
  of y-invariant runs, in which transverse modes are never excited.

**Second pass (`output/local_1d/cip_eigen_order_rule/yfilter/`, `EIG_OUT`).** It adds a y-structure
test: η and every `𝖴x` component sampled on `y = h/4` and `3h/4`. A mode is classed `transverse` if
its y-odd energy share or its `𝖴y` energy share exceeds ½. The `modes_*.csv` files carry both
fractions, so the threshold can be revisited without re-running.

**Second-pass results (2026-10-06).** The filter works:
* 34–58 transverse modes per case are now excluded;
* every remaining least-damped mid-band mode is a genuine x-wave (y-odd share 0, `𝖴y` share 0).

The table gives the energy damping (s⁻¹) of the weakest x-varying **mid-band** mode, 16 cells/λ:

| case | γ = 1e-4 | 3e-4 | 1e-3 | 2e-3 | 3e-3 | 1e-2 | 3e-2 | reading |
|---|---|---|---|---|---|---|---|---|
| Q3/Q2, order 1, both fields | 3.0e-7 | 9.1e-7 | 3.0e-6 | 5.9e-6 | 8.6e-6 | 1.9e-5 | 1.5e-5 | escape (`C¹` kernel), as in v1 |
| Q3/Q2, order ≤ 2, velocity only | 3.0e-7 | 9.1e-7 | 3.0e-6 | 5.9e-6 | 8.6e-6 | 1.9e-5 | 1.5e-5 | escape: the surface must be penalised |
| **Q3/Q2, order ≤ 2, both fields** | 0.049 | 0.13 | **0.37** | **0.65** | **0.91** | **2.8** | 8.3 | **window**, monotone in γ |
| Q4/Q3, order 1, both fields | 2.5e-7 | 6.9e-7 | 4.5e-7 | 2.3e-7 | 1.5e-7 | 4.5e-8 | 1.5e-8 | escape |
| Q4/Q3, order ≤ 2, velocity only | 2.5e-7 | 7.1e-7 | 5.1e-7 | 2.6e-7 | 1.7e-7 | 5.1e-8 | 1.7e-8 | escape |
| **Q4/Q3, order ≤ 2, both fields** | 3.2e-6 | 1.5e-6 | **5.0e-7** | **2.5e-7** | **1.7e-7** | **5.1e-8** | 1.7e-8 | **escape: no window**; damping falls like `1/γ` |
| Q4/Q3, ghost (order `p`), both fields | 2.7e-3 | 5.2e-3 | 6.0e-3 | 8.3e-3 | 0.011 | 0.034 | 0.10 | no escape (rises with γ), but weak at these γ |

The carrier damping at Q4/Q3 is below 1e-6 s⁻¹ in every case (ghost at γ = 0.03: 5.0e-7).

**Reading.**
* **(W), the weak-form-order rule, is refuted at the linear level.** At Q4/Q3 a penalty reaching
  order 2, the order of the weak form, on *both* fields leaves an x-varying mid-band mode (n = 8,
  ω = 7.55 rad/s) undamped at every γ. Its damping *falls* like `1/γ`, the signature of a mode
  retreating into the penalty's kernel, exactly as the first-order penalty behaves at Q3/Q2. At
  Q3/Q2 the same order-≤2 penalty opens a clean window. The difference between the two pairings is
  the velocity order: at Q3, order 2 is `p − 1`; at Q4 it is `p − 2`.
* **(P), the `p − 1` rule, is consistent with every case, but not yet proven sufficient.**
  * Every escape happens with the velocity penalised below `p − 1`, or with the surface unpenalised.
  * The one window (Q3/Q2, order ≤ 2) has the velocity at `p − 1` and the surface at `p`.
  * Whether `m = p − 1` on *both* fields suffices (velocity order 3 at Q4, surface order 1 at Q3)
    cannot be tested with `:jumpgrad`: Gridap's order-2 limit at Q4, and one shared `cip_order`
    for both fields.
* **The ghost penalty (order `p`) has no escape at Q4/Q3,** but its weakest mid mode is a
  well-resolved n = 5 wave (≈ 13 nodes per wavelength). Its damping is small at γ ≤ 0.01 (0.034 at
  γ\* = 0.01, against 0.85 at Q3/Q2 in v1) and reaches 0.10 at γ = 0.03. Because the carrier is
  barely touched (≤ 5e-7), γ can be raised much further at Q4/Q3 than at Q3/Q2. **γ must be
  calibrated per pairing**, as the dimensionless scaling already implied.
* **The v1 window is confirmed and corrected.** Q3/Q2 order ≤ 2: 0.37 / 0.65 / 0.91 at
  γ = 1e-3 / 2e-3 / 3e-3 (v1 recorded 0.27 / 0.55 / 0.82), and **2.8 at γ = 0.01**, where v1's
  zero was the classifier artefact. The window is therefore wider than v1 stated: it extends at
  least to γ = 0.01, at a carrier cost of 7.5e-4 s⁻¹ (11 % of the energy over 160 s).

**What this establishes, and what it does not.** On the linear box operator, the jump order a
penalty needs is set by the **element order**, not by the order of the weak form. This is the
kernel argument of §2, observed. Two things remain:
* the nonlinear closed-box run, which decides stability (the linear operator has none to show);
* a test of whether `m = p − 1` on every field suffices, which needs a penalty that stops at
  `p − 1` without derivatives. A truncated ghost penalty, dropping the order-`p` part of
  `G_jk` in (★), would do it.

The next step would be the closed box at Q4/Q3 (nonlinear pressure on, Crank–Nicolson) with
`:jumpgrad` order 2 against the ghost penalty. The prediction from (P) is that order 2 fails late,
as the first-order penalty did at Q3/Q2, and the ghost penalty holds.

---

## 5. Why it matters

If (P) holds, the right penalty is not the ghost penalty (order `p`, trivial kernel, a γ window
against locking). It is a penalty that stops at `p − 1`, whose kernel is the rich, optical-branch-free
maximally smooth spline space. That is the "locking-free" design of the ghost-penalty literature
(Burman, Hansbo & Larson), which needs a deliberately rich kernel. Such a penalty would:
* need no upper limit on γ;
* allow a Q1 surface to go unpenalised (`p − 1 = 0`);
* reduce the stabiliser to a statement about the element, independent of the model.

If (W) holds, the needed order is a property of the PDE, and `:jumpgrad` with order 2 suffices for any
pairing. Gridap's derivative limit would then never bind.

---

## 6. Related v2 evidence (unstabilised box ladder, Crank–Nicolson, 2026-10-06)

The model without `𝓝` (`nl_pressure=false`; the v1 twin was the partial `:native`), t = 60–160 s, mid
band, surface velocity:
* Q3/Q2, 16 cells: +0.031 s⁻¹;
* Q2/Q1, 32 cells: +0.024 s⁻¹;
* **Q2/Q1, 64 cells: +0.08 s⁻¹ (×2 100).**

All three are bounded at 160 s. The rate roughly triples from 32 to 64 cells, which is the grid-scale
signature (rule 38b) in a model with no nonlinear pressure at all. It is 2.5–3× slower than v1
`:native`. The full-pressure arms are still running. Record: LaTeX chapter 8, §"Results of the v2
campaign".
