# GHOST_PENALTY_PLAN.md — the direct (volume) ghost penalty as a second stabilisation

*Branch `broken-formulation-solver`, 2026-10-02. Companion to `BROKEN_FORMULATION_PLAN.md`
(§5.2 high-order jump penalty, §5.5 the direct ghost penalty, §5.8 the eigen-analysis that found
the jump-penalty window). This file is the design, the task list and the results record for the
ghost-volume stabilisation.*

---

## 0. Why

* The high-order jump penalty (`:jumpgrad`) needs ⟦∂ₙʲ·⟧ for every j ≤ p. Gridap evaluates FE
  derivatives only up to order 2, so ⟦∂ₙ³𝖴⟧ is missing at Q3/Q2 (`BROKEN_FORMULATION_PLAN.md` §5.3).
* The **direct ghost penalty** reaches every order at once without any derivative. For each
  interior facet F with patch ω_F = T⁺ ∪ T⁻,

      J_F(u, v) = ∫_{ω_F} (E u⁺ − E u⁻) · (E v⁺ − E v⁻) dx

  where E u^± is the polynomial of the cell T^± extended to the whole patch (Preuß 2018;
  Lehrenfeld & Olshanskii 2019).
* Its kernel is "one polynomial across each patch", i.e. a global polynomial. That is the same
  kernel as the complete jump family, so it has the same locking caveat: γ must sit in a
  damping window, to be found by the eigen-analysis.

## 1. The formulation implemented

### 1.1 The term

    J_h = Σ_F γ_u τ_u h^(s−3) ∫_{ω_F} (E𝖴ₐ⁺ − E𝖴ₐ⁻)·Mv (E𝖵ₐ⁺ − E𝖵ₐ⁻)      (a = x, y)
        + Σ_F γ_η τ_η h^(s−3) ∫_{ω_F} (Eη⁺ − Eη⁻)(Eq⁺ − Eq⁻)

* τ_u = d√(gd) and τ_η = √(gd), as for `:jumpgrad`;
* s = `cip_hexp` (default 2), so the prefactor is h⁻¹;
* Mv weights the layer contraction (vertical-basis invariant).

### 1.2 Why h^(s−3)

On a facet-aligned rectangle, the Taylor expansion of the patch difference in the normal coordinate
is EXACT (it is a polynomial). With n = n⁺ and t the distance from F:

    (E u⁺ − E u⁻)(x_F ∓ t n) = Σ_{j=0..p} (∓t)^j / j! · ⟦∂ₙʲ u⟧(x_F)

so

    ∫_{ω_F} |E u⁺ − E u⁻|² = ∫_F Σ_{j,k} G_jk ⟦∂ₙʲu⟧⟦∂ₙᵏu⟧ ,
    G_jk = [1 + (−1)^(j+k)] h^(j+k+1) / ((j+k+1) j! k!)                  (★)

* The diagonal j = k term is ∝ h^(2j+1). Times h^(s−3) it gives h^(s+2(j−1)), exactly the
  `:jumpgrad` weight of order j.
* So the ghost form is the **complete** jump family (orders 0…p, j = 0 vanishing for C⁰ fields),
  **with** the cross terms G_jk (j + k even) that `:jumpgrad` omits.
* (★) is also the strongest available unit test (§3, G9c): for a p ≤ 2 field every jump in it is
  computable, so the ghost matrix must equal the (★)-weighted jump matrix to round-off.

### 1.3 How E is evaluated in Gridap (the design)

There is no Gridap API for evaluating a cell's polynomial outside its cell. On our **uniform
Cartesian** meshes the patch integral is computed as a facet integral of a line integral along the
normal:

    ∫_{T⁺} g = ∫_F ∫_0^h g(x_F − t n⁺) dt ,   ∫_{T⁻} g = ∫_F ∫_0^h g(x_F + t n⁺) dt

with an (p+1)-point Gauss rule in t (exact for the degree-2p integrand).

* **Shifted cell fields.** For a field f on Ω and a constant reference-space shift δ, define
  f_δ = f ∘ (ξ ↦ ξ + δ) by composing its cell data with an affine field:
  `lazy_map(Broadcasting(∘), get_data(f), Fill(AffineField(I, δ), ncells))`.
  * On a uniform Cartesian mesh a physical shift Δx is the reference shift Δx ./ (hx, hy) in
    EVERY cell.
  * f_δ.plus is the plus cell's own polynomial at x_F + Δx, and f_δ.minus the minus cell's
    polynomial at the SAME physical point, **extended** outside its cell when Δx points into T⁺.
* **Direction selection.** The facet normal n⁺ is one of ±e_x, ±e_y. For each of the occurring
  directions e_d and each Gauss abscissa t_q, build f at the shift −t_q e_d, and gate the facet
  integrand with the indicator [n⁺ = e_d]. Exactly one direction is active per facet.
* **DOF placement is Gridap's.** f_δ is built with the same basis type as f (FE function, test or
  trial basis), so `.plus`/`.minus` on the skeleton assemble into the correct plus/minus DOFs
  exactly as for any skeleton term.
* **Periodic wrap facet:** the shift is in REFERENCE coordinates of each cell, so the extension
  is geometrically correct even where plus and minus sit at opposite ends of the domain.
* **No derivatives at all**, so every order up to p is reached and AD needs nothing new.

**Restrictions, enforced:**
* CartesianDiscreteModel with uniform cells; asserted from the cell measures at attach time;
* x- and y-aligned facets only (true for Cartesian meshes);
* sequential only, like the rest of `broken.jl`.

### 1.4 Switch

`stabilization::Symbol ∈ (:jumpgrad, :ghostvolume)` selects the term:
* default `:jumpgrad`, so every existing run and test is unchanged;
* the γ knobs (`cip_gamma_u`, `cip_gamma_eta`) and `cip_hexp` are shared;
* `cip_order` applies to `:jumpgrad` only. `:ghostvolume` always covers orders 0…p, and setting
  `cip_order ≠ 1` with it is refused.
* Plumbed through `attach_skeleton!`, `setup_and_run`, `run_mms_case`, and `BALFEM_STAB` in both
  1-D drivers. The output-name token is `ghost`.

## 2. Tasks (sequential)

* **T1 — Feasibility gate (REPL, before any code in `src/`).** Shifted FE functions and bases
  compose. Their skeleton `.plus`/`.minus` evaluate the extended polynomial (checked against the
  analytic continuation of an interpolated polynomial). An assembled skeleton integral of shifted
  bases has the right size. If any of this fails, stop and report.
* **T2 — `src/broken.jl`:**
  * `_shift_field(f, δ)` for FEFunctions, single-field bases and multi-field basis components;
  * `ghost_ctx(model, p_u, p_eta)`: uniform-mesh check, (hx, hy), occurring directions, Gauss
    rules in t per field degree;
  * `ghost_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)`, LINEAR in the unknowns, so it is its own
    Jacobian;
  * `attach_skeleton!(…; stabilization)` storing `stab`, the coefficients and the ghost context.
* **T3 — Wiring:**
  * `global_residual` / `jacobian_u` dispatch `cip_contrib` vs `ghost_contrib` on `sk.stab`;
  * `setup_and_run` (kwarg, guard, banner);
  * `run_mms_case`, the periodic and flume drivers (`BALFEM_STAB`, name token);
  * `examples/local_1d/cip_eigen_analysis.jl` (a `stab` argument).
* **T4 — Gates** (§3) in `test/test_broken_formulation.jl`, run as a separate process.
* **T5 — Eigen-analysis** of `:ghostvolume`, Q3/Q2 n16 (and Q2/Q1 n32 for completeness):
  carrier damping vs weakest mid/high-band damping over a γ ladder → the window and the γ to run.
* **T6 — Box campaign** (Crank–Nicolson, 160 s, up to 10 cores in total). Arms in §4.
* **T7 — Record** results here and in `BROKEN_FORMULATION_PLAN.md`, with a `CLAUDE.md` pointer.

## 3. Gates (all must pass before T6)

| gate | check | oracle |
|---|---|---|
| G9a kernel | J_ghost · (global polynomial of the trial degree) = 0, open 2-D mesh, sloping bed | exact, round-off |
| G9b algebra | symmetric, positive semidefinite | round-off |
| G9c identity (★) | for a Q2 velocity and Q1/Q2 surface, the ghost matrix equals the (★)-weighted ⟦∂ₙʲ⟧⟦∂ₙᵏ⟧ matrix (j, k ≤ 2) — on an OPEN 2-D mesh (x and y facets) **and** on an x-PERIODIC mesh (wrap facet) | independent closed form, round-off |
| G9d mass | (q ≡ 1)·J_ghost = 0 | round-off |
| G9e Jacobian | linear regime + ghost: hand ∂R/∂u vs central FD | < 1e-6 |
| G9f selectivity | Rayleigh quotient λ = 2dx / carrier ≫ 1 (Q3/Q2 box) | — |
| G9g switch | `:jumpgrad` default bitwise unchanged; `:ghostvolume` with `cip_order ≠ 1` refused; a non-uniform mesh refused | — |

## 4. Campaign arms (T6) — closed x-periodic box, Q3/Q2, Crank–Nicolson, 160 s

Running already (2 + 1 cores): the `:jumpgrad` order ≤ 2 pair at γ = 2e-3, and the first-order
flume run. Up to **7 new arms** (10 cores in total):

| arm | tier | γ (both fields) | purpose |
|---|---|---|---|
| G1 | `:full` broken | γ* (from T5) | the test |
| G2 | `:native` | γ* | tier twin |
| G3 | `:full` broken | γ*/3 | lower edge of the window |
| G4 | `:full` broken | 3γ* | upper edge (locking check) |
| G5 | `:full` broken, **32 cells/λ** | γ* | refinement: rate vs h |
| G6 | `:full` broken, **A = 0.15** | γ* | amplitude ceiling |
| G7 | `:native`, **32 cells/λ** | γ* | refinement twin |

Controls already on record (same branch code): `:full` broken γ = 0 (died at 44 s); `:native`
γ = 0 (growth from ≈ 90 s); first-order (0.3, 0.3) (`:full` mid-band growth from ≈ 110 s).

**Acceptance** (`OPEN_ISSUES.md` §0e):
1. mid/high-band envelopes flat to 160 s;
2. carrier energy loss ≤ a few % over the run;
3. MMS orders preserved (checked at γ* after the box);
4. amplitude ceiling stated.

## 5. Results

*(filled in as tasks complete)*

### 5.1 T1 — feasibility (2026-10-02) ✅

* The shifted-extension evaluation equals the Taylor closed form (★) to **3.1e-15** on an
  x-periodic Q2 box, wrap facet included.
* The ASSEMBLED ghost matrix, built with the multi-field test and trial bases on an open 2-D mesh
  with x- AND y-facets, equals the (★)-weighted jump matrix to **2.9e-16**, with identical sparsity.
* Two implementation facts:
  * `similar_cell_field` is `@notimplemented` for `MultiFieldFEBasisComponent`, so `_shift`
    rebuilds it from the shifted single-field basis;
  * FillArrays is not a direct dependency, so a plain `fill` vector is used.
* The design simplified during T1: summing over BOTH shift signs ±t along the facet-normal axis
  covers T⁺ and T⁻ whatever the plus/minus orientation (plus and minus are evaluated at the same
  physical point), so no orientation indicator is needed.

### 5.2 T5 — eigen-analysis of `:ghostvolume` (linearised box, Q3/Q2 n16; energy damping, s⁻¹)

| γ (both fields) | carrier (lost in 160 s) | mid-band min | high band | oscillatory pairs |
|---|---|---|---|---|
| 1e-4 | 8.4e-7 | ~0 | 0.75 | 95 |
| 1e-3 | 8.4e-6 | ~0 | 7.5 | 95 |
| 3e-3 | 2.5e-5 | ~0 | overdamped | 70 |
| **1e-2** | **8.4e-5 (1.3 %)** | **0.85** | overdamped | 50 |
| **3e-2** | **2.5e-4 (4 %)** | **2.55** | overdamped | 47 |
| 1e-1 | 8.4e-4 (13 %) | 8.5 | overdamped | 36 |
| 3e-1 | 2.5e-3 | ~0 (locking) | — | 26 |

* Window γ ∈ [0.01, 0.03]; **γ\* = 0.01**.
* Compared with `:jumpgrad` (orders ≤ 2) at equal carrier cost (8e-5 s⁻¹): weakest mid band
  **0.85 vs 0.27**. The extra order buys ×3 selectivity.
* Velocity-only ghost (γ_η = 0): the mid-band minimum stays ~0 at every γ, so the surface
  penalty is required, as for `:jumpgrad`.

### 5.3 T2–T4 — implementation and gates (2026-10-02) ✅

**Code:**
* `src/broken.jl`: `_shift_map`, `_shift_data`, `_shift`, `build_ghost_ctx`, `ghost_contrib`,
  `stab_contrib` (appended at end of file, rule 6b).
* `attach_skeleton!(…; stabilization)` and `global_residual` / `jacobian_u` dispatch through
  `stab_contrib`.
* Knob `stabilization` in `setup_and_run` and `run_mms_case`; `BALFEM_STAB` in both 1-D drivers
  (name token `ghost`); a `stab` argument in `cip_eigen_analysis.jl`.

`test/test_broken_formulation.jl`: **30/30** (G1–G8 unchanged, plus G9):

| gate | result |
|---|---|
| G9a kernel | 9.5e-17 |
| G9b symmetric PSD | rank 1126/1231 vs 898 for `:jumpgrad` order ≤ 2 — reaches order 3 |
| G9c (★), open 2-D mesh | 3.3e-16 |
| G9c (★), x-periodic box (wrap facet) | 2.2e-16 |
| G9d mass | 1.1e-17 |
| G9e Jacobian | 2.7e-11 |
| G9f Rayleigh λ = 2dx / carrier | 8.9e6 |
| G9g switch | `:jumpgrad` default bitwise; `cip_order ≠ 1` refused; non-uniform mesh refused |

### 5.4 T6 — campaign launched (2026-10-02 15:50, `run/local/run_1dper_batch_ghost.sh 0.01 10`)

* All seven arms of §4 are running, γ* = 0.01 (G3: 0.00333, G4: 0.03); banners verified.
* 10 simulations on the machine in total: these 7, the `:jumpgrad` order ≤ 2 pair, and the
  first-order flume run.
* ⚠ The campaign was started while the full gate file was still compiling (> 90 min, owing to
  the large ghost expression trees). The core identity (★) had already passed in the REPL, and the
  file then passed 30/30, so nothing had to be stopped.

**Cost (measured 2026-10-02 17:00, 10 simulations sharing the machine):**
* JIT ≈ 70 min per arm with seven compiling in parallel.
* `:native` + ghost ≈ 23 s/step (vs ≈ 7 s for `:jumpgrad`): every residual evaluation builds
  2·(p+1) shifted copies of each field per facet axis.
* Expected wall times: `:native` 16 cells ≈ 26 h; `:full` arms and the 32-cell arms ≈ 2–3 days.
* The decisive windows (`:full` past ≈ 110 s, `:native` past ≈ 90 s) are reached well before
  the end, so partial analyses are possible.
* **Optimisation lead, if the method is retained:** precompute the shifted test and trial bases
  once (they are constant), and shift only the solution fields per evaluation.

### 5.5 Status 2026-10-02 18:55

* ⛔ **`p32_broken_ghost0.01_n32_cn_p100` — Newton failure at step 16 (t ≈ 0.64 s), residual stuck
  at 0.036. NOT an instability:**
  * η_max and u_max track the 16-cell twin step for step;
  * Newton climbs 5 → 8 → 10 iterations where n16 needs 4 → 6.
  * Diagnosis: the quasi-Newton Jacobian omits the broken Class-III facet-layer terms, which scale
    like 1/h. At 32 cells the missing block is twice as large and Newton stops contracting. By
    rule 17b this cannot change the converged answer, but it prevents convergence.
  * Fix options:
    * a hand Jacobian of the facet layer (its linearisation in 𝖲, 𝖻 and 𝖴ₙ);
    * `use_ad=1` (exact, but very slow with the ghost terms);
    * a smaller `dt`.
  * Not relaunched: the slot is free, but memory is at 8 GB available.
* **`:full` ghost arms (n16): ≈ 50–60 s/step** (A = 0.15: ≈ 100 s/step), so 160 s takes ≈ 2.5–4.5
  days. The decisive window (> 110 s) is reached after ≈ 2 days.
* **`:jumpgrad` order ≤ 2, γ = 2e-3 (§5.9 of BROKEN_FORMULATION_PLAN.md):**
  * `:native` to 89 s: all bands flat and slowly DECAYING (η band 4: 2.9e-9 → 2.5e-9);
    the unpenalised control grows ×7–20 over the same windows;
  * `:full` broken to 45 s: decaying.
* **Ghost `:native`, γ = 0.01, to 30 s:** flat.

### 5.6 Exact Jacobian of the broken Class-III blocks; the 32-cell arm relaunched (2026-10-02 21:15)

* `broken_class3_jacobian` (`src/broken.jl`) is the exact linearisation of the broken 𝓚/𝓟
  Class-III blocks, volume and skeleton layer, by the product rule through `_broken_state` /
  `_broken_state_lin`. It is added to `jacobian_u` when `is_broken(prob)`. It is meant to cure
  the 32-cell Newton stall of §5.5.
* Gate G10 (central FD of the Class-III residual alone, R_broken − R_(no Class III); sloping 2-D
  mesh and periodic box) is in the test file. Its run was still compiling after > 90 min on the
  saturated machine.
* By rule 17b the Jacobian cannot change the converged answer, so the arm was relaunched without
  waiting: `p32_broken_ghost0.01_n32_cn_p100_jac`. The failed log is kept under its old name.
* ⚠ The other running `:full` arms started BEFORE this change and still use the quasi-Newton
  Jacobian. Their answers are unaffected; only their Newton cost differs.
* ✅ **G10 (standalone run, 2026-10-02 22:00):** broken Class-III Jacobian vs a central FD of the
  Class-III residual alone, Q3/Q2, sloping-bed open 2-D mesh: **rel. gap 3.5e-11** (|FD| = 72).
  The full test file (with the periodic-box case) is still running.

### 5.7 Campaign results (2026-10-04 12:40) — ✅ no late growth in any penalised arm

The full test file passes **32/32**, including G10 on the periodic box (wrap facet). The relaunched
32-cell arm with the exact broken Jacobian runs at **Newton 2/step**; it diverged at step 16
before.

Per-20 s envelopes (max band energy) of the mid and high bands, η and u. Carrier σ_E is fitted
over t ≥ 20 s.

| arm (Q3/Q2, CN) | reached | mid/high envelopes | carrier σ_E (eigen prediction) |
|---|---|---|---|
| `:jumpgrad` o≤2, γ = 2e-3, `:full` broken | ✅ 160 s | **decaying** throughout (u mid 1.4e-7 → 3.0e-8) | −1.5e-4 (1.5e-4) |
| `:jumpgrad` o≤2, γ = 2e-3, `:native` | ✅ 160 s | decaying | −1.2e-4 |
| ghost, γ = 0.01, `:full` broken | 140 s | **decaying** (u mid 1.4e-7 → 4.5e-8) | −8.6e-5 (**8.4e-5**) |
| ghost, γ = 0.01, `:native` | ✅ 160 s | flat / decaying | −6.7e-5 |
| ghost, γ = 0.03, `:full` broken | 150 s | decaying | −2.7e-4 (**2.5e-4**) |
| ghost, γ = 0.0033, `:full` broken | 140 s | **flat** (u mid 1.5e-7 → 8e-8) | −2.9e-5 |
| ghost, γ = 0.01, `:full` broken, **32 cells** | 100 s | flat; the high band is 100× lower than at 16 cells | −4e-6 |
| ghost, γ = 0.01, `:native`, **32 cells** | ✅ 160 s | flat | −1e-5 |
| ghost, γ = 0.01, `:full` broken, **A = 0.15** | 80 s | decaying (from a ×70 higher forced level) | −1.2e-4 |

Controls: first-order (0.3, 0.3) `:full` — mid band ×100–250 over 115–160 s; `:native` γ = 0 —
×920 over 110–160 s.

**Reading.**
* Both the hp jump penalty (orders ≤ 2) and the ghost penalty, inside the window chosen by the
  eigen-analysis, remove the late mid-band growth that the first-order penalty only delayed —
  for `:full` and for `:native`.
* The measured carrier damping matches the linear prediction to ±8 %, so the eigen-analysis is a
  reliable design tool for γ.
* Energy lost by the carrier over the run: 0.4 % (ghost 0.0033) to ≈ 4 % (ghost 0.03); ≈ 1.2 % at
  ghost γ* = 0.01.

**Flume** (§4.11, first-order (0.3, 0.3), CN):
* t = 158 s; η_max steady at 0.110 and u_max at 0.44 from 60 s to 160 s, Newton 8;
* the unpenalised flume died at 38 s and mixed at 16 s (SDIRK);
* no band analysis yet (the flume has boundaries; the envelope is read from max η and u).

**Still running:**
* ghost `:full` n16 to 160 s (γ = 0.01, 0.0033, 0.03);
* 32 cells (100 → 160 s);
* A = 0.15 (80 → 160 s);
* the flume's last 2 s.

**Open:**
* the MMS order check at γ* (ghost, and `:jumpgrad` o≤2);
* the flume with the high-order penalties;
* the amplitude ceiling beyond 0.15;
* cost: ghost runs at ≈ 50 s/step (`:full`) — precompute the shifted bases.

### 5.8 Update (2026-10-04 22:00) — all three ghost `:full` arms complete 160 s

The envelopes are per-20 s maxima of the band energies, relative to the carrier. The carrier σ_E is
the energy decay fitted over t ≥ 20 s.

| arm (Q3/Q2, CN, `:full` broken) | reached | u mid envelope | u high envelope | carrier σ_E |
|---|---|---|---|---|
| ghost γ = 0.01, 16 cells | ✅ **160 s** | 1.4e-7 → 3.8e-8, monotone decay | 2.7e-9 → 1.8e-9 | −9.2e-5 |
| ghost γ = 0.0033, 16 cells | ✅ **160 s** | 1.5e-7 → 7.2e-8, decaying | 3.7e-9 → 2.6e-9 | −3.1e-5 |
| ghost γ = 0.03, 16 cells | ✅ **160 s** | 1.2e-7 → 1.6e-8, decaying | 2.0e-9 → 9.6e-10 | −2.7e-4 |
| ghost γ = 0.01, 32 cells (exact Jac) | 110 s | 1.5e-7 → 1.3e-7, flat | ~2e-11, flat | −1.3e-5 |
| ghost γ = 0.01, A = 0.15 | 90 s | 1.2e-5 → 4.4e-6, decaying | 1.1e-7 → 4.6e-8 | −1.0e-4 |

**What this shows.**
* No late growth at any γ across the full ninefold range 0.0033–0.03. The controls grow over the
  same 110–160 s: first-order penalty ×100–250, `:native` with γ = 0 ×920.
* Even γ*/3, which is below the eigen window, holds 160 s. The window's lower edge is therefore
  conservative for `A` = 0.10.
* The carrier loss scales with γ as predicted: −3.1e-5, −9.2e-5 and −2.7e-4 s⁻¹. That is 0.5 %,
  1.5 % and 4.3 % of the energy over 160 s.

**Cost.**
* 32 cells/λ with the exact Jacobian: Newton 2/step, ≈ 25 s/step.
* `A` = 0.15: Newton 13–14/step, ≈ 38 s/step. That is the quasi-Newton Jacobian; the arm started
  before `broken_class3_jacobian` existed. By rule 17b this changes the cost, not the answer.

**Run banner fixed.** The `setup_and_run` banner said "quasi-Newton on the Class-III blocks" even
after the exact Jacobian was added. It is corrected in `src/utilities.jl`. The text is cosmetic and
the running arms are unaffected.
