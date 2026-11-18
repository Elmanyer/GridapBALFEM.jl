# MIXED_VMS_PLAN.md — the mixed-formulation solver with residual-based (VMS) stabilisation

*Branch `mixed-vms-solver`, cut from `v2-solver` @ `1be9c39` on 2026-11-16. Author's decision of
that date. This file is the plan and, as steps finish, the record. Status marks: ☐ open, ◐ partly
done, ☑ done (with the gate that says so).*

> **Where it stands (2026-11-18).** The mixed solver (called **v3-mixed**) works: gates 29/29, m5/m6
> 200/200 steps in the short comparison. Without stabilisation it diverges in the closed box (46.7 s
> at 16 cells/λ, 19.9 s at 32). The ASGS form of §3.3 fails (indefinite effective mass) and is set
> aside; OSS (§3.5) fixes the mass and runs models 3–4. On models 5–6 the OSS stabiliser itself
> **destabilises**, because it acts on the 𝓝 part of the residual, which has no partner in its test
> operator (§9.5). Leaving 𝓝 out of the stabilised residual turns it into a working stabiliser in the
> frozen spectrum. The design change is the author's decision (§9.6).

---

## 0. Why this branch exists

The delivered solver stays the v2 one: Class III cellwise plus facet terms, with the mandatory
C⁰-CIP (`v2-solver`). It is cheap, but it is still unstable in long runs, and the time-integrator
study showed a steady (ω ≈ 0) vertical-structure growth that no integrator can damp.

This branch builds a **second, independent solver for the same model**, as a research instrument:

1. **Mixed formulation.** The two Class III objects `∇𝖲` and `∇𝖻` become unknowns. Every unknown
   then carries at most one derivative in the weak form, so the weak form is conforming on `C⁰`
   Lagrange elements. There are **no skeleton terms and no penalty** on this path.
2. **Residual-based stabilisation (VMS).** The only stabiliser. No C⁰-CIP, no ghost penalty.

The question it answers: *with the non-conformity removed by construction, does a consistent,
residual-based stabiliser make the full nonlinear model stable, and at what carrier cost?*

What is already known and is **not** re-tested here: the mixed formulation alone delays the
instability but does not remove it (v1, `archive/HISTORY_V1.md` §5.2e–f), and model 3 (no 𝓝, already
conforming) is unstable under Crank–Nicolson. Both say the cause is not the non-conformity.

## 1. Scope

| in scope | out of scope (refused by the code, not silently ignored) |
|---|---|
| models 5–6 in mixed form, sequential | C⁰-CIP / ghost penalty with `formulation = :mixed` |
| VMS on models 3–4 (3 fields, already conforming) | VMS on the primal 𝓝 solver (its strong residual has third derivatives) |
| VMS on models 5–6 mixed | VMS on linear models 1–2 (stable; nothing to stabilise) |
| exact Jacobians of the new terms | distributed (MPI) mixed runs: **deferred**, §7 step 9 |
| the closed x-periodic box driver | LaTeX chapter (author's, in Overleaf) |

The v2 primal path must stay **bit-identical**: every new term sits behind `formulation === :mixed`
or an attached VMS record.

## 2. The mixed formulation

Notation as in `MODEL.md`: `H = d + η`, `𝖲 = ∇·(H𝗎)` and `𝖻 = 𝗎·∇H` are stacked `VectorValue{Nσ}`
fields.

**Unknowns** (7 fields; the first three keep their indices, so the physics rows need no reindexing):

```
[η, 𝖴x, 𝖴y, 𝖦x, 𝖦y, 𝖥x, 𝖥y]        𝖦ₐ ≈ ∂ₐ𝖲,   𝖥ₐ ≈ ∂ₐ𝖻
```

**Auxiliary rows**, with the integration by parts built in and the boundary term assembled (v1
measured a relative error of 22 in the boundary cells without it):

```
∫ 𝖦ₐ·Ψ + ∫ 𝖲·∂ₐΨ − ∮ 𝖲·Ψ nₐ = 0        ∫ 𝖥ₐ·Ψ + ∫ 𝖻·∂ₐΨ − ∮ 𝖻·Ψ nₐ = 0
```

They carry no time derivative: the system is an index-1 DAE. A non-rest initial state needs
consistent auxiliary fields (`mixed_consistent_ic`).

**Class III blocks** (surface-slope 𝓚 and leading 𝓟 halves of components {1,2,4,5}), through the same
reduced contractions the v2 path uses, fed with the unknowns instead of cellwise Hessians:

```
GU = Σₐ 𝖦ₐ⊗𝖴ₐ     SD = 𝖲⊗(∇·𝗎)     N4 = Σₐ 𝖴ₐ⊗𝖥ₐ
R_K = −∫ H Σₐ ∂ₐH (Wₐ · NK)        NK = W_K ⊙ GU + 𝓚⁽²⁾ ⊙ SD + 𝓚⁽⁴⁾ ⊙ N4
R_P = −∫ H² (NP · ∇·W)             NP = W_P ⊙ GU + 𝓟⁽²⁾ ⊙ SD + 𝓟⁽⁴⁾ ⊙ N4
```

Everything else is the v2 residual unchanged: components {3,6,7,8} direct, the bed-slope half by
exact integration by parts.

**Design decisions**

* **D1. Both arms, always.** v1 could drop `𝖥` by a mask; v2 has no masks, so the layout is always 7
  fields.
* **D2. Auxiliary space.** `VectorValue{Nσ}` Lagrange of order `p_aux`, default `p_u` (v1's
  default), no Dirichlet data. `p_aux = p_u − 1` stays available as a probe.
* **D3. No `p_η ≥ 2` requirement from the weak form.** That requirement came from the cellwise
  Class III terms and the η penalty. The VMS strong residual still wants it (§3.4), so Q3/Q2 remains
  the recommended pair.
* **D4. Not chosen:** taking `𝖲` and `𝖻` themselves as unknowns (2 auxiliary fields instead of 4) and
  differentiating them. It is cheaper, but it is "project, then differentiate", and v1's record is
  on the `𝖦`, `𝖥` form.
* **D5. One residual entry point.** `global_residual`, `jacobian_u` and `jacobian_u_t` add the mixed
  and VMS terms themselves when the problem carries them. Every caller (time loop, residual checker,
  eigen-analysis, AD) then sees the complete operator, and none can assemble half of it.

## 3. The residual-based stabilisation

### 3.1 Strong residuals

`R_c` (scalar) and `R_m,a` (one `VectorValue{Nσ}` per horizontal direction `a`) are the strong form
of the equations the weak form discretises, obtained by undoing each integration by parts:

```
R_c   = η_t + Φ·𝖲 − S_src + μη                                   (+ relaxation, − MMS forcing)
R_m,a = H M u̇ₐ                                                    acceleration
      + g H ∂ₐη Φ                                                  gravity
      + ∂ₐ( H² Σ_c P⁽ᶜ⁾ 𝓛⁽ᶜ⁾(u̇) )                                   leading pressure R_P
      − H ( ∂ₐh Σ_c A⁽ᶜ⁾𝓛⁽ᶜ⁾ + ∂ₐH Σ_c K⁽ᶜ⁾𝓛⁽ᶜ⁾ )                    slope package
      + μ M uₐ                                                     sponge (+ relaxation)
      + H 𝓜 ⊙ (u_b ⊗ ∂_b uₐ) + 𝓖 ⊙ (𝖲 ⊗ uₐ)                        advection
      + 𝓝ₐ                                                         models 5–6 only
```

`𝓝ₐ`, with `∇𝖲 → 𝖦` and `∇𝖻 → 𝖥` everywhere:

```
− (∂ₐh NA₆₈ + ∂ₐH NK₆₈) − H (∂ₐh NA₃ + ∂ₐH NK₃) + ∂ₐ( H NP₆₈ + H² NP₃ )        components 3,6,7,8
− H ∂ₐh ( W_A ⊙ GU + 𝓐⁽²⁾ ⊙ SD + 𝓐⁽⁴⁾ ⊙ N4 )                                   bed-slope half
− H ∂ₐH NK + ∂ₐ( H² NP )                                                       Class III
```

### 3.2 Subscales

Quasi-static algebraic subscales: `η′ = −τ_η R_c`, `𝗎′ₐ = −T (H M)⁻¹ R_m,a`.

The campaign's finding is that what grows is vertical structure, not the depth average. A single
scalar `τ` cannot serve both: the depth-averaged (barotropic) part travels at the wave speed, the
vertical-structure part only with the fluid. So `T` splits along the `M`-orthogonal projector onto
the depth-uniform velocity, `P = e Φᵀ` (`e` = all ones; `M e = Φ`, `Φ·e = 1`):

```
T M⁻¹ = τ_w e eᵀ + τ_a (M⁻¹ − e eᵀ)            both parts symmetric positive semi-definite
τ_w = C_w h_e / (2 c₀)        τ_η = C_η h_e / (2 c₀)        c₀ = √(g d),  h_e = h / p_u
τ_a = C_a [ (2/Δt)²·[dt term] + (2/h_e)² (u₀² + 𝗎·M𝗎) ]^(−1/2)        u₀ = 10⁻³ c₀  (floor)
```

`τ_a |u| ≤ h_e/2` always, so the floor never produces a large term. The `Δt` term is **off by
default**: with it the semi-discrete operator depends on `Δt`, which would confuse the `Δt` ladder
and the eigen-analysis.

### 3.3 The stabilisation terms

Inserting the subscales into the weak form and keeping the first-order transport terms, with the
derivative moved onto the test function:

```
+ ∫ Σₐ (N_b ∂_b Wₐ) · Zₐ          Zₐ = τ_w e (e·R_m,a) + τ_a (M⁻¹ − e eᵀ) R_m,a      SUPG, momentum
+ ∫ τ_w Σₐ ∂ₐq (e·R_m,a)                                                              PSPG, continuity
+ ∫ τ_η (ū·∇q) R_c                                                                    SUPG, continuity
+ ∫ τ_η g H (Φ·∇·W) R_c                                                               LSIC, momentum
```

`N_b[i,j] = Σ_k 𝓜[i,k,j] u_{b,k}` is the layer-coupled advection matrix. `𝓜` is fully symmetric, so
`N_b` is symmetric and the momentum term is `τ ∂Wᵀ N M⁻¹ N ∂U ≥ 0` on the advective part of the
residual: a streamline diffusion that acts on every layer velocity. On the linear shallow-water part
the PSPG and LSIC terms give `τ_w g|∇η|²` and `τ_η g d² (∇·ū)²`.

**Not included, stated:** the adjoint of the vertical-advection tensor `𝓖` (not symmetric; it would
break the positivity above), the zero-order adjoint terms, the cross and Reynolds subscale
stresses, the time derivative of the subscales.

Each of the three coefficients `C_a`, `C_w`, `C_η` switches its terms off at 0.

### 3.4 Limits of consistency, stated

* `R_P` needs cellwise second derivatives of `u̇` and of `η`. At Q2/Q1 the η Hessian is zero in 1-D,
  so those terms are truncated: use Q3/Q2.
* Models 5–6: `∂ₐ(H² NP)` needs `∂ₐ(∇·𝗎)`, a cellwise second derivative of `𝗎`. Everything else in
  `𝓝ₐ` is first order in the mixed unknowns.
* Model 6: `∂ₐ(H² NP₃)` needs the third derivative of the bathymetry. If it cannot be evaluated,
  VMS on model 6 is **refused**, never assembled without the term.

### 3.5 Orthogonal subgrid scales (option `vms_oss`, added 2026-11-17, author's request)

*Why:* the ASGS form above weights the full residual, including its time-derivative part (the
acceleration and the dispersive `−d²B̃∇∇·𝗎̇`). Momentum SUPG then adds an O(C_a), non-symmetric
block to `∂R/∂u̇` and the effective mass turns indefinite once the state moves (§9.3): exact
Newton fails. OSS (Codina 2000, 2002) removes that by construction.

*What:* the subscale is the part of the residual **orthogonal to the FE space**,
`𝗎′ₐ = −T(HM)⁻¹ P⊥(Rˢ_m,a)`, `η′ = −τ_η P⊥(Rˢ_c)`, `P⊥ = I − P_h` (L² projection), where `Rˢ` is the
**static residual**: `R` with every `u̇` term removed (acceleration, the whole `𝓛(u̇)` pressure
package, `η_t`). The four terms of §3.3 keep their shape with `R → Rˢ − P_h(Rˢ)`; the τ's are
unchanged.

*Why removing the u̇ part is the OSS argument, not a truncation:* `P⊥` annihilates whatever lies
in the FE space, so `P⊥(∂ₜu_h) = 0` for a constant mass. Here the time part carries the weights
`H` and the dispersive second derivative, so its projection is not exactly zero but `O(hᵖ)`, which
is the order to which OSS is consistent anyway (OSS is weakly consistent: `P⊥(R(u_h)) ≠ 0`).

*How (monolithic, exact Newton):* the projections are unknowns, appended after every other field:
`Πx, Πy ∈ VectorValue{Nσ}` (velocity order), `δ` scalar (surface order); rows
`∫(Πₐ − Rˢ_m,a)·χ = 0`, `∫(δ − Rˢ_c)ϑ = 0`. No time derivative (index-1 DAE, like 𝖦, 𝖥): a hot
start needs consistent projections (`oss_consistent_ic`; Crank–Nicolson does not damp an
inconsistent algebraic start). Layout: 6 fields on models 3–4, 10 on mixed models 5–6.
Consequences, gated (G6): `∂R/∂u̇` is **bitwise the Galerkin mass**; the residual does not read
`u̇`; the Jacobian of the terms and rows is exact (jets).

*Usage:* `setup_and_run(…; vms=true, vms_oss=true)`, env `BALFEM_VMS_OSS=1`; output names carry
`vms-oss-…`; refused with Dirichlet wave generation (no IC path yet) and with all coefficients 0.

## 4. Exact Jacobians without a second copy of the formulas

The new terms are polynomials in the fields. They are written once, on a small value-and-tangent
pair type (`Jet`, `src/jets.jl`): the value part is the residual integrand, the tangent part its
exact directional derivative. Seeding the tangents with the trial functions of `u` gives `∂R/∂u`;
seeding those of `u̇` gives `∂R/∂u̇`. Residual and Jacobian cannot drift apart, and identically-zero
branches (flat bed, switched-off coefficients) are pruned by a zero sentinel instead of assembled.

Rule 6 is kept: `∇` is only ever applied to raw fields and raw bases.

The v2 hand Jacobians of the three physics rows are unchanged (quasi-Newton in the same places).

## 5. Code layout

| file | change |
|---|---|
| `src/jets.jl` | **new.** `Jet`, `JZero`, the overloaded algebra |
| `src/mixed.jl` | **new** (ported from `v1-solver:src/mixed.jl`). `attach_mixed!`, `mixed_contrib`, `mixed_n_aux`, `make_initial_conditions_mixed`, `mixed_consistent_ic` |
| `src/vms.jl` | **new.** `attach_vms!`, `vms_fields`, `strong_residuals`, `vms_contrib` |
| `src/problem.jl` | `BALFEMProblem` gains `formulation`, `mixed`, `vms`; `build_problem(…; formulation)`; the three entry points add the new terms; the CIP requirement applies to `:primal` only |
| `src/horizontal.jl` | `build_fe_spaces(…; n_aux, p_aux)` |
| `src/utilities.jl` | `setup_and_run(…; formulation, p_aux, vms, vms_c_adv, vms_c_wave, vms_c_eta, vms_dt_term)`, the gates of §1, consistent IC |
| `src/monitor.jl` | diagnostics on a 7-field space |
| `src/configuration.jl`, `src/checkpoint.jl` | the new keys |
| `examples/local_1d/run_periodic_1d.jl` | `BALFEM_FORMULATION`, `BALFEM_P_AUX`, `BALFEM_VMS`, `BALFEM_VMS_CADV/_CWAVE/_CETA`; name tokens `mixed`, `vms…` |
| `test/test_mixed_vms.jl` | **new.** The gates of §6 |

## 6. Gates

| gate | what it shows | pass |
|---|---|---|
| G0 no leakage | the v2 primal residual and Jacobians are unchanged by this branch | regression snapshot 12/12; bitwise residual with everything off |
| G1 auxiliary rows | `𝖦`, `𝖥` from `mixed_consistent_ic` satisfy their rows | relative ≤ 1e-10 |
| G2 mixed ↔ broken | on a smooth state both Class III treatments approximate the same operator | difference falls under refinement |
| G3 weak ↔ strong | `∫ v·R_strong` equals the assembled weak residual up to the inter-element jumps of a `C⁰` field; this checks every sign of §3.1 independently of the weak form's own tests | difference falls under refinement, models 3, 4, 5 (6 if available) |
| G4 Jacobians | `Jet` tangents against a central finite difference of the assembled mixed and VMS terms | relative ≤ 1e-6 |
| G5 knobs live | each `C` changes the residual; all zero ⇒ no contribution (rule 38d) | exact |
| G6 smoke | the box driver runs models 3 and 5 with VMS, Newton converges | completes |

## 7. Steps

1. ☑ Branch, this plan.
2. ☑ `src/jets.jl`.
3. ☑ `BALFEMProblem` fields, `build_problem(…; formulation)`, FE spaces with auxiliary fields. G0.
4. ☑ `src/mixed.jl` and its wiring into the three entry points. G1, G2, G4 (mixed part).
5. ☑ `src/vms.jl`: strong residuals. G3.
6. ☑ `src/vms.jl`: τ, the four terms, wiring. G4 (VMS part), G5.
7. ◐ Driver, diagnostics, configuration card, box example. G6 9/9 with OSS (2026-11-18): models 3–4
   primal + OSS run 200/200 steps at ≈ 3 Newton iterations. Models 5–6 mixed + OSS: §9.5–9.6. Open.
8. ☑ Record: `STATUS.md`, `CLAUDE.md`, `INDEX.md`, `ARCHITECTURE.md`, `MODEL.md` §10, this file (§9).
9. ☐ **Deferred, needs a decision:** distributed mixed runs (the residual is partition-agnostic by
   construction; `mixed_consistent_ic` needs a distributed solve; the Jacobi preconditioner on the
   auxiliary mass blocks is untested).
10. ☐ **Not part of the implementation:** the stability study itself (box runs under
    Crank–Nicolson against the CN twins of `output/1D_stability_campaign/`, acceptance by
    `STABILITY.md` §6, carrier cost by rule 15).

## 8. Risks

* **Cost.** 7 fields instead of 3: `(1 + 6Nσ)` unknowns per node instead of `(1 + 2Nσ)`, and a
  denser Jacobian. Sequential LU only.
* **Operator-tree depth.** The VMS Jacobian is a long sum of products; Gridap compiles one kernel
  per tree. First assembly may take many minutes.
* **Auxiliary pairing.** `𝖦`, `𝖥` in the velocity space is v1's untested pairing; a negative
  stability result stays ambiguous between the stabiliser and this choice (`p_aux` is the probe).
* **τ is a model.** The barotropic/vertical-structure split and the constants `C` are choices; none
  is calibrated. The study, not the implementation, decides them.

---

## 9. Execution record (2026-11-16/17)

### 9.1 What was built

| file | contents |
|---|---|
| `src/jets.jl` | `Jet(value, tangent)` / `JZero`: one expression gives the residual integrand (value) and the exact Jacobian integrand (tangent, seeded by `du` or `dut`); zero tangents are pruned at build time (§4) |
| `src/mixed.jl` | ported from `v1-solver:src/mixed.jl`: `attach_mixed!`, `mixed_contrib` (auxiliary rows + the Class III blocks through `𝖦`, `𝖥`), `make_initial_conditions_mixed`, `mixed_consistent_ic` (projects `𝖦`, `𝖥` from the initial `η`, `𝖴`) |
| `src/vms.jl` | `attach_vms!`, strong residuals (`vms_strong_residuals`, `vms_strong_nlp`), τ (`_tau_a`, `τ_w`, `τ_η`), the four VMS terms (`vms_contrib`), `vms_strong_contrib` (the G3 check), `_bed_third` |
| `src/problem.jl` | `BALFEMProblem` fields `formulation`, `mixed`, `vms`; `build_problem(…; formulation)`; `global_residual`, `jacobian_u`, `jacobian_u_t` add `mixed_contrib` / `vms_contrib` (D5), so the ODE operator, the residual checker and AD see the full operator |
| `src/horizontal.jl` | `build_fe_spaces(…; n_aux, p_aux)`: the four auxiliary `VectorValue{Nσ}` fields |
| `src/broken.jl` | CIP is mandatory only on the primal 𝓝 path; a mixed problem refuses a skeleton penalty |
| `src/utilities.jl`, `monitor.jl`, `checkpoint.jl`, `configuration.jl` | driver knobs, refusals, mixed IC, diagnostics padding, checkpoint keys, `CONFIGURATION.txt` lines |
| `examples/distributed/_dist_common.jl`, `examples/local_1d/run_periodic_1d.jl` | environment knobs; the box driver passes them through |
| `test/test_mixed_vms.jl` | gates G0–G5 |
| `examples/local_1d/compare_mixed_vms.jl` | models 1–6, v2 against the new arms (§9.3) |
| `examples/local_1d/vms_newton_probe.jl` | which VMS term stops Newton (§9.3) |
| `run/local/run_mixed_vms_compare.sh` | the comparison arms in parallel, one process per arm, memory-gated |
| `examples/local_1d/vms_jacobian_gap.jl` | hand vs AD `A = ∂R/∂u`, `B = ∂R/∂u̇` per VMS variant and amplitude; quasi-Newton contraction ρ; dominant gap block (§9.3) |
| `run/local/run_mixed_box_controls.sh` | the unstabilised v3-mixed closed-box controls, n16 and n32 (§9.3) |
| `examples/local_1d/oss_newton_trace.jl` | Newton trace and energy of the first steps of a VMS variant (`BALFEM_PROBE_MODEL`, `_STEPS`, `_VARIANT` = all1, adv_eta, wave1, wave0.1, wave0.3, adv1, eta1, all0.1; `BALFEM_VMS_DT`) (§9.5) |
| `examples/local_1d/oss_cw_spectrum.jl` | frozen spectrum of m5 mixed + OSS per variant and amplitude, with the `_noN` discriminator (𝓝 left out of the stabilised residual) and the character of the fastest mode (§9.5) |

### 9.2 How to use it

`setup_and_run` keywords (and the environment variables read by `run_periodic_1d.jl`):

| keyword | env | values | default |
|---|---|---|---|
| `formulation` | `BALFEM_FORMULATION` | `:primal` \| `:mixed` (needs `nl_pressure=true`; CIP defaults to 0 and γ > 0 is refused) | `:primal` |
| `p_aux` | `BALFEM_P_AUX` | order of `𝖦`, `𝖥` | `p_u` |
| `vms` | `BALFEM_VMS` | `true` on models 3–4 (primal) or 5–6 (mixed); refused on 1–2 and on primal 𝓝 | `false` |
| `vms_c_adv`, `vms_c_wave`, `vms_c_eta` | `BALFEM_VMS_CADV/_CWAVE/_CETA` | `C_a`, `C_w`, `C_η` ≥ 0; all zero ⇒ no VMS record | 1, 1, 1 |
| `vms_dt_term` | `BALFEM_VMS_DT` | add `(2/Δt)²` inside `τ_a` | `false` |

Output names carry `mixed` and/or `vms` (`formulation_tags()`); `CONFIGURATION.txt` states the
formulation, the auxiliary order, the Class III treatment and the VMS constants. Mixed + Dirichlet
wave generation (`wave_bc`) is refused (no consistent auxiliary boundary data). Sequential only.

Cost: the first assembly of the mixed residual compiles for ≈ 8 min and its Jacobian ≈ 4 min; the
gate suite takes 45–60 min cold.

### 9.3 Results

**Gates** (`test/test_mixed_vms.jl`, full run 2026-11-17, `output/mixed_vms/_logs/test_mixed_vms_full.log`):
**29 passed, 0 failed.**

| gate | result |
|---|---|
| G0 | 11/11 (no leakage, every refusal). Regression snapshot 12/12 (`_logs/regression_snapshot.log`): models 1–4 bitwise identical, 𝓝 configurations ≤ 1.6e-15 |
| G1 | auxiliary rows 1.06e-15 |
| G2 | mixed ↔ broken 5.77e-3 → 1.17e-3 (nx 8 → 16, ratio 4.9) |
| G3 | weak ↔ strong, ratio 3.3 in models 3, 4, 5, 6 (fine mesh ≈ 1.8e-2). The test must refine **both** directions (`ny = nx ÷ 3`): refining only `nx` left a mesh-independent 5.8 % from the y-variation |
| G4 | jet Jacobians against central FD ≤ 6e-11 (models 3, 5 mixed, 6 mixed; ∂/∂u and ∂/∂u̇) |
| G5 | each coefficient live (C_a 8.6e-2, C_w 2.65e-1, C_η 1.29e-1); C_w term linear in C_w |

**Comparison box** (`examples/local_1d/compare_mixed_vms.jl`; 1 λ, 16 cells/λ, Q3/Q2, P1LFE-2,
d = 3.5 m, T = 1.6 s, kA = 0.08, CN Δt = 0.04 s, 5 periods; output `output/mixed_vms/comparison/`,
table in `summary.txt`):

* m1, m2, m3, m4 primal: 200/200 steps, E/E₀ = 1.000000 / 1.000003 / 1.000305 / 1.000091, Newton 1 / 1 / 3 / 3.
* m5, m6 mixed and mixed + VMS (C = 1, 1, 1): see below. Final arm table (`summary.txt`): m1–m4
  primal 200/200; m3/m4 primal + VMS FAILED (step 2, hand) and FAILED (steps 16–24, AD); m5/m6 primal
  + CIP and mixed 200/200 (mixed ≈ 2× faster, gauge difference ≈ 0.3 %); m5/m6 mixed + VMS NaN
  (hand) and stopped (AD).
* **Every VMS arm with the hand Jacobian fails at step 2** (C = 1, 1, 1): m3 stalls at 3.6e-5, m4 at
  1.2e-6 (50 iterations); m5 and m6 mixed give NaN (39 / 32 iterations). Step 1 converges, slowly
  (6–9 iterations against 3 without VMS).
* **First reading (superseded below): the quasi-Newton core Jacobian** (`vms_newton_probe.jl`, m3, 3 steps,
  `output/mixed_vms/newton_probe/`):

  | variant | Newton per step | final residual |
  |---|---|---|
  | all C = 1, hand Jacobian | 6, then stall | 3.6e-5 |
  | **all C = 1, exact AD Jacobian** | **3 3 3** | 1e-14 – 7e-11 |
  | C_a = 1 alone | 3 4 3 | ≤ 9e-9 |
  | C_w = 1 alone | 6 6 6 | ≤ 6e-9 |
  | C_η = 1 alone | 3 3 3 | ≤ 7e-11 |
  | C_w = 0.1 / C_η = 0.1 / all 0.1 | 3 3 3 | ≤ 5e-9 |

  The VMS Jacobian itself is exact (G4). The core `∂R/∂u` omits, by design, some η-derivatives
  (`MODEL.md` §8, rule 5); without VMS that costs nothing (3 iterations), but the barotropic VMS terms
  (C_w: momentum SUPG barotropic part + PSPG) change the iteration matrix enough that the omission no
  longer contracts. With the exact Jacobian Newton converges quadratically **for the first steps
  only** (see below: it fails at steps 16–24, so the Jacobian is not the root cause). ~~Decision
  needed: AD or complete the core Jacobian~~ — superseded 2026-11-17 by the measurements below.
* The VMS arms are rerun with `use_ad = true` (arms `*_vms_ad`). **⛔ Superseded: the exact Jacobian
  fails too.** m3 primal + VMS (AD) steps cleanly to step 15 (3–4 Newton, 1e-13), then Newton blows up
  (‖r‖ 3e198) between steps 16 and 19; m4 likewise after step 20 (1e206). The same models without VMS
  run 200/200. AD cost: 57–100 s/step against 4 s/step by hand (15–25×, not 3–4×). m5/m6 mixed + VMS
  (AD) did not finish step 1 in 5 h 38 min and were stopped by hand (author, 2026-11-17).
* **What fails, measured** (`examples/local_1d/vms_jacobian_gap.jl`, model 3 box state, a = amplitude /
  run amplitude; `output/mixed_vms/jacobian_gap/report.txt`, `spectrum.log`):
  * **No wiring bug.** The VMS part of `A_hand − A_AD` is ≤ 7e-17 for every variant and amplitude;
    `B_hand = B_AD` exactly. The only gap is the core's documented one (rows 𝖴x × columns η, the
    η-derivatives rule 5 leaves out): ‖ΔA‖/‖A‖ = 0.32 at a = 1.
  * **Quasi-Newton contraction ρ = spectral radius of I − J_hand⁻¹ J_AD** (CN, Δt = 0.04):

    | a | no VMS | C_a = 1 | C_w = 1 | C_η = 1 | all 1 | all 0.1 |
    |---|---|---|---|---|---|---|
    | 0.25 | 0.0008 | 0.0034 | 0.48 | 0.0008 | 0.17 | 0.05 |
    | 1 | 0.0033 | 0.020 | **0.91** | 0.0033 | 0.58 | 0.21 |
    | 2 | 0.0066 | 0.034 | **0.98** | 0.0066 | 0.86 | 0.39 |

    The **C_w** terms (barotropic SUPG + PSPG) make the iteration matrix amplify the core gap: the
    hand-Jacobian stall at step 2.
  * **Frozen linearised spectrum** `B u̇ = −A u` (exact matrices, a = 1): without VMS B (= ∂R/∂u̇) is
    positive definite, max |λ| 8.18 s⁻¹. With **C_a = 1** (momentum SUPG, advective part) **B acquires
    8 negative eigenvalues** (all 1: 6), min |eig B| falls 2.7e-4 → 4.9e-5, max |λ| rises to 21–24 s⁻¹.
    At rest (a = 0) every variant leaves B and the spectrum unchanged. C_w and C_η keep B definite.
    Reading: the SUPG weight `τ_a (N·∇)W` also tests the **time-derivative part** of the momentum
    residual, which contains the dispersive `−d² B̃ ∇∇·𝗎̇` (second derivatives of 𝗎̇). Its perturbation
    of the mass is O(C_a), not small, and it makes the effective mass indefinite once the state
    carries velocity. An indefinite B lets the CN matrix A + B/(θΔt) pass through singularity as the
    state evolves: a step with no nearby solution, which is how the exact-Newton runs fail.
  * ⚠ Frozen growth rates are not predictive (`STABILITY.md` §2); the definiteness of B is a
    structural fact, not a rate.
* Models 5–6 without VMS: m5 mixed 200/200 (E/E₀ 1.000081, 3 Newton), m6 mixed 200/200 (1.000042, 3);
  v2 reference m5 primal + CIP 200/200 (1.000075). Gauge differences: `summary.txt`.
* **Unstabilised controls: v3-mixed without VMS is unstable, and refinement makes it worse**
  (`run/local/run_mixed_box_controls.sh`, the chapter-8 campaign environment copied verbatim with
  `BALFEM_FORMULATION=mixed BALFEM_VMS=0`: m5 flat bed, Q3/Q2, A = 0.10, CN Δt = 0.04, seed 1e-8,
  1 λ, 200 periods; output `output/mixed_vms/box_controls/`, band analysis
  `_logs/periodic_growth_n{16,32}.log`; 2026-11-17):

  | cells/λ | E/E₀ > 1.01 | Newton NaN (end) | σ_E mid band | σ_E HIGH band | σ_E sub-element | carrier σ_E |
  |---|---|---|---|---|---|---|
  | 16 | 38.0 s | **46.7 s** (step ≈ 1167) | +0.37 s⁻¹ | +0.26 s⁻¹ | +0.28 s⁻¹ | −0.0007 s⁻¹ |
  | 32 | 16.2 s | **19.9 s** (step ≈ 497) | +0.71 s⁻¹ | +0.95 s⁻¹ | +0.95 s⁻¹ | +0.001 s⁻¹ |

  (σ_E = band energy growth rate of η, fit over the run; the u rates agree within 10 %.) The growth is
  at the grid scale and its rate rises ×3.6 in the HIGH band from n16 to n32: a discretisation
  instability (rule 38b), not under-resolution. The carrier is untouched until the blow-up.
  Comparison, same environment: v2 primal **without** CIP 81 s / 27 s (n16 / n32, `STABILITY.md` §3b);
  v2 primal **with** CIP clean at n16, burst at 245 s at n32; v1 mixed (component 4 off) Q3/Q2 n16
  26 s. So conformity alone does not stabilise: v3-mixed fails earlier than the v2 primal without
  CIP, and it needs a stabiliser as much as v2 does. VMS (§9.4) is that stabiliser and is not yet
  usable.

### 9.4 Open decisions (author, 2026-11-17)

The VMS failure has two separate causes; neither is a wiring bug (§9.3).

1. **The mass problem (exact Newton fails; the blocker).** Momentum SUPG (`C_a`) weights the
   time-derivative part of the momentum residual, which carries the dispersive `−d² B̃ ∇∇·𝗎̇`. That
   adds an O(C_a), non-symmetric, second-derivative term to `∂R/∂u̇`, so the effective mass becomes
   indefinite once the state moves. Candidate fixes:
   * **(a, recommended first, a few lines in `vms.jl`)** leave the dispersive `𝗎̇` term out of the
     residual weighted by SUPG/PSPG; the mass stays the Galerkin mass. Cost: an explicit loss of
     consistency in the dispersive term only (to state in §3.4). Check: the `spectrum.log`
     measurement (B definite at a = 1, 2) before any run;
   * **(b)** orthogonal subgrid scales (OSS, Codina): the subscale is the residual's component
     orthogonal to the FE space, so `∂ₜ` drops out by construction and consistency is kept; needs an
     L² projection of the residual (more fields or a projection solve per Newton iteration);
   * **(c)** the quasi-static shortcut: drop the whole `∂ₜ` part from the weighted residual (simpler
     than (b), less consistent than (a)).
2. **The stall problem (hand-Jacobian Newton stalls at step 2).** The `C_w` terms push the
   quasi-Newton contraction ρ to 0.91–0.98 at run amplitude. Either complete the core η-derivatives
   of the momentum rows (rows 𝖴 × columns η; a rule-5 decision for the author) or lower `C_w`
   (ρ ≈ 0.2 at `C_w` = 0.1). AD is not an option for production: 15–25× per step (§9.3).

No change to `vms.jl` or the core Jacobian has been made; both wait for the author.

**Decided 2026-11-18 (author):** (1) → option (b), OSS, implemented as §3.5 (`vms_oss=true`); the ASGS
form stays available but is set aside. (2) → the η-derivatives of the pressure terms in `jacobian_u`
are completed: models 3–4 hand ≡ AD to round-off (ρ = 0), `test_jacobians_ad` 14/14; rule 5 now
covers only the 𝓝 blocks. The m5 mixed + OSS failure is the same with the old and the completed
Jacobian (rule 17b).

### 9.5 Models 5–6 mixed + OSS: what fails and why (2026-11-18)

Box of the comparison (m5 flat bed, 1 λ, 16 cells/λ, Q3/Q2, P1LFE-2, kA = 0.08, CN Δt = 0.04).
Time probes `oss_newton_trace.jl` (≈ 60 s per step, ≈ 12 min compilation); logs
`output/mixed_vms/_logs/oss_trace_m5_*.log`, runs `output/mixed_vms/newton_probe/oss_m5_*/`.

**C_a: τ_a unbounded.** Without the Δt term τ_a reaches ≈ 7 s ≫ Δt at velocity nodes and the run
blows up at step 2 (adv1, all1). With `vms_dt_term` (τ_a ≤ Δt/2) C_a and C_a + C_η run 25 / 40 steps
at 3 Newton iterations. ⚠ **But that run is inert:** its energy matches the unstabilised m5 mixed run
to ΔE/E₀ ≤ 2.2e-7 over 40 steps (E/E₀ − 1 = 2.14489e-3 vs 2.14467e-3 at step 40). It is not evidence
of a working stabiliser.

**C_w: a growth riding the crest.** C_w = 1 alone (Δt term on): flat to step 16, E/E₀ = 1.015 at step
24, then ×1.85 per step in energy (≈ 15 s⁻¹ in energy, ≈ 8 s⁻¹ in amplitude); Newton fails at step 35.
max|η| follows the crest. All three coefficients follow the C_w-only run step for step. **C_w = 0.1**
(100 steps): the energy excess over the unstabilised run grows from 1.3e-4 (step 55) to 3.8e-3 (step
80); E/E₀ = 1.044 and max|η| = 7.7e-2 (+52 %) at step 100, Newton up to 4 iterations. The same
instability, ≈ 2.8× later. A smaller C_w delays it and does not remove it.

**Mechanism (hypothesis H1).** Inserting the subscales gives `∫ L*(v)·τ P⊥(Rˢ)`. The test operator
L*(v) holds the transport parts only (advection `N_b∂_bW`, PSPG `Φ∂q`, LSIC `∇·W`). When L*(u) matches
the part of Rˢ it multiplies, the term is `Rˢ·τ·P⊥Rˢ ≥ 0` in the energy pairing (test `(gη, 𝗎)`): a
dissipation. That holds for the continuity pair (`ū·∇η` + `HΦ·∇·𝗎` = `∇·(Hū)`) and for the
advection/PSPG pair. **But Rˢ also holds 𝓝**, whose strong form `∂ₐ(H² NP)` carries `∂𝖦 ~ ∂³𝗎`, and L*
has no partner for it. The product `τ ∇q·P⊥(e·𝓝)` (and `τ (N∂W)·P⊥(𝓝)`) has no sign. It is
proportional to the velocity and, at the grid scale, larger than the positive part by roughly
(U/c₀)(kH)², i.e. O(10²) here.

**Frozen spectrum** (`oss_cw_spectrum.jl`; report `output/mixed_vms/oss_spectrum/report.txt`). At the
hot-start eigenmode ×a: A = ∂R/∂u (hand core + exact jet VMS part), B = ∂R/∂u̇ (Galerkin mass),
algebraic fields (𝖦, 𝖥, Π, δ) eliminated by Schur complement; max Re λ of B u̇ = −A u (s⁻¹). `_noN`
evaluates the same stabilisation with 𝓝 left out of Rˢ (a copy of the problem with
`nl_pressure = false`, used only inside `vms_contrib` and `oss_consistent_ic`; nothing in `src/` changed):

| variant | a = 0 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| no stabilisation | 6e-13 | 0.067 | 0.30 | 0.77 |
| C_w = 1 | 1e-13 | 3.63 | **9.01** | 79.99 |
| C_w = 0.1 | 6e-13 | 1.15 | 2.70 | 5.72 |
| C_a + C_η (τ_a ≤ Δt/2) | 2e-13 | 0.064 | 0.37 | 1.65 |
| **C_w = 1, 𝓝 out** | 1e-13 | 0.049 | **0.099** | 0.19 |
| **all three, 𝓝 out** | 2e-13 | 0.032 | **0.069** | 0.15 |

* H1 holds: removing only the 𝓝 part of the stabilised residual takes C_w from 9.0 to 0.099 s⁻¹ at
  run amplitude. The frozen 9.0 s⁻¹ reproduces the ≈ 8 s⁻¹ of the time run.
* With 𝓝 in, C_a + C_η also destabilise at larger amplitude (×2 at a = 2, fastest mode at the grid
  scale, k h_e/π ≈ 0.85). The growth rate scales roughly as √C_w (0.1× C_w → 0.3× rate).
* With 𝓝 out, all three terms give 2–5× less growth than no stabilisation at every amplitude.
* The unstabilised fastest mode is η-dominated and oscillatory; the one left with 𝓝 out is η-free,
  non-oscillating (Im λ ≈ 0) and lives in the velocity (η share 0; vertical-structure share of its
  kinetic energy 0.27–0.85): the same kind
  as the v2 steady growth (`STABILITY.md`).
* ⚠ The core ∂R/∂u is quasi-Newton in the 𝓝 blocks (rule 5): the comparison between variants (same
  core) is the measurement, the absolute rates are indicative. Frozen spectra are not predictive of
  long runs (`STABILITY.md` §2).

### 9.6 Open decisions (author, 2026-11-18)

1. **Take 𝓝 out of the stabilised residual.** Consistency is kept at the order OSS already has: for a
   smooth solution P⊥ of any smooth term is O(hᵖ). Two forms:
   * **(a) minimal:** an option that evaluates Rˢ without 𝓝 in `vms_contrib` and `oss_consistent_ic`
     (what `_noN` measured);
   * **(b) term-by-term ("split") OSS** (Codina 2008): `Σ_k τ_k ∫ P⊥(L_k v)·P⊥(L_k u)` over the
     operators to stabilise (layer advection `N_b∂_b𝗎` with τ_a; the barotropic pair `Φ∂η`,
     `∇·(Hū)`; possibly a term for the highest derivative of 𝓝 with its own τ). Dissipative by
     construction, whatever the other terms of the model.
2. **Then the first stability test:** closed-box m5 at 16 and 32 cells/λ (chapter-8 environment)
   against the unstabilised controls (46.7 s / 19.9 s), with the carrier cost reported (rule 15).
3. **τ_a.** Whether the Δt term becomes the default, and whether τ_a ≤ Δt/2 is too weak for the
   non-oscillating vertical-structure mode that remains.

Defaults of `vms.jl` stay unchanged until these decisions.
