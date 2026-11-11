# STABILITY.md — the nonlinear grid-scale instability and its stabilisation: the condensed record

*Condensed 2026-10-09 from the v1 campaign records (now in [`archive/`](archive/README.md)) and the v2
box campaign. LaTeX: chapter 8 (`SolverValidation/StabilityAnalysis.tex`, the analysis and the
figures), chapter 9 (`Stabilisation.tex`), appendix C (jump penalties). Status and next step:
[`STATUS.md`](STATUS.md).*

Labels: **(v1)** = measured on the frozen v1 solver, where `:full` means Class III by the
mixed/projected/broken treatment with component 4 omitted (`C3_MASK=gs`) and `:native` means
`𝓝 = {3,6,7,8}` only. **(v2)** = all eight components (`𝓝` on) or none (`𝓝` off).

---

## 1. The phenomenon and how it was isolated

* Fully nonlinear runs grow an element-scale mode. **Refining the mesh advances the onset**, and so
  does raising the order. That is a discretisation-stability signature, not under-resolution
  (rule 38b).
* **Two earlier "instabilities" were artefacts, now closed:**
  * equal-order `Q2/Q2` pairings: an inf-sup checkerboard at `λ ≈ 2dx`, cured by Taylor–Hood (rule 12b);
  * `SDIRK_2_2` masking: Gridap's `DIRK22(1,0,1)` damps the fastest modes at 2–10 s⁻¹ and hid the
    growth. **The diagnosis (what grows, how fast) is made under Crank–Nicolson** (rule 15, rewritten
    2026-11-11: dissipation is a legitimate remedy, named and measured, §6b).
* **It is interior.** The flume always failed at the inflow, so the **closed x-periodic box** was built:
  * one model wavelength, `x_periodic=true`, one cell across with walls;
  * no inflow, relaxation, sponge or source;
  * the discrete linear eigenmode plus a deterministic 1e-8 seed on every harmonic;
  * sub-cell VTK, and band energies by exact DFT (`periodic_growth.jl`).

  The continuum null is known: nothing grows at `κa = 0.16` in one wavelength. The box diverges by
  itself, and the flume boundary only makes the failure come sooner (rule 38i).
* **Ruled out by measurement (v1):**
  * the Jacobian: hand and AD crash at the same time (rule 17b);
  * the operator: it collapses onto Yang & Liu to round-off;
  * quadrature, CFL / `dt`, sponge, domain length, the relaxation zone, Benjamin–Feir, and the
    integrator order (RK4 = SDIRK on the flume);
  * the auxiliary-field order in the mixed formulation;
  * the `C⁰` broken-`∂²` hypothesis: the mixed formulation, which cannot form `∂²`, kept the same
    signature (rule 39b).

## 2. Mechanism (chapter 8)

* **Class III carries the horizontal advection of the vertical acceleration.** The frozen spectrum
  gains a Doppler branch `ω ≈ kU + ω∞`, with spectral radius ∝ `A/h_e`. Without it, the radius
  saturates at `ω∞ = 7.86 rad/s`.
* **The carrier pumps energy up the wavenumber ladder** at a rate ∝ `kU`.
* **Nothing removes that energy.** `C⁰` elements represent it out to 3.5–4.4/Δx, beyond the node
  Nyquist, with no null and no dissipation. Yang & Liu's five-point finite difference is capped at
  1.37/Δx with a null at `2Δx`, **and they filter anyway** (nine-point Shapiro, every 20–200 steps, on
  the Ohyama bar).
* **The model without `𝓝` has the same instability, weaker.** It appears only once the effective
  resolution (cells × order) is high enough. Its mechanism is not identified; the advection block is
  the candidate. Its growing content is a vertical-structure velocity mode with almost no depth average,
  the same signature as the 𝓝 + CIP growth (§3b, 2026-11-07).
* ⚠ Frozen growth rates are **not** predictive (they rank `𝓝`-off above `𝓝`-on). Only the spectral
  radius and the box's Floquet-type rates are.
* The lower frequency branch of the frozen `𝓝`-on spectrum is unexplained.

## 3. Unstabilised evidence

**v2, closed box, CN, A = 0.10, `dt` = 0.04, 100 periods (2026-10-06…09).** Onset of divergence, or
the late growth rate if the run was still bounded at 160 s:

| | 8 cells/λ | 16 | 32 | 64 |
|---|---|---|---|---|
| Q2/Q1, `𝓝` on | bounded (≤ +0.02) | bounded | growing +0.067 | **81 s** |
| Q2/Q1, `𝓝` off | bounded | +0.03 | +0.024 | growing +0.08 |
| Q3/Q2, `𝓝` on | bounded (high band +0.002) | **82 s** | **27 s** | **9 s** |
| Q3/Q2, `𝓝` off | bounded | +0.03 | growing +0.07 | **44.6 s** |

* A = 0.15, Q2/Q1, 16 cells/λ: **57 s** with `𝓝`, bounded without.
* Runs: `output/local_1d/periodic/v2_p{21,32}_nlp{0,1}_n{8,16,32,64}[_A15]_cn/`. Logs: `…/periodic/_logs/`.
  Launchers: `run/local/run_1dper_v2_unstab{,_extra,_q32}.sh`. Each run directory holds `band_energy.csv`.
* The Q3/Q2 `𝓝`-on ladder (bounded → 82 → 27 → 9 s) is the Q2/Q1 ladder shifted one refinement
  earlier (rule 2d).
* v2 fails 1.4–3× later than v1. Restoring component 4 did not create the instability, and did not
  remove it.

**v1 for comparison (mixed `:full`, CN, Q2/Q1):** onset 157 / 104 / 90 / 31 s at 8 / 16 / 32 / 64
cells/λ. Q3/Q2 at 16 cells: 26 s. Projected treatment: 8.8 s. `:native` diverged at Q2/Q1 64 cells
(146 s) and grew at Q3/Q2 16 cells (+0.10 s⁻¹).

## 3b. The chapter-7 solver (𝓝 + mandatory first-order CIP γ = 0.3): the 1-D campaign (2026-10-31 → 11-07)

Closed box, Q3/Q2, CN, A = 0.10, dt = 0.04, **200 periods (320 s)**. Model 6 has a periodic cosine
bed (slope 1:10, one undulation per λ). Full record: `output/1D_stability_campaign/README.md`
("FINAL RESULTS").

| cells/λ | m3 (no 𝓝, no CIP) | m5 (𝓝 + CIP, flat) | m6 (𝓝 + CIP, cosine bed) | m5 *without* CIP (2026-10-07) |
|---|---|---|---|---|
| 64 | 44.6 s (reference, §3) | **clean** | **clean** | 9 s |
| 32 | **246 s** (mid band, σ_E 0.11–0.13 /s from t ≈ 20 s) | burst from 245 s, saturates against the CIP; carrier −12 % | onset ~290–300 s (cut off at 320 s) | 27 s |
| 16 | growing, σ_E ≈ 0.02 /s, still small at 320 s | **clean** | **clean** | 81 s |
| 8 | grid-locked u growth, then E/E0 → 1.11 from 300 s | long-wave u growth from ~140 s; carrier −49 % | the same from ~180–200 s; carrier −38 % | bounded to 160 s (end of run) |

"Clean" means every band flat relative to the carrier for 320 s, carrier −0.4 %, CIP dissipation rate
1e-9 to 4e-6 /s.

* **The CIP removes the fast 𝓝 instability, but not every growth.** The n32 onset is late and the
  n8 one is slow. The dx dependence is **not monotone** (n16 and n64 clean, n32 not). Cause unknown.
* **New: the growing content is a vertical-structure velocity mode, in every case including model 3.**
  The layer velocities u₁, u₂, u₃ carry comparable energy while the depth average ū = Σ D_j u_j
  (the only thing continuity sees) stays far below the physical surface-wave share. So η barely
  moves at onset.
  * This is the non-hydrostatic branch (ω∞, §2).
  * The CIP, which penalises ⟦∂ₙ𝖴⟧, cannot reach its long-wave form (n8: n = 2–4).
  * A stabiliser must act on every u_j, not on η or ū.
  * Evidence: `output/1D_stability_campaign/_logs/vertical_structure_2026-11-07.log`.
* ⚠ Two axes separate m3 from m5 (𝓝 and the CIP, rule 14c). The prepared no-CIP controls
  (`run/local/run_1d_nocip_controls.sh`) and an m3 + CIP arm have not been run.

## 4. What was tried, and what each did (v1 unless stated otherwise)

All arms below are the closed box at Q3/Q2, 16 cells/λ, CN, A = 0.10, 160 s, unless stated otherwise.

| treatment | result |
|---|---|
| broken Class III, no penalty | exact (layer identity holds to 7e-12) but **not a stabiliser**: box 44 s, flume 38 s |
| **first-order C⁰-IP** `⟦∂ₙu⟧⟦∂ₙv⟧`, (γ_u, γ_η) = (0.3, 0.3): **the standard method** | removes the grid-scale growth, but **late mid-band growth from ≈ 110 s at +0.13 s⁻¹** (λ ≈ 4–10 nodes, ×100–250 by 160 s). It **delays** and does not cure. The Q3/Q2 flume with the same penalty held η_max at 0.110 to 158 s (unpenalised: 38 s), but that was read from max η only, with no band analysis |
| first order, velocity only | delays only (Q2/Q1: ≈ 50 → 80 → 110 s for γ_u = 0, 0.03, 0.3); no effect on `:native` |
| first order, surface on Q1 (Q2/Q1) | **locks the carrier** (σ_E = −0.021 s⁻¹); Q2/Q1 has no γ window at any order |
| **hp jump penalty, orders ≤ 2** (`:jumpgrad`, `cip_order=2`), γ = 2e-3 | ✅ 160 s, every band decaying, carrier −1.5e-4 s⁻¹ (`:full` and `:native`) |
| **ghost-volume penalty** (`:ghostvolume`, every order 0…p), γ ∈ [0.0033, 0.03] | ✅ 160 s, decaying; carrier loss 0.5 / 1.5 / 4.3 % of the energy. 32 cells reached 110 s (flat); A = 0.15 reached 90 s (decaying); ≈ 25–50 s/step |

* **The linear eigen-analysis predicts the carrier damping to ±8 %**, so it is a reliable γ design
  tool (`examples/local_1d/cip_eigen_analysis.jl`).
* **Why the first-order penalty escapes: the kernel.** A jump penalty of orders `1…m` vanishes on
  `C^m` piecewise polynomials. Wave modes retreat into that kernel and stay undamped, with damping
  ∝ 1/γ. A penalty is a **constraint**, not a band-limited damper.
* **Consistency terms are not stabilisers.** The broken Class-III layer completes the weak form and
  has no sign in the energy. The penalties are separate, sign-definite additions. `:jumpgrad` order 1
  is only the penalty half of a classical C⁰-IP.

## 5. Why the penalties are not the way forward (2026-10-06…09)

1. **The order needed grows with the element order, and Gridap stops at 2.** The linear
   eigen-analysis (v2 operator, 16 cells/λ, `examples/local_1d/cip_eigen_order_rule.jl`, output in
   `output/local_1d/cip_eigen_order_rule/yfilter/`) gives the damping of the weakest x-varying
   mid-band mode, in s⁻¹:

   | case | γ = 1e-3 | 2e-3 | 3e-3 | 1e-2 | reading |
   |---|---|---|---|---|---|
   | Q3/Q2, order ≤ 2, both fields | 0.37 | 0.65 | 0.91 | 2.8 | window |
   | Q3/Q2, order 1, or velocity only | 3e-6 | 6e-6 | 9e-6 | 2e-5 | escape |
   | **Q4/Q3, order ≤ 2, both fields** | 5e-7 | 2.5e-7 | 1.7e-7 | 5e-8 | **escape, no window** |
   | Q4/Q3, ghost (order p) | 6e-3 | 8e-3 | 0.011 | 0.034 | no escape, weak; needs a larger γ |

   **The order is set by the element (`p − 1`), not by the order of the weak form.** `:jumpgrad`
   therefore cannot serve Q4 and above.
   ⚠ The v1 window values in appendix C were read through a centre-line mode classifier that cannot
   see transverse modes. The corrected window is wider and extends to γ = 0.01.
2. **The ghost penalty works, but it does not generalise yet.** It is restricted to uniform Cartesian
   meshes, sequential runs, AD-incompatible (`_shift`, OPEN_ISSUES §2), costly (≈ 50 s/step), and
   needs a carrier-locking γ window calibrated per pairing.
3. **No penalty has been shown beyond the v1 1-D flat box and a short flume.** The variable-bed
   lee-shoulder mode (OPEN_ISSUES §3), 2-D, MPI and v2 (all eight components) are all untested.
4. **The v2 ladder shows growth in both models**, so whatever is adopted must work for `𝓝` off as well.

**Decision (author, 2026-10-09): pursue filtering, as Yang & Liu do.** The options are in
[`STATUS.md`](STATUS.md) §4. The penalties stay in the code as the comparison arm.

## 6. Acceptance criteria for any stabiliser (v1; adapted 2026-11-11 for integrator dissipation)

Box first, then flume. A weak-form stabiliser is measured under Crank–Nicolson. A dissipative time
integrator is measured with itself, paired with the CN run of the same case (rule 15):
1. mid-, high- and sub-element-band σ_E ≤ 0 on the discriminating cells:
   * Q3/Q2 at 16 and 32 cells/λ;
   * Q2/Q1 at 32 and 64 cells/λ;
   * A = 0.10 and 0.15;
   * `𝓝` on and off;
2. the carrier unattenuated over 100 periods. The unstabilised CN trace has |σ_E| < 4e-5 s⁻¹;
   state the loss as a percentage of the energy. For an integrator this is the deciding cost; the
   accepted budget is the author's to set (§6b gives each candidate's figure);
3. the MMS orders of [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) preserved;
4. the amplitude ceiling it buys measured and stated, since the transfer rate grows with A.

Then: the flume factorial, the bar (where Yang & Liu filter and limit the slope), and a 100-period
CN box regression gate in the suite, which today has no long-duration stability gate.

## 6b. Next: dissipation through the time integrator, before new weak-form terms (author, 2026-11-11)

**The decision.** CN has shown that the discrete nonlinear operator is unstable (every box run so far).
Before adding stabilisation terms (filtering F1–F4, SUPG, VMS, residual-based) to an already complex
weak form, measure what numerical dissipation through the time integration buys, and at what cost.

**Why first.**
* No code change: every Gridap tableau runs through `BALFEM_SOLVER=sdirk BALFEM_TABLEAU=<name>`, and
  `BALFEM_DT` sets the second knob. The weak form, its verification (MMS, AD, parity) and the
  residual stay untouched.
* It is the reference model's own practice: Yang & Liu integrate with RK4.
* It measures how much damping, at which frequencies, the CN-unstable modes need. A later weak-form
  stabiliser has to deliver the same, so the numbers carry over.

**What it can and cannot do (predictions to test).**
* An integrator damps by **frequency**. Per step it removes `≈ c y^q` of a mode's amplitude,
  `y = ωΔt`. The **selectivity** between a growing mode at ω_g and the carrier at ω₀ is `(ω_g/ω₀)^q`,
  fixed by the scheme; `dt` moves only the overall level.
* A steady member (ω = 0) cannot be damped by any integrator (`R(0) = 1`). The m3 decay triad (steady
  shear + 2ω₀ partner) is reachable only through the partner: SDIRK_2_2 damps 2ω₀ 16× more than the
  carrier, RK4 64× more, but RK4 at Δt = 0.04 only at ≈ 1.6e-4 s⁻¹ against a CN growth of ≈ 0.075 s⁻¹.
* Growth at frequencies near ω₀ (the broadband n64 band at T = 1.4/1.5 s, CH8_CAMPAIGN_TODO 2026-11-09
  interim) has selectivity ≈ 1: no integrator can remove it without removing the carrier. That case
  is the expected limit of the approach.
* Explicit schemes amplify every mode faster than `y = √3` (RK3) or `2√2` (RK4): the fastest discrete
  frequency × Δt must stay below that.

**Candidates** (Gridap fork, read from the tableaus 2026-11-11; amplitude loss per step at Δt = 0.04,
T = 1.6 s, and the carrier amplitude lost over 100 periods):

| tableau | stages / order | L-stable | carrier y = 0.157 | 2ω₀ y = 0.314 | y = 1 | carrier, 100 periods |
|---|---|---|---|---|---|---|
| `SDIRK_2_2` (default) | 2 / 2 | no (R(∞) = −½) | 2.2e-4 | 3.0e-3 | 9.9e-2 | 58 % |
| `SDIRK_4_3` | 4 / 3 | yes | 1.3e-5 | 2.0e-4 | 1.7e-2 | 5 % |
| `SDIRK_Crouzeix_3_4` | 3 / 4 | no (0.63) | 6.2e-6 | 3.1e-4 | 4.7e-2 | 2.4 % |
| `DIRK_TRBDF_3_2` (TR-BDF2) | 3 / 2 | yes | 2.2e-6 | 3.5e-5 | 3.1e-3 | 0.9 % |
| `EXRK_SSP_3_3` (RK3) | 3 / 3, explicit | — | 2.5e-5 | 3.9e-4 | 2.8e-2 | 10 % |
| `EXRK_RungeKutta_4_4` (RK4) | 4 / 4, explicit | — | 1.0e-7 | 6.6e-6 | 6.1e-3 | 0.04 % |
| CN / `SDIRK_QinZhang_2_2` | — | — | 0 | 0 | 0 | 0 |

⚠ Excluded: Gridap's `SDIRK_3_3` (first order in fact) and `SDIRK_LobattoIIIC_2_2` (explicit Heun).

**Plan.**
0. *No runs.* Measure the frequencies of the growing modes in every CN-unstable run (temporal
   spectrum: the m3 triad, the m5/m6 n32 burst, the n8 long-wave growth, the n64 broadband band), and
   the fastest discrete frequency (`stability_eig.jl`) for the explicit limits. Tabulate each
   candidate's damping at those frequencies against its carrier cost.
1. *m3 n32 box* (cheapest CN-unstable case; CN diverges at 246 s): every candidate at Δt = 0.04, plus
   a Δt ladder (0.02, 0.08) on the most promising. `SDIRK_2_2` is running (2026-11-10): the CN-unstable
   modes n = 15–16 stay at seed level to 254 s, at a carrier cost of 93 % of the energy.
2. *The hard, cheap case:* m3 T = 1.4/1.5 s n64 (CN diverges at 8.7/12.9 s), the predicted limit.
3. *The full model:* m5/m6 n32 (CN burst ≈ 245 s) with the best one or two; then a steeper carrier
   (k₀A = 0.25, `SDIRK_2_2` launched 2026-11-11).
4. Judge by §6 (σ ≤ 0 over the target duration; carrier loss within the author's budget; MMS orders
   of the chosen scheme; the amplitude ceiling). Only if no integrator passes, or the accuracy cost is
   too high, move to weak-form terms (STATUS §4).

## 7. Tools

| tool | use |
|---|---|
| `examples/local_1d/run_periodic_1d.jl` | the closed box (`x_periodic=true`); env `BALFEM_*` |
| `postprocessing/examples/periodic_growth.jl <run> t0 t1` | band energies and growth rates; `band_energy.csv`. Rank by gain, not amplitude (rule 38f) |
| `examples/local_1d/cip_eigen_analysis.jl`, `cip_eigen_order_rule.jl` | linear damping budget of a penalty; transverse modes filtered by y-structure |
| `examples/local_1d/stability_eig.jl` | frozen-state spectrum (coloured-FD `J*`, skeleton-aware sparsity) |
| `run/local/run_1dper_v2_unstab*.sh` | v2 box launchers (MAXP, `OPENBLAS_NUM_THREADS=2`, logs in `output/…/_logs`) |

Cost: a 100-period box run takes 11–35 h of wall time (Q3/Q2 8 cells: 35 h at 31 s/step with
eight runs sharing the machine) and 2.5–3 GB of RAM.
