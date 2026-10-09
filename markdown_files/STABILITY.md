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
    growth. **Stability claims are made under Crank–Nicolson only** (rule 15).
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
  the candidate.
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

## 6. Acceptance criteria for any stabiliser (unchanged since v1)

All measured under Crank–Nicolson, box first, then flume:
1. mid-, high- and sub-element-band σ_E ≤ 0 on the discriminating cells:
   * Q3/Q2 at 16 and 32 cells/λ;
   * Q2/Q1 at 32 and 64 cells/λ;
   * A = 0.10 and 0.15;
   * `𝓝` on and off;
2. the carrier unattenuated over 100 periods. The unstabilised CN trace has |σ_E| < 4e-5 s⁻¹;
   state the loss as a percentage of the energy;
3. the MMS orders of [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) preserved;
4. the amplitude ceiling it buys measured and stated, since the transfer rate grows with A.

Then: the flume factorial, the bar (where Yang & Liu filter and limit the slope), and a 100-period
CN box regression gate in the suite, which today has no long-duration stability gate.

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
