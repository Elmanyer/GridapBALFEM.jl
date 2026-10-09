# STATUS.md — where the project stands (2026-10-09)

**Read this first.** One page on the state of the solver, the open problem, and what comes next.
Evidence for every claim is in [`STABILITY.md`](STABILITY.md) (the stability record) and
[`OPEN_ISSUES.md`](OPEN_ISSUES.md) (everything else that is open).

---

## 1. In one paragraph

The **v2 solver** (branch `v2-solver`) is correct and verified as far as the linear and the
nonlinear-without-`𝓝` models go. Its **full nonlinear model** (`nl_pressure=true`, all eight `𝓝`
components, Class III by the broken formulation) is **consistent but unstable at the grid scale**.
**The model without `𝓝` (`nl_pressure=false`) is unstable too** once the effective resolution is high
enough. The instability is a transfer of carrier energy into wavenumbers that `C⁰` elements represent
with no null and no dissipation, beyond what the reference finite-difference scheme can even
represent. **The standard remedy, the C⁰ interior penalty (first-order gradient jumps), is not
enough:** it removes the fastest grid-scale growth, but nonlinear spurious modes keep growing in the
mid band and the run fails later. Higher-order skeleton penalties held the v1 closed box, but they
are non-standard. Their required order grows with the element order (Gridap stops at order 2), and
they were never shown on a flume, a sloping bed, 2-D or MPI. **Direction adopted 2026-10-09:
filtering, as Yang & Liu do** (nine-point Shapiro filter of `η` and `u` every 20–200 steps).
Options are in §4.

## 2. What works and is verified

| area | status |
|---|---|
| v2 code base | Boolean `nl_pressure`; broken Class III only (exact Jacobian); v1 mixed, projected and mask code deleted; v1 knobs refused; `run_manifest.toml` per run |
| v1 → v2 equivalence | `test/v2_migration/regression_snapshot.jl` 12/12: residual and both Jacobians identical entry by entry (≤ 2e-16) |
| linear models 1–2 | optimal order in both fields, 1-D and 2-D, sequential and distributed, five vertical bases |
| nonlinear models 3–4 (no `𝓝`) | MMS at theoretical order to `nx = 32`; `p_η` degrades on finer meshes at Q3/Q2 (OPEN_ISSUES §1) |
| external oracle | at `p = 1` BALFE-M collapses onto Yang & Liu LFE-M coefficient by coefficient, `𝓝` included |
| Jacobians | `test_jacobians_ad` 14/14 on the patched Gridap; `∂R/∂u̇` exact, `∂R/∂u` quasi-Newton (rule 5) |
| suite | fast 9/9, `test_broken_formulation` 35/35, medium 20 pass (known failures: OPEN_ISSUES §2) |
| tooling | closed x-periodic box + `periodic_growth.jl` band rates; linear eigen-analysis of the penalised box operator (`cip_eigen_analysis.jl`, `cip_eigen_order_rule.jl`) |

**Not verified:** models 5–6 (`nl_pressure=true`) are outside the verified MMS scope (opt-in
`MMS_NL_P=1`; one exploratory run, `p_u` falls 3.94 → 1.97 on fine meshes, cause unknown).

## 3. The stability problem in numbers (v2, closed box, Crank–Nicolson, A = 0.10, 100 periods)

Divergence onset, or the high-band growth rate if the run stayed bounded for 160 s:

| | 8 cells/λ | 16 | 32 | 64 |
|---|---|---|---|---|
| Q2/Q1, `𝓝` on | bounded | bounded | growing +0.067 s⁻¹ | **81 s** |
| Q2/Q1, `𝓝` off | bounded | +0.03 | +0.024 | growing +0.08 |
| Q3/Q2, `𝓝` on | bounded (+0.002) | **82 s** | **27 s** | **9 s** |
| Q3/Q2, `𝓝` off | bounded | +0.03 | growing +0.07 | **44.6 s** |

Q2/Q1 at 16 cells and A = 0.15: 57 s with `𝓝`, bounded without. In every row, onset comes earlier
with refinement. That is the grid-scale signature: under-resolution would improve with refinement
(rule 38b). v2 fails 1.4–3× later than v1, but it does fail. Stabilisation history and analysis:
[`STABILITY.md`](STABILITY.md).

## 4. Next: filtering — the options

Requirement, unchanged since v1: a **wavenumber-selective energy sink** that removes the mid, high
and sub-element bands faster than the nonlinear transfer feeds them (≈ 0.1–0.5 s⁻¹ measured). It must
leave the carrier unattenuated over 100 periods, preserve the MMS orders, and work on both models.

| option | what it is | for | against |
|---|---|---|---|
| **F1 Shapiro-type nodal filter** | Yang & Liu's choice: a discrete low-pass stencil on nodal values every `n` steps | the reference scheme's own remedy; trivial to explain | needs a structured node lattice; Q_p nodes are not equispaced across cell types (vertex, edge, interior), so the stencil is not uniform; Cartesian only |
| **F2 modal / interpolation filter** (Fischer & Mullen 2001, spectral elements) | per cell, `u ← (1−α)u + α I_{p−1}u` (or damp the top Legendre modes) | made for `C⁰` high-order elements; keeps `C⁰` (traces interpolate consistently); local, cheap, any mesh, MPI-trivial; strength `α`, interval `n` | acts on the element's top mode, not on a wavenumber band, so mid-band (4–10 nodes/λ) content may need a stronger or repeated filter; must be checked against the carrier |
| **F3 differential (Helmholtz) filter** | solve `(I − δ²∇²)^r ū = u` (or a higher-order variant) every `n` steps | variational, mesh-agnostic, a sharp cutoff set by `δ`, `r`; distributed with the existing GMRES | one extra linear solve per filtered step; boundary conditions of the filter must be chosen |
| **F4 relaxation term in the residual** (time-relaxation, Stolz–Adams type) | `+χ (u − G u)` with a filter `G` (F2 or F3) inside the residual | continuous in time, no post-step hook, exact Jacobian if `G` is linear | changes the operator at every stage, with no interval to tune; `χ` must be calibrated |

Shared implementation points:
* F1–F3 need a **post-step hook**. Gridap's ODE `solve` is an iterator over `(t, uh)`, so filtering
  means restarting `solve(solver, op, t, tF, u_filtered)` every `n` steps, or wrapping the stepper.
  Under θ/CN the state is `u` alone. Generalised-α carries a velocity too.
* **Test bed:** the v2 unstabilised box runs above are the baseline: same launchers
  (`run/local/run_1dper_v2_unstab*.sh`), same seeds, `periodic_growth.jl`. Acceptance is the four
  criteria of [`STABILITY.md`](STABILITY.md) §6, all under Crank–Nicolson (rule 15).
* Then the flume, the bar (the lee-shoulder mode of OPEN_ISSUES §3 is exactly where Yang & Liu
  filter), MMS at the chosen strength, and a 100-period box regression gate in the suite.

The skeleton penalties (`:jumpgrad`, `:ghostvolume`) stay in the code as the comparison arm.

## 5. Other open items (details in [`OPEN_ISSUES.md`](OPEN_ISSUES.md))

* Nonlinear `p_η` order reduction at Q3/Q2 on fine meshes (advection block).
* Ghost penalty under AD (`_shift` loses the minus side); `test_selfconsistency` step 12.
* Distributed broken formulation and stabilisers (v2 plan step 10; the distributed driver refuses
  `nl_pressure=true`).
* Gridap fork commit `3928b98b9` not pushed; Manifest uses `Pkg.develop`.
* Cluster suite not re-run; sysimage stale.

## 6. Housekeeping state

* Branch `v2-solver`. The 2026-10-09 documentation restructure, the eigen order-rule driver and the
  v2 box launchers are committed; nothing is pushed.
* LaTeX v2 (`latex_docs/BALFEM_models_v2`, Overleaf): chapters 6–9 restructured. Chapter 8 has the v2
  Q2/Q1 rows; **the Q3/Q2 rows of §3 above are still to be added**. Chapter 9 describes the penalties
  and lists filters as "not tried". It will need rewriting around the filter.
