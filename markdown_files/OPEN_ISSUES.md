# OPEN_ISSUES.md — everything open except stabilisation, each with its next step

*Rewritten 2026-10-09. The stability problem is the first open item; it has its own file,
[`STABILITY.md`](STABILITY.md), and the decision is in [`STATUS.md`](STATUS.md). The full v1 version
of this file, with every measurement, is [`archive/OPEN_ISSUES_v1.md`](archive/OPEN_ISSUES_v1.md);
the section names below give its matching section.*

---

## 1. 🔴 Nonlinear `p_η` loses an order at Q3/Q2 on fine meshes (v1 §0b)

* The pairwise rates for `nx` = 4…64:
  * nonlinear: 2.971 → 2.922 → 2.757 → **2.450**;
  * linear: 3.000.

  The error still falls; the *rate* decays with refinement. `p_u` rises at the same time
  (2.77 → 3.69).
* Reproduced in 2-D (2.5565). Invisible at Q2/Q1, where it is degenerate with the optimal 2, and at
  Q4/Q3, where the ladder is too short.
* **Eliminated:** algebra (Newton to 1e-15), `𝓝`, `∇h`, the linear core, quadrature, 1-D posing.
  That leaves the **nonlinear advection block** (`𝓜₃/𝓖`).
* **Next:** derive that block term by term against the LaTeX. Cheap discriminators: amplitude
  scaling (`a_eta` = 0.4 and 0.2), or extend Q4/Q3 to `nx = 64`.
* The verified scope stands as measured (to `nx = 32`), with this qualifier.

## 2. Known test failures

| test | state | next step |
|---|---|---|
| `test_skeleton_ad` S2 `:ghostvolume` (0.75), S3 (Newton NaN) | `_shift` in `src/broken.jl` has no method for `SkeletonCellFieldPair` or transient wrappers, so AD loses the minus side. Our code, not Gridap; hand-Jacobian runs are unaffected | add the methods (needs approval: existing library) |
| `test_selfconsistency` step 12 | Newton converges (‖r‖ 1.7e-15) to a **different discrete root**, 27 % off, at Q2/Q1, large amplitude, `𝓝` on. New in v2 | decide the test configuration (Q3/Q2, or a smaller amplitude) |
| `test_mms_convergence` G7 | temporal window contaminated at both ends (`TEST_SUITE.md` §6); a specification decision, not a bug | re-specify (≈ 3× runtime) or gate the fields separately. **Do not widen ±0.3** |
| `test_bc_generation` case C | 10.1 % against a 10 % gate | carried from v1; unexamined |

## 3. ⛔ Any variable bed grows a lee-shoulder mode (v1 §0d)

* Tested with `𝓝` off (v1 `:native`), A = 0.10, a bar of height 2.0 m on `d` = 3.5 m. Onset 37 /
  58 / 94 / 137 s for `max|∇h|` = 2.0 / 1.0 / 0.67 / 0.5; the flat bed reaches 100 periods.
* It is a **rate** set by the slope, not a threshold. The growth is pinned at the downwave shoulder,
  with `u_max` 2–3 against 0.42.
* **Next:**
  * halve `dx` on one bar case. Earlier onset means a grid-scale problem; later onset means
    under-resolution (rule 38b);
  * then test the filter there: Yang & Liu filter, and limit the slope, on exactly this case.
* Never call a bar run stable before 100 periods (the 1:2 bar "stable at 80 s" died at 137 s).

## 4. Full nonlinear model (`nl_pressure=true`) outside the verified scope

* MMS models 5–6 are opt-in (`MMS_NL_P=1`). In one exploratory Model-5 run, `p_η` was optimal
  (2.997) but **`p_u` fell 3.94 → 1.97** for `nx` = 8–32. The cause is not established: either the
  trial Hessians are inherently O(h^{p−1}), or there is a defect.
* **Next:** after the filter is chosen, run the models 5–6 MMS with and without it.

## 5. Distributed path behind the sequential one

* **Step 10 of `V2_SOLVER_PLAN.md`:** the broken Class-III formulation and the stabilisers (and now
  the filter) on `DistributedTriangulation`. `setup_and_run_distributed` refuses `nl_pressure=true`.
* No MPI tier in `runtests.jl`. The four distributed tests are run by hand; this has cost three stale
  constants. Fix: a tier that runs when `mpiexecjl` resolves and **skips loudly** otherwise.
* **Preconditioner** (the biggest performance item). GMRES+Jacobi needs 300–770 iterations; the
  cheap drop-ins are measured dead (`CONFIGURATION.md` §5). Next to try: field-split/Schur, then
  geometric multigrid.
* **Cluster memory.** `--mem-per-cpu=4G` is required. The leading explanation is the baseline
  footprint (≈ 1.4 GB per Julia+Gridap process), but a slow growth over hours is also documented
  (rule 41). Unresolved.

## 6. Infrastructure

* **Gridap fork commit `3928b98b9`** (skeleton AD on transient fields) is **not pushed**. The local
  Manifest loads it via `Pkg.develop(path="Gridap.jl")` (backup:
  `output/v2_logs/Manifest.toml.before_gridap_develop`). Push it, restore the URL pin, and describe
  the patch in `CONFIGURATION.md` §2.
* **Sysimage stale; cluster production suite not re-run** since v1.
* `examples/local_1d/run_flume_1d.jl` ends with an undefined `tag`, so successful runs exit 1
  (rule 38h). One-line fix.
* Record the v2 migration values in `V2_SOLVER_PLAN.md` §3 and `TEST_SUITE.md`.
* Run-output gaps:
  * the in-run discrete-equation check runs under θ only;
  * no per-Newton-iteration linear counts;
  * no wall-time breakdown;
  * the achieved linear tolerance is not recorded.
* GitHub repository still named `GridapLFEM.jl`.

## 7. Verification gaps (lower priority)

* Q2/Q1 velocity rate 2.003 against an optimal 3, on every model. The vertical-basis analogue was
  pre-asymptotic; the horizontal pairings were never re-run on an extended ladder.
* Error versus DOF at production resolution: a prerequisite for changing the default pairing
  (`Q3/Q3` was 40× more accurate than `Q3/Q2` at `nx = 24`, before Taylor–Hood was enforced).
* Only P1LFE-2 at `kd = 5.5` has been exercised in long runs.
