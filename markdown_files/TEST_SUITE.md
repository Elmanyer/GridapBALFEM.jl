# TEST_SUITE.md — the test campaigns and what each one can and cannot detect

> **Scope — what belongs here.** The *gate inventory*: every test file, its measured score, and —
> the part that matters — **what it would fail to notice**. Answers **"what is actually checked,
> and what would slip through?"**
>
> **What does NOT belong here:** whether the model is correct → [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md);
> defects → [`OPEN_ISSUES.md`](OPEN_ISSUES.md); future runs → [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md).
>
> ⚠ **"The suite passes" is not "the model is verified."** Most gates here are self-consistency
> checks, which pass for *any* residual. Only the analytic MMS is an oracle.


**What this file is.** The gate inventory: every test file, its current measured score, and — the
part that matters — **what it would fail to notice**. Scores below were measured in the 2026-08-18/19
full-suite run, not carried over.

Related: [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) (why the MMS tier is different in kind),
[`CONFIGURATION.md`](CONFIGURATION.md) (the settings under test).

**Headline: sequential 20/20 files · distributed 13/13 gates · Jacobian-vs-AD 17/17 over 8 models ·
nonlinear MMS 8/8 · `test/local/` 50/50.** One failure in the whole campaign
(`test_mms_convergence` G7), which is a **gate-window specification defect, not a solver defect**.

---

## 1. Tiers, and how to run them

```bash
# Sequential suite (per-file subprocesses, verdict from GATE OUTPUT)
julia --project=. test/runtests.jl                 # BALFEM_TESTS=fast|default|all|<files>

# A single file
julia --project=. test/test_basic.jl

# Distributed trio (4 ranks) — NOT reached by any runner; run by hand
~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. test/test_basic_distributed.jl
~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. test/test_nlpressure_distributed.jl
~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. test/test_bc_generation_distributed.jl

# Quasi-1D machinery suite (~45 min at JOBS=2)
bash test/local/run_local_tests.sh
```

> **`runtests.jl` takes its verdict from PASS/FAIL gate output, never from the exit code.** A file
> emitting no gate lines is reported **BLANK and counted as a failure**. This is not fussiness: a
> clean exit code is not evidence a test ran. One file was wrapped in
> `if abspath(PROGRAM_FILE) == @__FILE__`, so under `include()` it printed three header lines,
> returned "OK", and executed nothing — twice. Another carried a syntax error and could not parse,
> while its documented score sat in the docs for weeks. Both are fixed; the runner's design is what
> prevents a recurrence. `test/local/run_local_tests.sh` prefers each test's own `Results:` line,
> falling back to raw counting, with the BLANK check keyed off raw counts.

---

## 2. Sequential suite — 20/20 files

| file | score | what it gates |
|---|---|---|
| `test_vertical` | 15/15 | σ-tensor assembly: `A` symmetric, `Σ D = 1`, `w_j(0)=0`, `P[:,:,3] = −B` |
| `test_primitives` | 9/9 | index order, matvec, `⊗`, `⊙`, `double_contraction`, `∂x/∂y` orientation |
| `test_taylor_hood` | **13/13** | the horizontal pairing gate (added 2026-09-06). `Q2/Q1`, `Q3/Q2`, `Q4/Q3` accepted; **equal order rejected at three orders**, plus gap-of-two, η-above-velocity and `p_eta=0`; the rejection reachable through `build_fe_spaces` on the real path; and the error message checked for naming the CAUSE (inf-sup, the `2·dx` checkerboard) and a fix. ⚠ This gate exists because equal order caused the year-long "nonlinear instability" — `CLAUDE.md` rule 12b |
| `test_basic` | 6/6 | integration regression — references **bit-identical**: max η 0.00410, gauge 0.00212, Newton 240 |
| `test_dispersion` | 1/1 | phase speed vs model dispersion, 0.93 % |
| `test_dispersion_curve` | 9/9 | `C(kd)` across the band |
| `test_dispersion_nonlinear` | 3/3 | full-NL ⇒ Airy across the band |
| `test_shallow_water` | 6/6 | `kd→0` ⇒ `√(gd)`, 0.05 % against a 5 % gate |
| `test_sloshing` | 2/2 | standing-wave period vs BALFE-M theory, 1.93 % |
| `test_conservation` | 2/2 | closed-basin mass drift **7.8e-16** |
| `test_energy` | 3/3 | `ΔE/E₀ = 3.7e-14` (**pins `:theta`**) |
| `test_convergence` | 2/2 | temporal order `q ≈ 1.73` (**pins `:theta`**) |
| `test_nlpressure` | 10/10 | 𝓝 blocks — see §3 |
| `test_vertical_profile` | 7/7 | reconstructed `w(σ)`, `p_nh(1)=0` |
| `test_waveinput` | 30/30 | `WaveInput` construction, polarizations, dispersion solve |
| `test_bc_generation` | 11/11 | generated wave amplitude/phase (**pins `:theta`**) |
| `test_bc_spectrum` | 8/8 | multichromatic amplitude transfer, Goda–Suzuki (**pins `:theta`**) |
| `test_linear_newton_gate` | 10/10 | **1 Newton iteration per implicit stage on a SLOPING bed** — §4 |
| `test_jacobians_ad` | 17/17 | hand Jacobians vs AD, 8 models, amplitude-scaled — `VERIFIED_SCOPE.md` §5 |
| `test_mms_forcing` | 5/5 | forcing identities, no FE solve |
| `test_mms_forcing_nonlinear` | 20/20 | 𝓝 forcing: all tier×bed combinations finite, `flat_bed` control point, non-triviality, `O(A²)` excess |
| `test_mms_convergence` | G6 ✓ / **G7 ✗** | spatial `2.995`/`3.770` reproduced bit-identically; temporal gate mis-specified (§6) |
| `test_mms_convergence_nonlinear` | 8/8 | four studies — Models 3,4 `:none` and both `:native`, all at theoretical order |
| `test_selfconsistency` | 3/3 | **code support, not model validation** — see §5 |
| `test_equivalence` | RETIRED | its external per-layer reference predates the completion of the weak form (no `R_P`), so its result measured the reference's age. Correctly not counted |

`test_mms_convergence_nonlinear` runs **~2.5 h** — the `𝓝` forcing costs ~4× a `:none` evaluation.
Shorten with `MMS_NL_LEVELS`. The two `:full` studies are **deliberately not enabled** and the file
says why at length (`VERIFIED_SCOPE.md` §4).

## Distributed — 13/13 gates on 4 ranks (+ a fourth file added 2026-08-21)

| file | score | agreement with sequential |
|---|---|---|
| `test_basic_distributed` | 6/6 | 2.5e-7 |
| `test_nlpressure_distributed` | 3/3 | 1.6e-5 (`REF_EMAX = 0.0028640`) |
| `test_bc_generation_distributed` | 4/4 | 6.6e-7 |
| `test_mms_distributed_parity` | added 2026-08-21, not yet run | gates the A2 fix — see below |

**`test_mms_distributed_parity.jl` (4 ranks)** guards the 2026-08-19 defect: `run_mms_case_distributed`
hard-coded Model 1, so a distributed 8-model campaign returned 8 copies of it under 8 labels, all
passing. It compares the distributed and sequential branches on **Model 2** (moves `flat_bed` and
forces the general `mms_forcing` path) and **Model 5** (moves `regime` and `nl_pressure`), with the
reference recomputed in-process rather than pinned.

> ⚠ **Its G0 is a RESOLUTION gate, not a "the models differ" gate, and the distinction is the whole
> point.** Measured 2026-08-21 (nx=12, Q3/Q2, d=2.5): Model 1, 2 and 5 give
> `e_eta = 1.9844455e-4 / 1.9844463e-4 / 1.9844477e-4` — **separation ~4e-7**. **THE MMS L² ERROR IS
> A WEAK DISCRIMINATOR BETWEEN MODEL TIERS**, and structurally so: the MMS *forces* the same
> manufactured field for every model, so the exact solution is identical and the models differ only
> through the discretisation error of their different operators. `e_eta` must NOT be used to
> separate them; `e_u` separates by ~1e-4 and is the usable channel. G0 therefore asserts
> `separation > 30 × PARITY_RTOL`, i.e. that the test *can resolve* the defect it exists to catch.
> If that ever fails, the fix is a **more discriminating case** (`nl_pressure=:full`, whose `e_u`
> sits ~600× higher), never a looser parity tolerance.

## `test/cluster/`

`cluster_conservation` (2 ranks, drift 5.8e-10) and `cluster_selfconsistency` + a SLURM template.
The latter is named for what it is: **its forcing is the solver's own residual**, which is exactly
why `test_mms.jl` was renamed `test_selfconsistency.jl` — leaving it named "mms" reintroduced the
confusion that rename existed to kill.

## `test/local/` — the quasi-1D machinery suite, 50/50

| gate group | score | result |
|---|---|---|
| rest state (flat + sloping bed) | 8/8 | `max|η| = 0` **exactly**; mass drift 0 |
| sponge damping law | 18/18 | `R² = 0.9997 / 0.9996 / 0.9991`; reflection 1.0–1.4 % |
| relaxation zone | 9/9 | in-zone amplitude error ≤ 0.1 %, phase 0.9°, absorption **145×** its control |
| boundary modes | 16/16 | includes a **negative control that must fail** — and does, in 6 s of simulated time |
| `test_2d_reduces_to_1d` | 8/8 | a y-invariant 2-D run reproduces the quasi-1-D flume |

`test_2d_reduces_to_1d` deserves its own note: it tests the solver against a **symmetry the model
must respect** rather than against theory via its own residual, and it exercises the `Ey`/`𝖴y`
machinery that no 1-D case touches. **It is the sharpest check available short of the analytic MMS.**

---

## 3. What `test_nlpressure` proves, block by block

| gate | proves | to |
|---|---|---|
| G1 | the **∇h half** of `𝓝{1,2,4,5}` — exact-IBP identity against *analytic* second derivatives, on a deliberately asymmetric state that would catch a k/j slot swap | **4.0e-15** |
| G2 | every block's nonlinear **order** by amplitude scaling: ∇h-IBP → 4, ∇H-frozen → 8, 𝓟-frozen → 4 | ratios 4.0001 / 8.0002 / 4.0001 |
| G3 | dynamics over a tanh bar, `:nonlinear` + `:full` + `flat_bed=false` — the one **sequential** configuration that can see the bed-slope physics | `REF_EMAX = 0.0028640`, `REF_RTOL = 1e-4` |

G3 previously asserted only `emax < 20A` — **boundedness**. It passed 9/9 while the quantity it
computes moved **58 %** under a residual fix, because that gate had ~7× headroom. It now pins a
value, measured at 1 % first and tightened 100× once the `%.7f` print gave the digits.

What is still **not** value-checked: the ∇H-frozen and 𝓟-frozen halves, and any dynamic
confirmation at realistic amplitude. See §7.

---

## 4. Coverage rules — the part worth more than the scores

> **1. A BOUNDS CHECK ON THE RIGHT CONFIGURATION IS NOT A VALUE CHECK.** Running the right physics
> proves nothing if the assertion is only that nothing exploded.

> **2. For a guard that is a CONJUNCTION (`nonlinear ∧ ∇h ≠ 0`), the suite needs a case satisfying
> the conjunction AND asserting a value.** Conjuncts satisfied separately, or satisfied together but
> only bounded, both look green.

This is not abstract. Most of the suite structurally cannot reach the bed-slope terms:
`test_sloshing` is `:linear`; `test_shallow_water` is `flat_bed=true`; `test_conservation`,
`test_energy` and `test_basic` use a **constant** `h_bathy` (so `∇h = 0` numerically);
`test_linear_newton_gate` uses a sloping bed but is `:linear` by construction. A missing nonlinear
`∇h` term therefore sailed through 21 green files.

The configuration is now covered **quantitatively** by three independent things — the analytic-MMS
Model 4 study, `test_jacobians_ad.jl` rows M4/M6/M8, and `test_nlpressure` G3 / its distributed twin.
**Do not delete any of them without replacing the coverage.**

> **3. `test_linear_newton_gate` is free, sharp, and needs no reference value.** Under
> `regime=:linear` the residual is affine in `(u,u̇)`, so exact hand Jacobians reach round-off in
> **one Newton iteration per implicit stage** from any guess, at any amplitude, over any bathymetry.
> Any excess is proof of a residual↔Jacobian inconsistency. It is run on a **sloping** bed, which
> closes the structural blind spot above. Measured: exactly 1 iteration/stage (θ: 1/step, SDIRK:
> 2/step), residual 1.9e-15…5.4e-15.

> **4. The default integrator is DISSIPATIVE — a test measuring a non-dissipative property must pin
> `solver_type=:theta`.** Measured difference:

| test | measures | SDIRK_2_2 | `:theta` (CN) |
|---|---|---|---|
| `test_energy` | non-dissipativity | `ΔE/E₀ = −1.34e-2` | **`+3.65e-14`** |
| `test_bc_spectrum` | amplitude transfer, 3 components | 12.5 / 21.9 / 32.2 % | **2.5 / 4.6 / 9.4 %** |
| `test_bc_generation` | generated wave amplitude | 23.1 % | **8.4 %** |
| `test_convergence` | temporal order (Richardson) | **0.01** | ≈2 |

**Do not remove those pins, and never "fix" such a failure by moving a threshold** — the test's
subject is a property CN has and SDIRK deliberately lacks.
**Recognising this failure mode:** amplitude damped while **phase is correct**, error **growing with
frequency**, and — decisively — **refining the mesh does not help** (one such test got *worse* at
high frequency when refined, 32.2 → 41.6 %). That last check separates it from genuine
under-resolution, which is a real and separate failure mode here: a celerity gate once read 6.47 %
because it meshed at 6 cells/λ, and refining **space** alone fixed it (→0.05 %) while refining
**time** alone changed nothing. **A physics gate is only as sharp as the discretisation feeding it.**

---

## 5. `test_selfconsistency.jl` — what it is and is not

Its forcing is `f = R(u*)` built with the solver's **own** residual, so writing `R = R_true + E`
the error cancels identically and it passes for **any** `R`; it also measures no convergence rate.

It is kept because it certifies what nothing else does — that the hand Jacobians are the exact
derivatives of the residual and the multi-step bookkeeping is consistent — but as **code support,
not model validation**. Its file header says so at length.

---

## 6. The one open gate: `test_mms_convergence` G7 (temporal order)

Measured `p_η = 1.382`, `p_u = 1.332` against a `2 ± 0.3` gate. **Not a solver defect**: G6
(spatial) reproduces its documented `2.995 / 3.770` bit-identically, and the verified-scope table is
entirely spatial.

**The window is contaminated at BOTH ends at once** — pre-asymptotic at the coarse end, spatially
floored at the fine end — which is why every intervention moved one field and not the other. Seven
configurations were measured (recorded in the test file). **`u` is second order beyond doubt**
(pairwise `1.992/2.056/2.016`, reproduced in a second window); **`η` is too, but its asymptotic
window is ONE refinement wide** — over-converging above `dt ≈ 0.0375` (rate 2.631) and saturating on
its spatial floor below `dt ≈ 0.01875` (rate 1.157).

> ⚠ The best window found (`T=2.4`, `nx=36`, `dt0=0.075`) reads `p_η = 1.974`. **Do not trust that
> number** — it hits the target only because coarse-end over-convergence and fine-end saturation
> cancel in the least-squares fit.

**The cause is structural, a cost of the Taylor–Hood pairing.** Under `Q3/Q2`, `η` sits in the
*lower*-order space, so its spatial floor is relatively higher and squeezes its temporal window from
below. **The pairing that optimises the spatial study pessimises the temporal one for `η`.**
Options (re-specify at ~3× runtime, gate the fields differently, or refine to `nx ≥ 48`) are left
open in the test file. **Do not widen the ±0.3 tolerance — the window is wrong, not the expectation.**

**G10 is also under-powered** — it is the guard that should have caught this. It measured **3.4 %**
contamination for `u` against **0.5 %** for `η`, a **sevenfold** asymmetry, and passed both on one
shared 10 % threshold. It already computes the two numbers; it just does not act on them separately.

> **General rule: a mirrored guard that collapses two per-field measurements into a single pass/fail
> can hide the very asymmetry it exists to detect.**

---

## 7. The resolution principle

> **A test validates a term only if the test can RESOLVE that term's contribution. Always ask the
> counterfactual: if this term were wrong, would this test have noticed?**

Worked example, `nl_pressure=:full`. Switching it off changes `max|η|` by **0.0129 % (1-D) /
0.094 % (2-D)**, so the local observation set establishes that it *runs, converges and stays
bounded* — useful, since it is the path that once NaN'd on the cluster — but carries **no
information about whether its coefficients are right**.

Say *"unvalidated by the local set"*, never bare *"unvalidated"*: the ∇h half **is** validated to
machine precision by G1, and G2 pins every block's nonlinear order.

**The sensitivity is tunable and the missing experiment is cheap.** The block is `O(A³)` against
`O(A)` leading terms, so its relative effect scales as `A²`: 0.013 % at `A = 1e-3`, **~1.3 % at
`A = 1e-2`** — comfortably measurable. The 2-D bar case already ran at `A = 0.01` with `:full`; its
`:none` counterpart was never run, and that single pair would both resolve the difference and
confirm the `A³` scaling dynamically.

---

## 8. Known gaps in the suite itself

> ### ⚠ THE SUITE CANNOT DETECT THE NONLINEAR GRID-SCALE INSTABILITY (added 2026-09-05)
>
> The equal-order mode documented in `CLAUDE.md` rule 12b (now cured by the Taylor-Hood
> requirement, but the coverage gap it exposed is real) emerges only after
> **50–80 s** of simulated time at production amplitudes, and *every* nonlinear test in the suite is
> far shorter. It therefore reached production runs with a fully green suite behind it.
>
> This is the **resolution principle** (§7) failing in the time axis rather than the space axis: the
> tests are not wrong, they simply cannot resolve the phenomenon. Ask of any stability test, *if the
> discretisation had a slowly growing mode, would this test have run long enough to see it?*
>
> Two further traps this exposes, both of which would make a new gate green and useless:
> * **A short nonlinear run is not evidence of stability** — the growth time at `A=0.05` is ~84 s.
> * **A coarse mesh flatters the result** — growth is *slower* at `dx=0.50` than at `dx=0.125`, so a
>   cheap regression on a coarse mesh is the least sensitive test that could be written. Any gate
>   must pin `dx` fine enough and `dt` small enough that integrator dissipation is not doing the
>   stabilising (rule 12c).



* **The MPI tests are not reached by any runner.** `runtests.jl` covers the sequential suite; the
  distributed trio is only ever run by hand. That cost three stale reference constants, all found
  the first time they were re-run. **Fix: give `runtests.jl` an MPI tier that runs them when
  `mpiexecjl` resolves and SKIPS LOUDLY when it does not** — a silent skip recreates the blind spot.
* Only **one vertical resolution and one `kd`** are exercised locally (P1LFE-2, `Nσ=3`, `kd=5.5`).
  P1LFE-3/4 and the rest of the band are untouched. The sponge in particular is verified at
  `kd = 5.5` only, while it is *required* to cover the **longest** component — a `T` sweep measuring
  reflection at `kd = 1, 3, 5.5` would test that directly.
* No run-and-reconstruct **pressure** profile test; no distributed-gauge utility.
