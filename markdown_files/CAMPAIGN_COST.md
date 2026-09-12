# CAMPAIGN_COST.md — measured run costs, memory bands, and scheduling lessons

**Written 2026-09-03**, from the vertical-basis campaign, the Phase-B MMS batch (paused at
28/42), and the small-domain cluster suite. Every number here is **measured**, not estimated —
each is traceable to a `diagnostics.csv`, a job `.out`, or a supervisor log. Where something is
an extrapolation it says so.

> **Purpose.** The Phase-B batch produced 28 studies in ~37 h of wall time on a 16-core
> workstation and had to be paused. This file exists so the next campaign is *scheduled* rather
> than discovered: what a study actually costs, how much memory it needs and for how long, where
> the cluster ceilings are, and which of our own decisions cost the most time.

---

## 1. Headline numbers

| quantity | value |
|---|---|
| MMS study, cheapest observed (`T9_tier2_1d`, Q2, 5 levels) | **24 min** |
| MMS study, dearest observed (`T8_full_projection`, P1LFE-4) | **32.0 h** |
| spread within one family | up to **11×** |
| cluster case, linear, 32 ranks, 8000 cells | **24–27 s/step** |
| cluster case, nonlinear `:full`, same mesh | **143–248 s/step** (≈ **9×**) |
| worker RSS, fresh → after 20 h | **0.8 → 4.6 GB** |
| cluster peak RSS per rank | **2501–3849 MB** ⇒ 4 GB/core required |
| useful concurrency, 16-core / 31 GB workstation | **6–8 heavy MMS workers** |

**The single largest avoidable cost in the whole campaign** was running the 1-D MMS studies with
`ny_1d = 3` instead of `1`. See §7 — it is worth ~3× in DOFs and far more in wall time, for
**bit-identical answers**.

---

## 2. MMS study cost (Phase B, sequential, one worker per study)

Wall-clock seconds per *study* (a study = a whole convergence ladder, 3–5 solves), from the 28
completed rows:

| family | n | min | mean | max |
|---|---|---|---|---|
| `T9_tier2_1d` (Q2/Q3/Q4, 1-D, 5 levels) | 6 | 1 465 s (24 min) | 12 141 s (3.4 h) | 31 206 s (8.7 h) |
| `T10_highorder_newnodes` (p≥2 bases, 4 levels) | 12 | 3 918 s (1.1 h) | 25 469 s (7.1 h) | 65 403 s (18.2 h) |
| `T7_extladder` (P1LFE-4, 5 levels to nx=128) | 4 | 17 748 s (4.9 h) | 30 605 s (8.5 h) | 41 088 s (11.4 h) |
| `T8_full_projection` (`:full`, 4 levels) | 6 | 28 723 s (8.0 h) | 55 530 s (15.4 h) | **115 387 s (32.0 h)** |

### Campaign C + C3 (2026-09-11/12, workstation, 16 cores / 31 GB)

Per *study*, `nx0=4`, `a_eta=0.8`, all P1LFE-2. 18 × 1-D (5 levels for Q2Q1/Q3Q2, 4 for Q4Q3) and
the 2-D Q4/Q3 tier (4 levels):

| family | n | mean | max | peak RSS |
|---|---|---|---|---|
| Q2/Q1, 1-D | 6 | 96 min | 164 min | ~1.6 GB |
| Q3/Q2, 1-D | 6 | 181 min | 347 min | ~1.7 GB |
| **Q4/Q3, 1-D** | 6 | 165 min | 342 min | **4 970 MB** |
| **Q4/Q3, 2-D** (short ladder `nx≤32`) | 2 so far | 263 min | 332 min | ~3.6 GB |
| Q4/Q3, 2-D, **LONG** ladder `nx≤64` | — | *never completed* | — | **12 100 MB** |

⚠ **`PLANNED_CAMPAIGNS.md` §5.4 said "1-D is the cheap end … minutes to ~1 h". That is wrong and it
mattered.** The true 1-D figure is **1.5–6 h per study**, and it is the *regime*, not the dimension,
that drives it: a LINEAR Q3/Q2 study runs in 30 min while the NONLINEAR one on the same ladder takes
183–347 min — a 6–11× factor, because every step costs a full Newton solve instead of one iteration.
Estimating a nonlinear campaign from a linear timing under-predicts by an order of magnitude.

⚠ **Q4/Q3 IS THE MEMORY CLIFF, AND ITS COST IS NOT IN THE DOF COUNT.** A 1-D Q4/Q3 study peaks at
**4 970 MB** on a mesh with ~4 000 DOFs — three times the Q3/Q2 figure on a problem of comparable
size. The cost is the assembly/specialisation machinery for high-order elements with the
`ThirdOrderTensorValue` contractions, which scales with element order and quadrature points, not
with the mesh. Consequences for sizing:

* **budget ~5 GB per concurrent Q4/Q3 worker, 1.6–1.7 GB for Q2/Q1 and Q3/Q2.** On a 31 GB box that
  is ~5 concurrent Q4/Q3 studies, not 8, whatever the core count says.
* **the 2-D LONG ladder (`nx0=8`, to `nx=64`) reached 12.1 GB** before it was killed — a further
  reason it is cancelled in place (`PLANNED_CAMPAIGNS.md` §1), independent of its saturation problem.

⚠ **FIRST-STEP JIT DOMINATES THE COARSE LEVELS AND LOOKS EXACTLY LIKE A HANG.** Measured: step 1 of
the first level costs **175 s** of `t_solve` on a *linear* Q3/Q2 case against 0.28–0.32 s for every
subsequent step; on the nonlinear residual with hand Jacobians at Q3/Q2 it is **10–15 minutes**, and
the only external sign of life is that `diagnostics.csv` sits at header-only. Four probes launched
together all showed empty diagnostics for 10–15 min at 100 % CPU. Before diagnosing a stall, check
`cputime ≈ elapsed` (compiling counts as work) and whether `src/*.jl` was edited recently — an edit
invalidates the precompile cache and every worker then pays its own rebuild.

### The cost drivers, in order of importance

1. **Pressure tier.** `:full` is the dearest family by a wide margin — mean 15.4 h against 3.4 h
   for the `:none` tier-2 studies. The `𝓝` forcing carries all eight components and the
   `Nσ²×8` block the outer gradient then differentiates.
2. **Vertical richness `Nσ`.** Forcing cost scales ≈ `Nσ²`. P1LFE-4 / P2LFE-2 (`Nσ=5`) studies
   are the slowest rows inside every family.
3. **Ladder length.** A 5th level at `nx=128` costs roughly as much as the preceding four
   combined — `h`-refinement in 2-D quadruples cells per level.
4. **Nonlinearity.** Newton iterations rise 2 → 6–12 per step (§4), and each is a full assembly.

### Planning rule of thumb

```
study_hours  ~  0.5  x  tier_factor  x  (Nsigma/3)^2  x  2^(levels-3)
    tier_factor:  :none 1,  :native 3,  :full 8
```
Checked against the table: `T8` P1LFE-4 (`:full`, `Nσ=5`, 4 levels) → `0.5·8·2.8·2 = 22 h`
against a measured 20.2–32.0 h. Good enough to schedule with; do not quote it as a result.

---

## 3. Cluster production-run cost (small-domain suite, Snellius `rome`)

All 32 ranks / 8000 cells unless noted; `dt=0.02`.

| case | ranks | cells | s/step | steps | wall | peak RSS |
|---|---|---|---|---|---|---|
| bcplane linear none flat | 32 | 8 000 | 23.9 | 1601 | 10:37 | 2625 MB |
| bcplane linear none bar | 32 | 8 000 | 26.6 | 1601 | 11:50 | 2501 MB |
| plane linear none flat A=0.001 | 32 | 8 000 | 27.4 | 1401 | 10:40 | 2607 MB |
| plane linear none flat A=0.1 | 32 | 8 000 | 27.3 | 1401 | 10:38 | 2635 MB |
| plane nl `:full` flat A=0.001 | 32 | 8 000 | 50.4 | 1401 | 19:36 | 2807 MB |
| directional linear none | 48 | 16 000 | 38.3 | 2600 | 27:41 | 2642 MB |
| bcplane nl `:full` flat A=0.1 | 32 | 8 000 | 142.8 | 1601 | **63:30** | 2882 MB |
| bcplane nl `:full` bar A=0.1 | 32 | 8 000 | 186.8 | 1555 | **80:42** | 3849 MB |
| plane nl `:full` flat A=0.1 | 32 | 8 000 | 248.5 | 488 † | 33:41 | 2779 MB |
| plane nl `:full` bar A=0.1 | 32 | 8 000 | 238.2 | 306 † | 20:15 | 3633 MB |
| directional nl `:full` | 48 | 16 000 | 245.6 | 157 † | 10:43 | 2875 MB |
| irregular nl `:full` flat | 32 | 8 000 | 233.5 | 690 † | 44:45 | — |
| irregular nl `:full` bar | 32 | 8 000 | 363.0 | 714 † | 72:00 | — |

† diverged before `T_final` — the wall time is what was *spent*, not what the case needs.

**Observations**

* **Linear vs nonlinear `:full` on the identical mesh: 24 → 143–248 s/step, a factor 6–10.** This
  is the dominant planning number for the cluster suite. Budget a nonlinear `:full` case at
  **~10× its linear twin**.
* **Amplitude costs nothing by itself.** `plane linear A=0.001` and `A=0.1` are 27.4 and 27.3
  s/step — identical. Cost comes from the *regime*, not the amplitude.
* **Bathymetry adds ~15–30 %** (`bcplane` flat 142.8 vs bar 186.8).
* **A diverging case is not cheap.** The bar case burned 20 h to reach 22 % of `T_final`. Diverging
  runs get *slower* as Newton struggles: 238 s/step average against 143 for the case that survived.
* **The 119 h wall-clock limit is the real constraint.** `bcplane nl :full bar` used 80:42 of it
  for 1555 of 1400 requested steps. A nonlinear `:full` case much larger than the small domain
  **will not fit in one job**.

---

## 4. Solver-level cost breakdown

Averages over each run's `diagnostics.csv`:

| case | Newton/step | GMRES/solve | s/step | `lin_sat` |
|---|---|---|---|---|
| linear (all four) | **2.0** | 319–387 | 24–38 | 0 |
| nl `:full` A=0.001 | **2.0** | 402 | 50 | 0 |
| nl `:full` A=0.1 (bcplane) | **6.4** | 360 | 143 | 0 |
| nl `:full` A=0.1 (plane, diverging) | **11.6** | 374 | 238 | 0 |
| nl `:full` directional (diverging) | 6.4 | 439 | 171 | 0 |

* **GMRES count is nearly constant (320–440) across every case** — linear or nonlinear, converging
  or diverging. The linear solve is *not* what makes nonlinear runs expensive.
* **Newton count is.** 2.0 → 6.4 → 11.6 tracks the cost almost exactly. Cost per step ≈
  `Newton_iters × assembly`, and each Newton iteration is a fresh Jacobian assembly plus a GMRES
  solve of ~370 iterations.
* **`lin_sat = 0` everywhere**: GMRES never hit its 1000 cap during normal operation. `ls_maxiter`
  is not a binding constraint; only diverging runs saturated it, and only *after* Newton failed.
* **A rising Newton count is the earliest divergence warning available** — it climbs several
  hundred steps before `eta_max` misbehaves. Worth an automatic flag.

### Linear-solver tolerance (fixed 30-step 2-D case)

| `ls_rtol` | GMRES/solve | s/step | final `eta` |
|---|---|---|---|
| 1e-6 | 299 | 55.3 | 2.6234e-03 |
| 1e-7 | 363 | 69.3 | 2.6234e-03 |
| 1e-9 | 497 | 103.0 | 2.6234e-03 |

**`ls_rtol` costs 86 % more wall time from 1e-6 to 1e-9 and changes the answer in no printed
digit.** Production default 1e-5 is right; tightening it is pure loss.

---

## 5. Memory bands

### Per-process, sequential MMS worker (workstation)

| age | RSS |
|---|---|
| fresh (post-JIT) | 0.8 – 1.5 GB |
| after ~4 h | 2.0 – 3.0 GB |
| after ~16–20 h | **3.0 – 4.6 GB** |

Matches the vertical-basis campaign's independent measurement (1.5 → 3.9 GB over ~14 h, up to
4.6 GB). **This drift is unbounded within a run** and is the reason for bounded worker lifetime.

### Per-rank, distributed cluster run

Peak RSS **2501 – 3849 MB**, i.e. 61–94 % of a 4 GB/core allocation. Two distinct facts, often
conflated:

* **Within a run, RSS rises for ~100 steps then PLATEAUS** (plane 1894→~3230 MB then flat;
  directional 1546→~2800 MB then flat). There is **no per-step leak**.
* **Across a long-lived worker, RSS drifts without plateau.** That is a different phenomenon and
  applies to the workstation MMS workers, not to a single cluster job.

**Consequence: 4 GB/core is required and permanent.** 2 GB/core was tested and OOM-killed;
3 GB/core would OOM the 3849 MB case. Baseline Julia+Gridap is ~1.4–1.5 GB before solving
anything, so 2 GB/core leaves ~0.5 GB for the computation.

### Per-cell marginal cost (for sizing new meshes)

Subtracting the ~1450 MB baseline:

| tier | MB per cell per rank |
|---|---|
| linear (directional, 333 cells/rank) | **≈ 3.6** |
| nonlinear `:full` (bcplane, 250 cells/rank) | **≈ 5.7** |

Sizing formula: `RSS_MB ≈ 1450 + per_cell × (cells / ranks)`. Checked: 32000 cells / 40 ranks
linear → 1450 + 3.6·800 = **4310 MB**, which is why that case needed 5 GB/core.

---

## 6. Cluster queue and scheduler limits (Snellius `rome`)

Measured by submission, not read from documentation:

| request | GB/node | outcome |
|---|---|---|
| 32 × 4 GB | 128 | ✅ accepted |
| 48 × 4 GB | 192 | ✅ accepted |
| 40 × 5 GB | 200 | ✅ accepted |
| 64 × 4 GB | 256 | ❌ **refused** — "Requested node configuration is not available" |

* **A rome node advertises 256 GB but SLURM can allocate only ~224 GB.** Any request built as
  "N ranks × 4 GB = the node's nominal RAM" fails. Eleven launchers had this bug.
* **This site additionally does not admit >48 cores per node together with the raised per-core
  memory.** So the usable envelope is **≤48 cores/node AND ≤224 GB/node**.
* Raising `mem-per-cpu` therefore *forces the rank count down*: 48 × 5 GB = 240 GB is refused, and
  40 × 5 GB = 200 GB is the largest that fits at 5 GB/core.
* Wall-clock limit **119:59:00**. A nonlinear `:full` case at 143–248 s/step fits ~2000–3000 steps
  in one job. Budget accordingly or checkpoint.

---

## 7. Core saturation and concurrency (16-core / 31 GB workstation)

### The measurement that matters most

`ny_1d` for 1-D MMS studies, static Q3/Q2 at `nx=16`:

| | free DOFs | `e_eta` | `e_u` | s/step |
|---|---|---|---|---|
| `ny=3` | 2817 | 3.365575179792e-05 | 1.302487479913e-06 | **630.3** |
| `ny=1` | 957 | 3.365575179792e-05 | 1.302487479915e-06 | **3.4** |

**Identical to 14 and 12 significant figures, at 2.94× fewer DOFs — and ~180× less wall time**
in that instance (a 3-cell-wide `Q3` mesh widens the LU front far more than the DOF ratio
suggests). The `:d1` manufactured field has `ky=0`, so `u*_y ≡ 0` and the solution is y-invariant;
`mms.jl` already said "refine nx only; refining ny changes nothing".

**36 of the 42 Phase-B studies are `:d1`.** Running them at `ny=3` was the campaign's single
largest self-inflicted cost.

### 1-D MPI is a loss, decisively

Same 240×3 mesh, 20.2k DOFs:

| ranks | 1 (LU) | 2 | 4 | 6 | 12 |
|---|---|---|---|---|---|
| s/step | **1.36** | 10.18 | 8.46 | 8.23 | 18.27 |
| GMRES | — | 758 | 758 | 758 | 715–774 |

GMRES count is **rank-independent**, so the cost is the weak Jacobi preconditioner, not the
decomposition. **Direct LU beats every MPI split.** 2-D scales, but modestly: 98.2 → 77.7 → 67.6
s/step for 2 → 4 → 6 ranks (1.45× for 3× the ranks).

### Useful concurrency

| workers | outcome |
|---|---|
| 12 | ❌ RSS drift → swap 7850/8191 MB → all workers GC-bound; **worse than a third as many** |
| 10 | ⚠ 4.7 GB free, watchdog fired 23×, **2 rows in 13.5 h** |
| 8 | ✅ 9.3 GB free, stable |
| 6 | ✅ comfortable |

**6–8 is the usable band for heavy MMS workers.** The binding resource is memory, not cores:
at 8 workers load was ~12 on 16 cores, so cores were *not* saturated — memory was.

Cheap 1-D cases (`ny=1`) are different: 16 concurrent flume runs sit at 11.2 GB total and load 12.7.

---

## 8. Where the time actually went — churn accounting

Phase B, 19:08:39 → 08:27:47 (≈37 h), 28 rows of 42:

| cause | worker launches |
|---|---|
| **exhausted-shard spin (a bug)** | **219** |
| watchdog kills → relaunch | 23 |
| legitimate recycling (1 per study) | ~18 |
| **total** | **258** |

* **219 of 258 launches were one shard relaunched after exiting instantly.** The supervisor picked
  the lowest free id; a shard whose queue was empty exited in seconds and was free again next
  cycle. It also **starved shards 8–11**, which held real work. Fixed with exhausted-markers.
* **23 watchdog kills destroyed 23 in-flight studies.** The watchdog assumed killing the fattest
  worker restores headroom faster than work is lost. That holds for short studies; with the
  remaining families at 1.5–32 h it inverted, and between 10:44 and 17:08 exactly **one** row
  landed while 22 kills happened.
* **A watchdog that fires constantly is a sizing signal, not a solution.** The fix was 8 workers,
  not a better kill policy.

---

## 9. Recommendations for the next campaign

**Before launching**

1. **Set `ny_1d = 1`** for every `:d1` study. Free 3× DOFs. *(Now the default.)*
2. **Cost-tag every job honestly** and schedule longest-first (LPT) if the deliverable is the whole
   matrix, shortest-first (SJF) only if partial coverage has value. Note SJF schedules the cheap
   high-value studies last.
3. **Budget from §2's rule**, then **cap concurrency by memory**: `workers ≈ (RAM_GB − 4) / 3`.
   For 31 GB that is 8, not 12.
4. **Split `:full` out.** It is 8× the `:none` tier and 15.4 h mean. Run it as its own batch so it
   cannot starve the cheap studies.

**Launcher design**

5. **Resume by content key, never by shard index** — it is what makes a kill or a restart cost one
   JIT instead of a shard.
6. **Mark exhausted queues.** Otherwise a lowest-free-id scheduler spins forever and starves the
   ids that still have work.
7. **Append new jobs at the END of a cost-sorted queue.** Inserting renumbers every later job and
   moves it to a different shard while workers are mid-study on it.
8. **Take the verdict from gate output, never exit codes.**
9. **`PHASEB_MAX_STUDIES`**: 1 bounds RSS drift but pays a full JIT per study — right for a 9 h
   study, wasteful for a 20 min one. Consider 2–3 once concurrency is memory-safe.

**Cluster**

10. **Validate every launcher against ≤48 cores/node and ≤224 GB/node before submitting.** A
    mechanical sweep (`ntasks == nodes × tasks-per-node == balfem_run N == PX·PY`) catches these in
    seconds and caught 11 broken launchers here.
11. **Never lower `mem-per-cpu` to get past a queue rejection** — 4 GB/core is measured-required.
    Lower the rank count instead.
12. **Set `BALFEM_OUTDIR` explicitly** whenever a case varies something the driver's tag does not
    carry (depth, duration, tolerance, geometry). Silent overwrite is otherwise certain.

**Instrumentation gaps worth closing**

13. Flag a **rising Newton count** automatically — it precedes divergence by hundreds of steps.
14. Record **per-Newton-iteration** linear counts, not per-stage min/max.
15. Split `t_solve` into assembly vs linear solve; a 55 s/step cost cannot currently be attributed.
