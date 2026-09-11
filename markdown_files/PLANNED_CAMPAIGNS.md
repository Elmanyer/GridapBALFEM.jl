# PLANNED_CAMPAIGNS.md — the MMS convergence tests still owed

> **Scope — what belongs here.** *Everything designed but not yet run*: the next MMS convergence
> campaigns, their ladders, the code they need, and other unblocked runs. Answers **"what runs
> next, and what must be true before it does?"**
>
> **What does NOT belong here:** defects in completed work → [`OPEN_ISSUES.md`](OPEN_ISSUES.md);
> what is proven → [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md); the gate inventory →
> [`TEST_SUITE.md`](TEST_SUITE.md).


Written 2026-09-11, after Phase-B reached 32/34 valid studies.
Companion documents: `TEST_SUITE.md` (what exists), `COMPLETED_VBASIS_STUDY.md` (studies not begun),
`OPEN_ISSUES.md` (verification gaps). Data: `output/local/mms_phaseB/`.

---

## 0. THE DEFAULT LADDER — `nx = 4, 8, 16, 32` (+ 64 for Q2Q1 and Q3Q2)

| pairing | ladder | levels | why |
|---|---|---|---|
| `Q2Q1` | 4, 8, 16, 32, **64** | 5 | low order, errors stay far above the floor — the extra level is free information |
| `Q3Q2` | 4, 8, 16, 32, **64** | 5 | same |
| `Q4Q3` | 4, 8, 16, 32 | 4 | **stop at 32.** Fifth-order convergence puts `nx=64` near the algebraic floor, and `nx=128` on it (`e ≈ 8e-11`, measured) |

**Why it starts at 4.** A refinement study measures the rate of whichever error *dominates*, so the
fine end must stay above the algebraic floor (rule 32). Starting lower buys margin at the coarse
end, where the discretisation error is largest and cleanest. Four levels give **three pairwise
rates**, so a saturated finest level can be discarded and `4, 8, 16` still leaves a sequence —
which is what rule 33 says to read.

**Why Q4Q3 stops at 32.** Measured on the long ladder: at `nx=128` both `e_η` and `e_u` collapse to
~9e-11, and the reported `p_u` (1.08 / 2.51) is **unmeasurable, not defective**. Spending two more
levels there buys a number that cannot be read.

### The floor this is protecting against

Measured on the Q4/Q3 1-D long ladder (`T9_tier2_1d`, to `nx = 128`):

| | `e_η` | `e_u` |
|---|---|---|
| Q4/Q3, `nx = 128` | 8.96e-11 | 8.02e-11 |
| Q4/Q3, `nx = 32` (short ladder) | 2.29e-08 | 3.29e-09 |

Two independent quantities collapsing to the *same* ~9e-11 is the signature of a common floor, not
of convergence. The measured `p_u` there was 1.08 and 2.51 — **unmeasurable, not defective**.

### ✅ `nx = 32` IS CLEAN FOR Q4/Q3 — measured 2026-09-11, keep it

The suspicion that `e_u = 3.3e-09` might already be floor-contaminated is **refuted**. Extended
ladder, P1LFE-2 Q4/Q3 1-D, model 1, `nl_tol = 1e-14`:

| `nx` | `e_η` | ratio | `p_η` | `e_u` | `p_u` |
|---|---|---|---|---|---|
| 4  | 9.3412e-05 | --- | --- | 4.0567e-05 | --- |
| 8  | 5.8620e-06 | 15.9 | 3.994 | 2.0142e-06 | 4.332 |
| 16 | 3.6667e-07 | 16.0 | 3.999 | 7.9639e-08 | 4.661 |
| 32 | 2.2921e-08 | 16.0 | **4.000** | 3.2925e-09 | 4.596 |

`e_η` falls by exactly `16 = 2⁴` at **every** level including the last, and `p_η` converges to
**4.000**, the optimal order for Q4/Q3. A level sitting on an algebraic floor cannot produce an
exact fourth-order rate — saturation flattens the rate, it does not sharpen it. So `nx = 32` is
genuine discretisation error and **stays in the ladder**. `p_u` runs 4.33 → 4.66 → 4.60, still
climbing toward its optimal 5, i.e. pre-asymptotic rather than saturated (rule 33: a rising
sequence with the error still dropping).

The floor itself lies near **1e-10** — where the `nx = 128` long ladder put both errors — roughly
33x below the finest point used here. That is the margin the new ladder preserves.

⚠ The *mechanism* of that 1e-10 floor is still unidentified. `ls_rtol` is ruled out by construction
(the sequential path is a direct `LUSolver`; `ls_rtol` only binds the distributed GMRES path),
leaving `nl_tol = 1e-14` — at the edge of achievable, since `run_mms_matrix` died against a residual
of 1.4e-14 — and LU round-off amplified by conditioning. A confirming `nl_tol` sweep at fixed `nx`
was still running when this was written. It does not affect the conclusion above, which rests on the
rate sequence rather than on the floor's origin.

---

## 1. Q4/Q3 — re-run in 1-D now; 2-D deferred

| | status | why |
|---|---|---|
| **Q4/Q3, 1-D** | re-run | ran at `nx0 = 8` (`T9b`); needs `nx0 = 4` |
| **Q4/Q3, 2-D** | **never run** — deferred, not scheduled | too expensive for now; also blocked, §4 |

The 1-D short ladder (`T9b_tier2_1d_short`) already demonstrated that the re-specification works:
`p_η = 3.9997 / 4.0020` (optimal 4, exact) and `p_u = 4.60 / 4.30` climbing toward 5, against the
long ladder's meaningless 1.08 / 2.51. Re-running from `nx0 = 4` adds the coarse level and the
fall-back rate pair.

**2-D Q4/Q3 has no data at all, and remains the largest single gap** — but it is **deferred**
(2026-09-11): the 2-D matrix is 12 solver runs per model and too expensive to take on now. Recorded
here so it is not mistaken for covered. Its original specification carried the *long* ladder
(`levels = 4, nx0 = 8` → `nx` to 64), which would have saturated for the same reason 1-D did, at
roughly 20 h per study — its third level alone measured 3 h 20 m. **When it is picked up it must be
re-specified short; never launch it as it stands.**

---

## 2. THE CAMPAIGN TO RUN — 1-D only, four models, `:full` excluded

**`:full` is excluded by decision, not by omission.** Models 7 and 8 can never reach optimal order:
the `{1,2,4,5}` blocks are assembled in the residual but frozen (projected) in `jacobian_u`, so
Newton converges to the fixed point of a different map and the study measures the projection lag
rather than the discretisation. `T8_full_projection` already measured that floor (§3) and it grows
with `Nσ`. **Running convergence studies on `:full` is not useful and is not planned.**

**2-D is deferred.** The full matrix is 3 pairings × 4 ladder levels = 12 solver runs per model per
domain, which is too much to launch at once. 1-D first; a 2-D campaign follows once these land.

### The four models

Tier 2 (`T9`/`T9b`) swept only models **1 and 3** — both flat-bed, both `:none`. So every existing
pairing result is blind to bathymetry and to the `𝓝` blocks. By rule 4 a flat-bed study **cannot**
exercise `∇h` code, so the Q2/Q1 one-order velocity shortfall — the headline tier-2 result — is
currently established for `:none` physics on a flat bed only.

| model | `regime` / `bed` / `nl_pressure` | status |
|---|---|---|
| 1 | linear / flat / `:none` | ✅ done (Q2/Q1, Q3/Q2, Q4/Q3) |
| **2** | **linear / variable / `:none`** | ⛔ **run** |
| 3 | nonlinear / flat / `:none` | ✅ done (Q2/Q1, Q3/Q2, Q4/Q3) |
| **4** | **nonlinear / variable / `:none`** | ⛔ **run** |
| **5** | **nonlinear / flat / `:native`** | ⛔ **run** |
| **6** | **nonlinear / variable / `:native`** | ⛔ **run** |
| 7 | nonlinear / flat / `:full` | ✗ excluded — projection, see above |
| 8 | nonlinear / variable / `:full` | ✗ excluded — projection, see above |

Models 5 and 6 are `:native`, **the production tier**, so its pairing behaviour must be measured
rather than extrapolated from `:none`.

### Counting, and one correction

A **study** is one `run_conv_study` call: it runs the whole ladder internally and emits **one CSV
row**. So the unit that goes in the job queue is the study, not the individual solver run.

| item | models | studies | solver runs |
|---|---|---|---|
| **A. new campaigns** — `linear_var_none`, `nonlinear_var_none`, `nonlinear_flat_native`, `nonlinear_var_native` | 4 | 12 | 56 |
| **B. re-runs** — `linear_flat_none`, `nonlinear_flat_none` (per-level data was never recorded) | 2 | 6 | 28 |
| **TOTAL** | **6** | **18** | **84** |

Per model: `Q2Q1` 5 levels + `Q3Q2` 5 + `Q4Q3` 4 = **14 solver runs, 3 studies**.

### B. Why the two "already measured" models must be re-run

`linear_flat_none` and `nonlinear_flat_none` *were* measured by Phase-B tier 2 — but **only the
finest level's error was ever persisted.** `run_phaseB_shard.jl` calls `run_conv_study(verbose=false)`,
so per-level errors were never printed, and `run_mms_case` writes to `mktempdir()` with
`save_every=0`, so no diagnostics file exists either. The coarser errors were computed, consumed by
`convergence_rate`, and discarded at process exit.

⚠ **They are not recoverable.** Every file under `output/` was searched: 37 logs hold per-level MMS
lines, none from these studies. The nearest candidates (`conv_1d_seq.log`, `mms_matrix.log`) agree
on `e_η` to three figures but differ on `e_u` by a **factor of 10** at the same `nx` — a different
configuration, not the same study with a longer ladder.

⚠ **They must NOT be reconstructed from the recorded rate** (`e(nx/2) = e(nx)·2^p`). That assumes
the quantity being measured, so the fabricated points would "confirm" the rate they came from, and
it would erase the distinction between a pre-asymptotic and a saturated sequence.

Re-running is the only route and it is affordable: the old ladders reached `nx=128` and cost
**22.2 h** for eight studies; the new ones stop at 64 (32 for Q4Q3), removing the dominant levels.

**2-D is deferred — too expensive for now (decided 2026-09-11).** 2-D Q4/Q3 has never been run and
remains the largest single gap (§1), but it stays out of this campaign. When it is picked up it
must be re-specified short (`levels=4, nx0=4`) and never launched on the original `nx0=8`
specification.

⚠ For the record, the earlier count of "49" undercounted the Q4/Q3 item: a study is a whole ladder,
so the 1-D re-run alone is 2 studies = 8 solver runs. 48 + 8 = **56 solver runs across 14
studies**.

## 3. `:full` models — NOT to be re-run; the recorded floor, for reference

`T8_full_projection` ran models 7 and 8 at `a_eta = 0.4` on four bases. Results:

| basis | `p_η` | `p_u` | `e_u` floor |
|---|---|---|---|
| P1LFE-2 | 2.56 | 1.63 | 1.4e-06 |
| P1LFE-3 | 2.45 | 1.88 | 5.4e-06 |
| P1LFE-4 | 2.41 | 2.10 | 9.4e-06 |
| P2LFE-1 | 2.57 | 2.65 | 2.4e-07 |

Neither order is optimal, and **the floor grows with `Nσ`**. This is the known quasi-Newton gap:
the `{1,2,4,5}` blocks are in the residual but absent from `jacobian_u`, so Newton converges to the
fixed point of a different map, and the gap has a **cliff in amplitude** — it stalls at ~9.2e-04 at
`a_eta = 0.8` and only falls to ~2e-09 at `a_eta ≤ 0.4`.

⚠ **`:full` floors are not comparable across amplitudes.** Any `:full` number quoted must state its
`a_eta`, and must not be set beside one at a different amplitude.

**These numbers are kept as the record of the projection floor, not as a study to repeat.** The
values above already establish what a `:full` convergence study measures — the quasi-Newton
projection lag, not the discretisation order — which is precisely why §2 excludes models 7 and 8
from the campaign. Closing this properly means completing `jacobian_u`, not refining a mesh.

---

## 4. Blocked, and why — the two 2-D Q4/Q3 jobs

They are owned by **shards 6 and 7**, whose `exhausted_06` / `exhausted_07` marker files are dated
**Sep 2** and therefore predate the Sep 8 edit. The supervisor skips marked shards permanently, so
those jobs will never be launched.

**Cause, recorded because it generalises.** Removing P2LFE-2 from `T8` and `T10` shortened the
queue from **42 to 34 jobs**. The queue is cost-sorted and dealt round-robin, so deleting from the
*middle* **renumbered every later job and reassigned its owning shard** — exactly the hazard the
`T9b` comment in `run_phaseB_shard.jl` already warns about for *insertion*. The `exhausted` markers
then referred to the old numbering.

**Rule to adopt: after any change to the job list, clear `output/local/mms_phaseB/exhausted_*`.**
Those are empty flag files; deleting them costs nothing and the supervisor re-derives them. Better
still, neutralise cancelled jobs in place (or append) rather than deleting from the middle, so
ownership is preserved.

Since both blocked jobs need re-specifying anyway (§1), the marker clearing should happen as part
of that change, not before it.

---

## 4b. PREREQUISITE CODE CHANGES — diagnostics and output directories

⚠ **None of the runs below can satisfy the output requirements as the code stands.** The MMS driver
actively disables every diagnostic. In `run_mms_case`:

```julia
diags = run_time_loop(op, solver, u0, t0, T_final;
                      output_dir=output_dir, save_every=0,   # ← no VTK
                      print_every=typemax(Int),              # ← silent
                      diag_every=-1, check_every=0)          # ← diagnostics OFF
                      #  and no `monitor=` ⇒ no SolverMonitor
```

That is why every Phase-B log shows `solve time 00:00 (0.0% of wall)` and `Newton iters: 0` — those
counters are never incremented. Three changes are required, in this order.

### C0.1 — `run_mms_case`: expose and forward the diagnostic knobs

Add `save_every`, `diag_every`, `diag_csv`, `check_every` and `monitor` as keyword arguments
(defaulting to today's silent values, so nothing existing changes), and forward them to
`run_time_loop` instead of the hard-coded `0 / -1 / typemax`. Install a `SolverMonitor` when one is
requested, so `t_solve` and `nl_iters` become real.

### C0.2 — `run_conv_study`: a real per-level output directory

It already takes `output_dir` but passes `save_every=0` and one directory for the whole ladder. It
must instead build `<output_dir>/nx<N>/` per level and pass the diagnostic settings down, so each
level writes where the tree expects it:

```
output/local_1d/mms_convergence_campaigns/<model>/<pair>/nx<N>/
```

Also **set `verbose=true`** in `run_phaseB_shard.jl`: that alone would have preserved the per-level
errors this campaign has to re-measure (§2B).

### C0.3 — widen the CSV schema to carry the per-level sequence

`run_phaseB_shard.jl` currently persists `r.pw_eta[end]`, `r.fit_eta`, `r.e_eta[end]` — the last
pairwise rate, the fit, and the finest error. **The pairwise *sequence* is discarded**, and that is
the thing rules 32/33 say to read: a rising sequence with the error still dropping is
pre-asymptotic; a flat sequence with a stalled error is saturated; the fitted slope cannot tell them
apart. `run_conv_study` already returns `h`, `e_eta`, `e_u` as arrays, so this costs only columns.

## 4c. THE DIAGNOSTIC SNAPSHOT POLICY FOR THIS CAMPAIGN

Every run must record solver behaviour along the simulation, **first and last step included**:

| setting | value | why |
|---|---|---|
| `diag_csv` | `true` | writes `<nx>/diagnostics.csv` |
| `diag_every` | `max(1, nsteps ÷ 20)` → **5** at `nsteps=100` | ~20 samples per run: enough to see a trend, small enough not to dominate runtime |
| first/last step | **must both appear** | `run_time_loop` samples on `step % diag_every`, so step 1 and the final step need explicit inclusion — verify, do not assume |
| `save_every` | **0** | VTK is not wanted here: an MMS run is 100 steps on a tiny domain and the errors are the product. Revisit only if a run misbehaves and needs visual inspection |
| `monitor` | `SolverMonitor()` | makes `t_solve` and `nl_iters` real instead of zero |
| `verbose` | `true` | per-level errors printed to the shard log — the redundancy that would have prevented §2B |

⚠ **Verify the policy is live before spending 84 solver runs on it** (rule 38d): run one study,
confirm `<nx>/diagnostics.csv` exists for every level, that it contains the first and last step, and
that `nl_iters` is non-zero. A dead knob yields clean, confident, entirely wrong output.

## 5. Execution plan

### 5.1 The job list

Append to `examples/local_mms/run_phaseB_shard.jl` **at the end**, at the highest `cost`, so no
existing index — and therefore no existing shard ownership — moves (§4):

```julia
#  ---- CAMPAIGN C (2026-09-11) ------------------------------------------------
#  Ladder is PAIRING-DEPENDENT (PLANNED_CAMPAIGNS.md §0):
#      Q2Q1, Q3Q2 -> levels=5, nx0=4  ->  nx = 4, 8, 16, 32, 64
#      Q4Q3       -> levels=4, nx0=4  ->  nx = 4, 8, 16, 32      (64 saturates)
#  Diagnostics per §4c; per-level output_dir per §4b.C0.2.
levels_for(pu) = pu == 4 ? 4 : 5

#  C1 — the four models tier 2 never covered
for pu in (2,3,4), m in (2,4,5,6)
    push!(jobs, Job("C1_pairings_1d","P1LFE-2",2,1,m,pu,:d1,levels_for(pu),4,0.8, 96))
end

#  C2 — RE-RUN the two models tier 2 did cover. Their per-level errors were never
#       recorded (verbose=false + save_every=0) and are NOT recoverable; see §2B.
for pu in (2,3,4), m in (1,3)
    push!(jobs, Job("C2_rerun_1d","P1LFE-2",2,1,m,pu,:d1,levels_for(pu),4,0.8, 97))
end

#  ⚠ NO 2-D (decided 2026-09-11): too expensive for now. 2-D Q4/Q3 has never been
#    run and is still the largest gap; when taken up it must use nx0=4, NEVER the
#    original nx0=8 spec (~20 h/study, and it saturates).
```

**18 studies, 84 solver runs.** `cost` orders C1 before C2, so the genuinely new physics lands
before the re-measurement of models whose *rates* are already known.

### 5.1b Output directories

Each study writes into the tree it belongs to, one directory per level:

```
output/local_1d/mms_convergence_campaigns/<model>/<pair>/nx<N>/
    config.json          # written by the driver, not reconstructed afterwards
    diagnostics.csv      # the per-step snapshots of §4c
    result.csv           # e_eta, e_u at this level
```

with `<model>` ∈ {`linear_flat_none`, `linear_var_none`, `nonlinear_flat_none`,
`nonlinear_var_none`, `nonlinear_flat_native`, `nonlinear_var_native`} and `<pair>` ∈
{`Q2Q1`, `Q3Q2`, `Q4Q3`}. The mapping from `(model, p_u)` to that path belongs in the shard, beside
the job definition, so a job carries its destination rather than having it inferred later.

⚠ The existing `nx8 … nx128` directories under `linear_flat_none` and `nonlinear_flat_none` are
from the **old** `nx0=8` ladders. Leave them: they hold each level's reconstructed `config.json`
plus a `NO_RESULT.md`, and they are the record that those levels ran on a ladder no longer used.
The new runs land in `nx4 … nx64`, which overlaps at 8/16/32 — so **suffix the new files by task**
(`result__C2_rerun_1d.csv`), exactly as the `Q4Q3` levels already do for T9 vs T9b.

### 5.2 Before launching — three things, in this order

0. **Make the §4b code changes and verify them on one study** (§4c) — per-level `output_dir`,
   diagnostics on, `verbose=true`, widened CSV. Everything else is wasted without this.
1. **Clear the stale ownership markers.** `rm output/local/mms_phaseB/exhausted_*`. They are empty
   flag files; the supervisor re-derives them. Without this the appended jobs may land on a shard
   that is permanently skipped (§4).
2. **Restart the shards.** They read the job list at startup, so workers already running hold the
   old queue. Kill them and let the supervisor relaunch; `PHASEB_MAX_STUDIES=1` plus the resume
   logic makes a restart cost one JIT.
3. **Verify the queue before committing hours to it.** Run the shard enumeration and confirm the
   new jobs appear, are owned by live shards, and that `levels=4, nx0=4` is what the driver
   receives. Rule 38d: verify a new knob is live before spending the run on it.

### 5.3 Order of value, if the batch has to be cut short

1. **C1 model 2** (linear, variable bed) at all three pairings — the cheapest way to put `∇h`
   under the pairing question, and linear models give the cleanest rates.
2. **C1 models 5, 6** (`:native`) — the production tier.
3. **C1 model 4** (nonlinear, variable bed).
4. **C2** (the two re-runs) — these recover *per-level* data for models whose fitted rates are
   already known, so they add resolution rather than new physics. Lowest value of the four, but
   still owed: without them half the tree has one data point per study.

⚠ **Do §4b first regardless.** Running C1 before the driver writes per-level output would reproduce
exactly the gap C2 exists to repair.

### 5.4 Cost

Measured from the 38 Phase-B studies: **310.8 core-hours**, i.e. ~8 h/study averaged, but the
spread is wide. This campaign is **1-D only**, which is the cheap end: comparable 1-D studies in
Phase-B ran from minutes to ~1 h. Run under `supervise_phaseB2.sh`, which bounds worker lifetime
(rule 41) and watches swap.

⚠ **Read the pairwise sequence, not the fitted slope** (rule 33), and check `e_u` magnitude against
the floor (§0) before trusting any Q4/Q3 rate.

## 6. Housekeeping owed alongside

* **`ny0` for 2-D studies** — `run_conv_study` defaults `ny_1d = 1` for `:d1` (correct: the
  manufactured field has `k_y = 0`, so the solution is y-invariant). Confirm the `:d2` ladder
  refines `ny` with `nx`, or the 2-D rates measure a 1-D refinement.
* **`pw_eta` / `pw_u` are written to the CSV as a single fitted number**, so the pairwise sequence
  — the thing rules 32/33 say to read — is recoverable only from the shard logs. Widen the CSV to
  carry the per-level errors, or the next campaign repeats this analysis by hand.
* **`solve time 0.0%` / `Newton iters: 0`** in every Phase-B log are instrumentation gaps, not
  stalls: `mms_driver.jl` installs no `SolverMonitor`. Either install one or stop printing the
  fields, because they currently read as a hang.

---

## 7. Other runs designed, unblocked, not yet executed

* **Re-run the small-domain suite with the current fixes in place.** Every cluster run on record
  predates at least one of: the surface-damping sponge, the GMRES configuration fix, and the
  `mu_max` 8→40 raise. **The archived outputs are not a baseline** — the one that NaN'd at t=13.8 is
  a diagnosed pre-fix failure, not a solver defect. What to check on the re-run: `gmres=` no longer
  pinned (now flagged automatically), Newton ~4/step, and — straight off `diagnostics.csv` —
  `x_at_max` staying interior with `eta_max_damped/eta_max_int < 1`. **Rebuild the sysimage first.**
  Cheap local proxies exist and should be run first: `run/local/run_2d_*.sh` are scaled-down
  siblings of eight of these cases.
* `Hs = 0.2 m` at `d = 3.5 m` with `nl_pressure=:full` is far outside the conservative
  `A ≤ 0.001` guidance and is the most aggressive case in the suite. If it still destabilises,
  the sponge/relaxation widths for the longest components are the next thing to size.
* **At-scale physical benchmarks** (Stokes harmonics, Dingemans bar) — scripts exist; quantitative
  overlays on the published data are pending.
* **Production-length (200+ `Tp`) irregular/directional sea runs** — scripts ready and unblocked;
  not yet run at length.

*(Moved here from `OPEN_ISSUES.md` §8 on 2026-09-11: these are plans, not defects.)*

