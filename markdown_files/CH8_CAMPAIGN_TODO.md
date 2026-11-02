# Chapter 8 campaign: the working to-do list (opened 2026-10-31)

The ordered path from the current state to a running 1-D stability campaign. Work top to bottom.
When a task finishes, tick it, add the date and a one-line result, and continue with the next open
task. A task marked **ASK** waits for the author. Do not skip it, and while it waits, keep working
on later tasks that do not depend on it.

Standing constraints:
* No commits or pushes unless the author asks in that message.
* Stop and report on any design error in existing code.
* Runs longer than 30 min go detached, logging into `output/`.
* At most about 6 Julia processes at once; give each one its own systemd scope (`workstation-ram-budget`).
* Any `src/` edit makes the sysimage stale (rule 25).

---

## 0. Already in flight (check these first)

- [x] (2026-10-31; P1 measured with the profiler on; P2 crashed because the driver had gained the
  checkpoint kwargs mid-run, which is harmless) **Profile of model 5 in the 8-cell box.**
  * Scope `prof-m5-box2`, log `output/profile_m5_box/profile.log`.
  * Outputs: `profile_flat.txt`, `profile_tree.txt`, and the `>>> SUMMARY` line (P1 with VTK and
    diagnostics vs P2 without).
- [x] (2026-10-31) PASS: 2 steps converge (2 Newton iterations each), bed d ∈ [3.437, 3.563] m, slope
  0.100, 480 free DOFs, about 29 s/step (the same assembly overhead as model 5). **Smoke test of model 6 with the cosine bed, 2 steps.**
  * Scope `stab1d-smoke`, log `output/1D_stability_campaign/_smoke/m6_n8.log`.
  * Pass: both steps converge, the `BED` banner line is printed, and there is no error.

## 1. The assembly overhead (from the profile)

Evidence: about 44 s of assembly per Newton iteration for 14 k DOFs (`krylov_conditioning/mpi_check`),
and 31 s/step in the 8-cell box. The step time is nearly independent of the cell count, so the cost
is fixed, not proportional to the problem size.

- [x] (2026-10-31) **Finding.** About 54 % of all samples are spent BUILDING Gridap's lazy-array
  caches (`array_cache`), not evaluating anything. That cost is per assembled term, not per cell,
  which is why the step time barely depends on the mesh. It has two parts:
  * ~30 %: `Gridap.Arrays._get_cache` (LazyArrays.jl:155). When a node's objectid is not in the
    cache dict, it scans the WHOLE dict with a structural `same_branch` comparison. The cost is
    O(nodes²) per term, and the node that keeps missing is `∇(d_cf)` of the depth CellField:
    `alg_dx(d_cf)` creates a new node at every call site, even on a flat bed (`0.0*alg_dx(d_cf)`).
  * ~20 %: Gridap's debug-mode `@check` (Helpers/Macros.jl:61). This is the default execution mode;
    `Gridap.Helpers.set_performance_mode()` removes it.
  Experiment `output/profile_m5_box/cache_experiment.jl` (no repository change): the same run with
  `_get_cache` reduced to the objectid lookup. **RESULT (2026-10-31):** 35.9 → 25.5 s/step
  (−29 %; 18.0 → 12.7 s per Newton iteration), and the answer is **bit-identical** (final max η
  equal to all 16 digits). The machine was loaded by the campaign's compilation, so the absolute
  numbers are noisy; the ratio was measured in one process under the same load. Next, while the
  campaign runs: Gridap performance mode in a separate environment copy (the @check share), and
  the number of DAG nodes per term.
  Read the profile and name where the time goes: residual, `jacobian_u` or `jacobian_u_t`;
  volume terms or skeleton terms; repeated compilation or construction; VTK and diagnostics.
- [ ] If the cause is fixable without changing validated maths:
  * fix it;
  * re-run the profile to measure the speed-up;
  * re-run `test_jacobians_ad`, `test_broken_formulation` and the regression snapshot. The answer must
    not change.
- [ ] If the fix touches validated code or its design: **ASK** before changing anything.
  The options for the author:
  * (a) Gridap performance mode, a preference with no code change;
  * (b) build the depth CellFields (d, ∇h, ∇∇h, and constant zeros on a flat bed) once per problem,
    in `problem.jl`/`broken.jl`/`nlpressure.jl`, then verify with the regression snapshot and
    `test_jacobians_ad`;
  * (c) make `_get_cache` objectid-only in the Gridap fork, a library change.

## 2. Checkpoint and resume (u_n only; the author chose this on 2026-10-31)

Design (agreed):
* A checkpoint holds t, the step number, the free DOF values of u_n, the git commit and a
  configuration fingerprint.
* Resuming rebuilds the problem, refuses a fingerprint mismatch, and continues from t0 to `T_final`.
* It does **not** store Gridap's Newton-guess buffers (u_{n−1} for CN, the stage slopes for SDIRK).
  A resumed run is therefore a valid continuation but not bit-identical: it differs at the Newton
  tolerance level. Measure that difference; never assume it.

- [x] **Format** (2026-10-31): raw little-endian Float64 plus a hand-written TOML header, the
  same no-new-dependency style as `run_manifest.toml`. Files: `checkpoints/ck_<step>.{toml,f64}`,
  or `.part<k>.f64` per rank. The header is written last and atomically, and it stores Σv and Σv².
- [x] (2026-10-31, `src/checkpoint.jl`, both time loops) **Writing (sequential).** `checkpoint_every` kwarg in `setup_and_run` and `run_time_loop`:
  * write `checkpoint_<step>.jld2` every N steps, at the end of the run, and at divergence;
  * keep only the latest K, to bound disk use.
- [x] (2026-10-31: `CHECKPOINT_HARD_KEYS`/`_SOFT_KEYS`, the free DOF count and a 7×3 bathymetry sample) **Fingerprint.** Every quantity that fixes the discrete problem:
  * M, p_vertical, c_bdy, p_u, p_eta, domain, partition, periodicity and BCs;
  * regime, nl_pressure, flat_bed and a bathymetry sample;
  * the CIP γ/order/hexp and the stabilization;
  * dt, solver_type, tableau and θ;
  * the free DOF count.
  `dt`, the integrator, `nl_tol` and the output cadence are allowed to change only with
  `restart_allow_change = true`, and the manifest then records that change.
- [x] (2026-10-31; mass0/energy0 are carried over, so diagnostics stay relative to the original t = 0) **Resuming (sequential).** `restart_from = <checkpoint>`:
  * build everything as usual, verify the fingerprint, and set u0 from the stored free values;
  * t0 = t_checkpoint, and the step counter continues from there;
  * forcing (source, BC generation, ramp) is evaluated at absolute t, so check that the time loop
    receives t0 everywhere.
- [ ] **Output continuity.**
  * Write to a new output directory (`unique_output_dir`, never overwrite).
  * The manifest gets a `[restart]` section: parent directory, checkpoint, step, t, and the commit
    of both the parent and this run.
  * `diagnostics.csv` and `band_energy.csv` continue in absolute t.
  * Teach `GridapBALFEMPost` / `periodic_growth.jl` to follow a chain of run directories. **(open)**
- [x] (2026-10-31; the part count is checked before MPI starts) **Distributed.** Each rank writes its own values plus the partition description. Resuming
  requires the same `cpu_grid` and refuses anything else. Re-partitioning is out of scope; write
  it down in OPEN_ISSUES.
- [x] (2026-10-31: `checkpoint_kw()` in `_dist_common.jl`, plus a `resume` name token) **Drivers.** `BALFEM_CHECKPOINT_EVERY` and `BALFEM_RESTART_FROM` in `_dist_common.jl`,
  `run_flume_1d.jl` and `run_periodic_1d.jl`.
- [x] **Gates:** `output/restart_gate/test_restart.log` (sequential, R1–R8): **16/16 PASS** (2026-10-31).
  Resumed vs straight at step 10: m3 CN 1.4e-8, m3 SDIRK 3.4e-6, m5 CN 5.5e-8. Every reloaded
  state is bit-identical, and all four refusals fire. The remaining open items: the runtests tier entry
  and the postprocessing chain.
  `test_restart_distributed.log` (2 MPI ranks, D1–D4): **4/4 PASS** (2026-10-31). The resumed run
  differs from the straight run by a relative 4.6e-6 at step 10 (model 3, forced flume, CN, `:schwarz`,
  nl_tol 1e-8 and ls_rtol 1e-5, so the GMRES tolerance sets the level). Then add `test_restart.jl` to the
  `runtests.jl` tier, and the MPI one to its list of MPI tests.
  **Gate `test/test_restart.jl`.** Cases: model 3 and model 5, in the box, under CN and SDIRK;
  sequential, and distributed on the debug backend.
  * One run of N steps straight through, against N/2 steps + checkpoint + resume for N/2 steps.
  * Report max|Δu| and |Δη|/|η| at the end; the expected size is the Newton-tolerance level.
  * Confirm that a resumed run with an unchanged state reproduces the time and step counters.
  * The refusals fire on a mismatched fingerprint and on a mismatched `cpu_grid`.
  * Verdict from the gate output, never from the exit code (rule 35).
- [ ] **Optional.** `restart_from_vtu` for old runs. Their sub-cell VTK holds every Lagrange
  node in Float64. Only if the gate shows that u_n-only resumption is acceptable.
- [ ] Docs: `README.md` (kwargs), `RUNNING.md` (how to resume), `CLAUDE.md` §3 features and §5.

## 3. Housekeeping before the rebuild

- [ ] Fast suite `test/runtests.jl`. It has not been run since the 2026-10-30 changes; read the
  verdict from the gate output.
- [ ] Re-run `kc4_ownblocklu.jl` if `timeloop_dist.jl` changed in task 2.
- [ ] **Docs for the solver changes of 2026-10-31:**
  * `CONFIGURATION.md` §5: `:schwarz` now works as `OwnBlockLUSolver`, with 15–19 GMRES iterations on
    the production case against Jacobi's cap. `:gs` is dropped.
  * `OPEN_ISSUES.md`: the preconditioner item.
  * `CLAUDE.md`: §3, rule 16, and the integrator list (`:gen_alpha` and `:rk3` removed).
  * `STATUS.md`.
  * `STABILITY.md` §3: model-3 Q3/Q2 n32 was **accelerating** at 160 s (+0.10–0.12 s⁻¹ over the last
    60 s), not just "growing +0.07". n16 was growing slowly (+0.02–0.03).
  * Rule 24 caveat: at Q3/Q2, Jacobi at `ls_rtol = 1e-5` leaves a 4–5 % solution error
    (`kc4_ownblocklu.log`).
- [ ] **ASK: distributed default preconditioner.** Keep `:jacobi` as the default, or switch to
  `:schwarz`? `:schwarz` is measured better on 1-D; on 2-D at more ranks it is unmeasured.
- [ ] **ASK: commit before the campaign?** The run manifests otherwise record a dirty tree (rule 2c
  reproducibility). The pending set:
  * the integrator removal;
  * `OwnBlockLUSolver` and the removal of `:gs`;
  * the periodic cosine bed;
  * checkpoint and resume;
  * `compile/compile_local.{jl,sh}` and the campaign launcher.

## 4. Local sysimage (optional, but it saves ~25 min of compilation per campaign case)

**Skipped for the campaign (2026-10-31, author: "prioritise launching").** `src/` changed that day,
and a rebuild costs about 4 h, more than the ~25 min of JIT each campaign slot pays. Rebuild after
the campaign, when `src/` has settled.

- [ ] Only once `src/` is final (tasks 1–3), and with nothing else heavy running.
- [ ] `compile/warmup.jl` dry run: no `gen_alpha`/`rk3`/`gs` cases; the `:schwarz` case is still there.
  Consider adding a model-6 x-periodic box case with the cosine bed, so the campaign path is traced.
- [ ] Build:
  `systemd-run --user --scope --unit=balfem-sysimage -p MemoryMax=24G nohup compile/compile_local.sh > output/sysimage_build/build.log 2>&1 &`
  * The 2026-10-30 build was OOM-killed at 16 G in the final compile, after 3 h 40 min of warmup.
  * If 24 G is not enough, reduce the warmup to the campaign path (sequential models 3, 5 and 6, CN, box).
- [ ] Pass: the warmup summary shows **0 FAILED**, `GridapBALFEM_sysimage_local.so` exists, and
  the `.src.sha256` stamp matches.
- [ ] Smoke: one 2-step box case with `-J`, checking that the first step takes seconds, not about
  25 min.

## 5. The stability campaign

- [x] **LAUNCHED 2026-10-31 14:29.**
  * Scope `stab1d-campaign`, logs in `output/1D_stability_campaign/_logs/`.
  * Source snapshot in `_logs/source_20261031_142905/` (commit 1ae2695 + tracked diff + untracked
    tarball), because the tree is uncommitted.
  * No sysimage; checkpoints every 250 steps, the 2 newest kept.
  * Resume a lost case with `BALFEM_RESTART_FROM=output/1D_stability_campaign/<case>`, deleting
    its log first (the launcher skips existing logs).
- [ ] Launch:
  `nohup run/local/run_1d_stability_campaign.sh 6 > output/1D_stability_campaign/_logs/_batch.out 2>&1 &`
  * 11 cases: m5 and m6 at n64/32/16/8, and m3 at n32/16/8.
  * Q3/Q2, CN, A = 0.10, 200 periods, seed 1e-8.
  * m6 has a cosine bed with slope 1:10.
  * Check `_batch.log`, which says whether the sysimage was used.
- [x] Checkpoints are on (launcher `CHECKPOINT_EVERY=250`).
- [ ] Monitor:
  * one status line per case (t, max η, Newton iterations, s/step, RSS);
  * divergences, with t and x_at_max;
  * memory against the 5 G per-case cap.
- [ ] Per run, when it ends:
  `PG_NCELL=<n> PG_PU=3 julia --project=postprocessing postprocessing/examples/periodic_growth.jl <dir>`.
  Report a growth rate per band, both early and late, read together with x_at_max (rules 38f and 40).
- [ ] Compare against `output/1D_stability_campaign/model3_reference/` and the m3 reruns, then
  write the results into `STABILITY.md` and report them to the author, for LaTeX chapter 8.

## 6. Parked (waiting for a decision)

- [ ] Re-launching the 1-D flume CIP-monitor runs (m5/m6, A = 0.02/0.05). Sequential LU, or
  4 ranks with `BALFEM_PRECOND=schwarz`, now that it works.
- [ ] 2-D with 𝓝: `:schwarz` at more ranks, which has no coarse space; check its iteration growth
  with the rank count.
- [ ] γ calibration of the CIP (chapter 8 route, step 2).
