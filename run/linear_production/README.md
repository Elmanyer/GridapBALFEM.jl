# `run/linear_production/` — linear-model production runs on Snellius (`rome`)

Seven SLURM launchers for the **linear** model: model 1 (flat bed) and model 2 (variable bed).
They run v2 (branch `v2-solver`). The linear residual and both Jacobians are identical to v1 entry by
entry (`test/v2_migration/regression_snapshot.jl`: `lin_flat`, `lin_bed`, `lin_bed_q21`, relative
difference 0.00e+00).

Domains are those of the v1 `run/dist_small/` suite. The directional sea uses the **wide_small**
domain (50 × 40 m, 2 m lateral sponges).

## Cases

| launcher | wave | bed | domain [m] | mesh (dx) | ranks · nodes | `--mem` / node | billed | duration |
|---|---|---|---|---|---|---|---|---|
| `run_lin_bcplane_flat.sh` | BC plane, A 0.1, T 2.0 | flat | 50 × 20, y-periodic | 200 × 40 (0.25) | 28 = 7×4 · 1 | 112G | 4/8 | 40 T = 80 s |
| `run_lin_bcplane_bar.sh` | BC plane, A 0.1, T 2.0 | bar | 50 × 20, y-periodic | 200 × 40 (0.25) | 28 = 7×4 · 1 | 112G | 4/8 | 40 T = 80 s |
| `run_lin_irregular_flat.sh` | JONSWAP, Hs 0.2, Tp 2.0 | flat | 50 × 20, y-walls | 200 × 40 (0.25) | 28 = 7×4 · 1 | 112G | 4/8 | 60 Tp = 120 s |
| `run_lin_irregular_bar.sh` | JONSWAP, Hs 0.2, Tp 2.0 | bar | 50 × 20, y-walls | 200 × 40 (0.25) | 28 = 7×4 · 1 | 112G | 4/8 | 60 Tp = 120 s |
| `run_lin_ring_flat.sh` | point source, A 0.1, T 2.0 | flat | 40 × 40, 4 sponges | 160 × 160 (0.25) | 42 = 7×6 · 1 | 168G | 6/8 | 30 T = 60 s |
| `run_lin_directional_flat.sh` | JONSWAP × cos-power, 7 θ, ±40° | flat | 50 × 40 wide_small | 200 × 160 (0.25) | 42 = 7×6 · 1 | 168G | 6/8 | 50 Tp = 100 s |
| `run_lin_directional_bar.sh` | JONSWAP × cos-power, 7 θ, ±40° | bar | 50 × 40 wide_small | 200 × 160 (0.25) | 42 = 7×6 · 1 | 168G | 6/8 | 50 Tp = 100 s |

**The bar** is the same in every variable-bed case:
* `d = 3.5`, height 1.5 (crest depth 2.0);
* centred at `x = 26`, half-width 6, `tanh` ramps 2 m wide;
* y-invariant, i.e. shore-parallel for the directional sea.

The irregular and directional pairs share `BALFEM_SEED`, so the flat and bar runs are driven by the
same sea.

## Shared settings (`_linear_common.sh`)

* `BALFEM_REGIME=linear`, `BALFEM_NL_PRESSURE=0` (v2 Boolean).
* P1LFE-2 vertical, Q2/Q1 Taylor–Hood, `dt = 0.02`.
* **Crank–Nicolson** (`BALFEM_SOLVER=theta`). In v1 these cases ran on SDIRK_2_2, which removes
  about 20 % of the kd≈6 tail per domain transit (rule 15) and costs two stages per step. The linear
  model has no instability for its dissipation to mask. To compare against v1, set
  `BALFEM_SOLVER=sdirk`.
* `BALFEM_STRICT_SYSIMAGE=1`: a stale sysimage aborts the job instead of silently running old code.

## Memory: `--mem` only

* Every launcher requests memory with `--mem=<per-node total>` and **never** with `--mem-per-cpu`.
* Slurm counts `--mem-per-cpu` against all *allocated* CPUs, and rounds that allocation up to
  16 cores. So `28 × 4G` became a request for `32 × 4 = 128 GiB`, which was billed 5/8 instead of
  4/8, or refused outright.
* With `--mem`, a tier is `k/8 = 28k GiB` for `7k` ranks at 4 GiB/rank.
* `--mem` is **per node**. Every job here is single-node.
* The directional pair runs on **one node, 1 × 42 at 168G**, not v1's 2 × 28 at 112G. The per-rank
  budget is the same, the bill is 6/8 instead of 8/8, and there is no inter-node traffic.
  v1's measured wide run (the 2.84 GB and 73.5 s/step below) was 1 node × 40 ranks on the same
  200 × 160 mesh, i.e. 800 cells/rank against 762 here. If it is OOM-killed, use 1 × 56 with
  `--exclusive --mem=0`, `PX=8 PY=7`.
* The measured billing and the full rule are in
  [`../SNELLIUS_ROME_LAUNCH_CONFIGS.md`](../SNELLIUS_ROME_LAUNCH_CONFIGS.md) §2 and §4.

## Cost estimates

The v1 linear runs on the same meshes and rank grids were measured under SDIRK_2_2. Under CN, expect
roughly 0.55–0.65 of those wall times.

| case | v1 s/step (SDIRK) | steps | worst case (SDIRK rate) | expected (CN) | v1 peak RSS/rank |
|---|---|---|---|---|---|
| bcplane flat / bar | 23.9 / 26.6 | 4000 | 27 / 30 h | ~15 / 17 h | 2.6 / 2.5 GB |
| irregular flat / bar | ~24–27 | 6000 | ~45 h | ~25–27 h | ~2.6 GB |
| ring | not measured (linear) | 3000 | uncertain, ≲ 90 h | ~40–50 h | ~2.9 GB (est.) |
| directional flat / bar | 73.5 (flat, v1 at 40 ranks) | 5000 | 102 / ~112 h | ~55–65 h | 2.84 GB |

The directional bar at SDIRK rates would come close to the 119:59 wall. That margin is a second
reason for CN.

## Launch

```bash
# on Snellius, after pulling v2-solver:
cd $HOME/GridapBALFEM.jl && git fetch && git switch v2-solver && git pull
cd compile && sbatch compile_snellius.sh                 # rebuild the sysimage (rule 25)
# when it has finished:
bash $HOME/GridapBALFEM.jl/run/linear_production/submit_all.sh            # all seven
bash $HOME/GridapBALFEM.jl/run/linear_production/submit_all.sh directional  # a subset
```

Outputs go to `output/<output_dir_name>/`, one per case, with a `run_manifest.toml` that records
the integrator and every knob. `unique_output_dir` never overwrites an existing directory.

## What to watch

* `diagnostics.csv`: `eta_max`, `x_at_max`, `rss_peak_mb` and `lin_last` (GMRES iterations).
* Directional: `eta_max_damped / eta_max_int`. A ratio above 1 with `x_at_max` ≈ 2–3 m is the inflow
  relaxation zone and is benign. If the maximum is inside the lateral sponges, try
  `BALFEM_MUMAX=16`.
* Variable bed: the lee-shoulder growth (`markdown_files/OPEN_ISSUES.md` §3) was observed in
  nonlinear runs only. These are the first long linear bar runs in 2-D, so read the bar's downwave
  shoulder before calling them clean.
