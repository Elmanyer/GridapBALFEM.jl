# `run/local/` — running GridapBALFEM on this workstation

The cluster loop is too slow to debug against: the four archived runs in `output/` burned 8 h,
44 h, 58 h and 72 h before failing. These launchers move the loop onto the local machine, where a
case finishes in **minutes** and can be watched in real time.

Everything here is sized for **6 cores**.

---

## Entry points

| | what it is |
|---|---|
| `balfem_local.sh` | the shared helper — resolves the project, caps the rank count, launches. `balfem_local_run <script.jl>` (sequential, **keeps point gauges**) and `balfem_local_mpi <n> <script.jl>` (MPI, `n ≤ 6`). Source it; don't run it. |
| `run_1d_<case>.sh` (7) | quasi-1D flume cases → `examples/local_1d/run_flume_1d.jl` |
| `run_2d_<case>.sh` (8) | small 2-D cases → `examples/local_2d/run_small_2d.jl` |
| `../../test/local/run_local_tests.sh` | the validation *gates* (sponge, relaxation zone, boundary modes, rest state) — pass/fail, not observation |

```bash
run/local/run_1d_lin_inner_flat.sh              # one 1-D case, sequential, ~8 min
run/local/run_2d_nl_full_plane_flat.sh          # one 2-D case, 6 MPI ranks
julia --project=. examples/inspect_run.jl output/local_2d/*   # judge the results
```

---

## The two things to know before reading a result

**1. Sequential runs have gauges; MPI runs do not.** The distributed driver evaluates no point
gauges — a point search would need inter-rank communication (`src/timeloop_dist.jl:21-22`). So:

* the **1-D** cases run **sequentially** by default (`BALFEM_MPI=0`) and return `diags` with
  `gauge_vals`, which is what phase- and amplitude-based analysis needs;
* the **2-D** cases run under **MPI** (`3×2 = 6` ranks) and are judged entirely from
  `<output_dir>/diagnostics.csv` — the machine-readable step log.

`BALFEM_MPI` flips either one, at the cost of the corresponding capability.

**2. The ~3.5 min JIT compile is paid once per process.** A case that runs for 4 minutes spends
nearly half its wall time compiling. Running several cases inside one Julia process amortises it;
running them as separate launchers does not. That is why the *tests* group several runs per file.

---

## What to look at

`diagnostics.csv` is written every `BALFEM_DIAG_EVERY` steps and carries the columns that the
archived cluster failures needed and did not have:

| column | question it answers |
|---|---|
| `eta_max`, `x_at_max` | how big, and **where** — a maximum pinned at a domain edge is the boundary mode, not a wave |
| `eta_max_int`, `eta_max_damped` | interior vs sponge/relaxation zone. **Read this as "≫1 means trouble", not "<1 means fine".** μ→0 at the sponge *entrance*, so once a wave train reaches the sponge it is at full amplitude there and the ratio sits naturally at ~0.6–0.8, oscillating with wave phase (measured on a healthy 2-D plane wave: 0.59–0.74). A boundary mode gives ~65. A ratio near 1 **combined with `x_at_max` pinned at a fixed x** is the thing to worry about — that pattern also appears if the interior source has been placed inside the sponge, which is a configuration error, not a solver fault |
| `u_max` | with `eta_max`, the kinematic ratio: a real wave gives `ω/tanh(kd)`; the η-dominated mode gives far less |
| `mass`, `mass_drift`, `energy` | the cheapest early warning of an instability |
| `lin_min`, `lin_max`, `lin_sat` | GMRES iteration range and whether it hit its cap. `lin_sat=1` means the solve was **truncated** and Newton is getting poor steps — the defect that hid for months |
| `nl_iters`, `nl_stages` | `nl_iters` is summed over the SDIRK stages; divide by `nl_stages` for the per-stage count |
| `rss_mb`, `rss_peak_mb` | per-rank memory. Flat is fine; a monotone climb is the OOM in progress |

`examples/inspect_run.jl` turns that file into a verdict; it is stdlib-only, so it also runs
against a cluster output directory. For fields, spectra and plots use
`postprocessing/GridapBALFEMPost`.

---

## Configuring a case

Both scripts take their whole configuration from `BALFEM_*` environment variables, so a launcher is
just the overrides that make it that case. The full list is in each script's header and in
`examples/distributed/_dist_common.jl`; the ones you change most often:

```bash
BALFEM_WAVE_GEN     1-D: inner|bc|sea          2-D: line|point|bc|sea
BALFEM_REGIME       linear | nonlinear
BALFEM_NL_PRESSURE  none | native | full
BALFEM_FLAT_BED     1 = flat bed | 0 = submerged bar
BALFEM_AWAVE        wave amplitude [m]         BALFEM_PERIODS  duration in wave periods
BALFEM_NX/BALFEM_NY   mesh                       BALFEM_DIAG_EVERY  diagnostics sampling
```

Two constraints the scripts enforce and will error on rather than let you discover later:
`BALFEM_NX % BALFEM_PX == 0` (and the same in y), and `BALFEM_NY ≥ 3` whenever the lateral BC is
`:periodic` — Gridap requires at least 3 elements in a periodic direction.

---

## Not the cluster

`run/balfem_env.sh` (the cluster helper) loads modules and resolves the prebuilt system image;
none of that exists locally, so `balfem_local.sh` deliberately shares no code with it. The cluster
counterpart of the 2-D cases is the suite in `run/dist_small/`; each 2-D launcher names its
cluster sibling in its header.

**1-D geometry default: `ny=1`, `y_wall_bc=:wall`, `Ly=dx`.** Gridap's 3-element minimum applies
to a *periodic* direction only; `:wall` and `:open` take `ny=1`. For a normal-incidence wave
`𝖴y ≡ 0` exactly, so a solid wall is consistent with the solution rather than an approximation to
it, while `:periodic` admits y-periodic modes a 1-D model does not have. Measured on the 240-cell
flume: **7215 free DOFs against 20202** for `ny=3`/`:periodic`.

**1-D cases are ALWAYS normal-incidence.** A 1-D domain has one propagation direction, so oblique
or short-crested content — transverse wavenumber `k_y = k sin θ` — cannot exist in it. Such a
request would be silently aliased onto a normal-incidence wave rather than rejected by the
mathematics, so `run_flume_1d.jl` refuses it explicitly: setting `BALFEM_WAVE_DIR` (non-zero),
`BALFEM_NTHETA`, `BALFEM_SPREAD_STD`, `BALFEM_THETA_MAX` or `BALFEM_DIRECTIONAL` is an error, not a
no-op. Directional seas belong in the 2-D driver. This is also *why* `:wall` is exact here rather
than an approximation: with normal incidence `𝖴y ≡ 0`, which is what the wall imposes.

**The 1-D cases have NO cluster counterpart, deliberately.** Cluster twins existed until
2026-09-02 and were deleted. A quasi-1D flume is a narrow y-periodic strip pinned at Gridap's
periodic minimum `ny=3`; direct-LU cost scales with the front width that cross-section fixes, so
the DOF count cancels from the serial-vs-distributed ratio and MPI buys nothing. Measured on this
workstation: 1 -> 2 ranks went 1.36 -> 10.18 s/step, i.e. NEGATIVE parallel efficiency. Run 1-D
here, sequentially; `run_all_1d.sh` runs the whole set side by side (`JOBS` processes at once),
which is the right kind of parallelism for cases this size.

**After editing `src/*.jl`, rebuild the cluster system image** (`compile/compile_snellius.sh`)
before submitting anything — the image bakes a compiled copy of the solver, and a stale one runs
old code. The launchers warn; `BALFEM_STRICT_SYSIMAGE=1` makes it an abort.
