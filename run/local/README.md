# `run/local/` — running GridapBALFEM on this workstation (v2)

`balfem_local.sh` is the shared helper: it resolves the project, caps the rank count, gates the
Taylor–Hood pairing and launches. Source it; don't run it.

* `balfem_local_run <script.jl>` — sequential (keeps point gauges);
* `balfem_local_mpi <n> <script.jl>` — MPI.

A launcher is a few `export BALFEM_*=…` lines followed by `balfem_local_run`. ⚠ Overrides placed
AFTER the run line are no-ops (rule 38h), and a launcher must never be edited while bash executes
it. Write logs to `output/…`, never to `/tmp` (rule 47).

**v2 (2026-10-05).** The 122 v1 launchers that lived here — the equal-order ladder, the
Taylor–Hood campaign, the 1-D production and refinement factorials, the mixed/projected flume
factorial, the periodic-box batches, the broken/ghost-penalty campaign — were removed from the v2
tree. They ran configurations that no longer exist (`:native`, `:full` projected or mixed, the
component mask). They remain in tag `v1_final_solver`, and their outputs are in `output_v1/`.
v2 launchers are written as the v2 campaign needs them (`markdown_files/V2_SOLVER_PLAN.md` §4).

v2 knobs that changed: `BALFEM_NL_PRESSURE` takes `0`/`1` (all eight 𝓝 components off/on).
`BALFEM_MIXED`, `BALFEM_P_AUX`, `BALFEM_C3_MASK`, `BALFEM_NLP_INLOOP` and `BALFEM_BROKEN` are
removed, and every driver REFUSES them if set (`check_v1_env`).
