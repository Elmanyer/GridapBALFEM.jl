# `run/integrator_study/` — the time-integrator stability study (Snellius `rome`)

**Question** (STABILITY.md §6b, author 2026-11-11): Crank–Nicolson has shown the discrete nonlinear
operator unstable in every box run. How much of that growth does numerical dissipation from the time
integration hold, and at what cost to the carrier? The study comes before any new weak-form term.
It leans on the regime that challenges the full nonlinear solver: steep waves (kA up to 0.32).

## The set-up (every case)

* **Box:** the closed x-periodic box (`examples/local_1d/run_periodic_1d.jl`), one model wavelength,
  flat bed (m6: cosine bed, max slope 0.1), seed 1e-8 on every harmonic.
* **Discretisation:** P1LFE-2 (Nσ = 3), Q3/Q2, `ncell` cells per wavelength. CIP is the driver
  default: γ = 0.3 with 𝓝, none without.
* **Time:** 200 periods, dt = T/40. ωΔt = 0.157 at every T, so each integrator removes the same
  share of the carrier per period whatever the period.
* **Steepness:** `BALFEM_KA` fixes k₀A with k₀ the model's own wavenumber. The amplitude follows:

| T [s] | kd (model) | 2ω/ω∞ | A at kA = 0.158 / 0.25 / 0.32 [m] | H/λ vs Miche | H/d at kA 0.32 |
|---|---|---|---|---|---|
| 1.6 | 5.54 | 1.00 | 0.100 / 0.158 / 0.202 | 35 / 56 / 72 % | 0.12 |
| 2.2 | 2.94 | 0.73 | 0.188 / 0.298 / 0.381 | 36 / 56 / 72 % | 0.22 |
| 3.0 | 1.68 | 0.53 | 0.329 / 0.521 / 0.667 | 38 / 60 / 77 % | 0.38 (Ursell 5.3) |

  All within the P1LFE-2 ranges (kd ≤ 10.84 for C, ≤ 6.42 for Cg). The kd 5.5 cases sit at the
  triad threshold (2ω = ω∞); kd 3.0 and 1.7 are below it.

## Integrators (`cases.py` INTEGRATORS)

Five in the matrix (author, 2026-11-11: TR-BDF2, Crouzeix 3/4 and generalised-α dropped, explicit
RK3 added):

| tag | scheme | order / stages | loss per step: carrier (ωΔt 0.157) / 2ω₀ / ωΔt = 1 | role |
|---|---|---|---|---|
| `cn` | θ-method θ = 0.5 | 2 / 1 | 0 / 0 / 0 | **diagnosis pair** (rule 15): what grows with no damping |
| `sdirk22` | `SDIRK_2_2` (DIRK22(1,0,1)) | 2 / 2, R(∞) = −½ | 2.2e-4 / 3.0e-3 / 9.9e-2 | strong damping at every frequency; the carrier pays (58 % per 100 periods) |
| `sdirk43` | `SDIRK_4_3` | 3 / 4, **L-stable** | 1.3e-5 / 2.0e-4 / 1.7e-2 | the only L-stable scheme: removes stiff modes outright, carrier cost 5 % per 100 periods |
| `rk4` | `EXRK_RungeKutta_4_4` | 4 / 4, explicit | 1.0e-7 / 6.6e-6 / 6.1e-3 | Yang & Liu's scheme: nearly non-dissipative; unstable if a discrete mode has ωΔt > 2.83 |
| `rk3` | `EXRK_SSP_3_3` | 3 / 3, explicit | 2.5e-5 / 3.9e-4 / 2.8e-2 | the dissipative explicit scheme: SDIRK_2_2's (ωΔt)⁴ shape at 1/9 of its level (10 % of the carrier per 100 periods); unstable if a discrete mode has ωΔt > √3 |

The two explicit schemes are a pair. Their stability limits differ (√3 against 2.83), so a mode
between the two limits blows up RK3 and leaves RK4 standing. Each explicit stage still solves a
system here, because the mass operator ∂R/∂u̇ carries the dispersion term B; the stages are not
cheap. To run a scheme in the pilot only, put it in `PILOT_EXTRA` in `cases.py`.

Not in the matrix:
* θ > 0.5: first order, 2.4e-3 per step at θ = 0.6, so the carrier is gone within a few hundred steps;
* the schemes that amplify oscillations: explicit first/second order, the DIRK Radau and Lobatto IIIC variants;
* Gridap's defective `SDIRK_3_3` and `SDIRK_LobattoIIIC_2_2` (rule 15).

## The matrix (`cases.tsv`, 120 cases in 18 packs)

| tier | what | cases |
|---|---|---|
| **P** pilot | m5 n32 T 1.6 kA 0.25, the 5 integrators, 100 steps: s/step per scheme on the node | 5 |
| **A** core | m5 n32 × T {1.6, 2.2, 3.0} × kA {0.158, 0.25, 0.32} × 5 integrators | 45 |
| **B** resolution | m5 T 1.6 × kA {0.25, 0.32} × n {16, 64} × 5 (rule 38b: does refinement help or hurt?) | 20 |
| **C** variable bed | m6 n32 T 1.6 × kA {0.25, 0.32} × 5 | 10 |
| **D** time step | m5 n32 T 1.6 kA 0.32 × 5 integrators × dt {T/20, T/80} | 10 |
| **E** null control | m3 (𝓝 off, no CIP) n32 × T {1.6, 3.0} × kA {0.158, 0.25, 0.32} × 5 | 30 |

Regenerate after any change: `python3 run/integrator_study/cases.py` (it writes both tables and
prints the per-tier relative cost). Names: `<model>_T<100T>_kA<1000kA>_n<ncell>_<integrator>[_dt20|_dt80]`.

## How it runs

* **Packs.** One array task runs one pack of up to 8 sequential cases side by side on a **rome 1/8
  node** (each tier is cut into the fewest packs, as equal in size as possible: 18 packs of 5–8) (16 cores, `--mem=28G`, ≈ 3.5 GiB per case; a 𝓝 box case peaked at 2.8–3.0 GB locally),
  one thread each, Julia loaded from the **sysimage** (`-J`, strict freshness check). Cases of one
  tier and model are cost-sorted and grouped, so a pack's cases end at similar times.
* **Segments.** `run_pack.sh` stops each case with `timeout` at `BALFEM_SEG_HOURS` (119 h). The next
  segment `<case>/seg<NN+1>` resumes from the newest checkpoint (every 250 steps). A pack with a
  case cut this way **resubmits itself** (up to `BALFEM_MAX_SEG` = 6 segments). Markers in
  `<case>/`:
  * `FINISHED`: the driver printed `done:` (completed, or stopped by the divergence guard; both are
    results);
  * `FAILED`: an error. A Newton `NaN` here is a **blow-up result**, not a fault; read the message.
  The stub test of this logic (finish, cut + resume, fail, resubmission limit) passed on 2026-11-11.
* **Output.** `$BALFEM_STUDY_OUT/<case>/seg<NN>/` (default `$BALFEM_PROJ/output/integrator_study`):
  `CONFIGURATION.txt`, `run_manifest.toml`, `diagnostics.csv`, VTK every 5 steps (about 160 MB per
  case, about 30 GB in all: check the quota), `run.log`, `exit_status`.

## Launch

```bash
# on Snellius, after pulling v2-solver
cd $HOME/GridapBALFEM.jl/compile && sbatch compile_snellius.sh      # sysimage: performance mode, Nσ = 3
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh P         # pilot: one pack, ~100 steps per case
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh status    # per-case state
# read the pilot's s/step per integrator (run.log), re-estimate the cost below, then:
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh E A       # cheap null control + the core
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh C B D
```

## Cost (to be replaced by the pilot's numbers)

The only measurements are local: m5 n32 under CN at ≈ 58 s/step on a loaded laptop, before the
depth-field memo (−29 %) and Gridap performance mode (≈ −20 %). m3 runs at ≈ 5 s/step. An
SDIRK_2_2 step costs ≈ 2.5 CN steps. In "m3 CN case" units (8000 steps, ≈ 11 h), `cases.py`
weights the tiers A 1314, B 584, C 321, D 365, E 88, about 2670 units in all.

* **Order of magnitude:** ≈ 3 × 10⁴ case-hours. Each case is billed 2 cores (1/8 node per 8
  cases), so ≈ 6 × 10⁴ SBU before the speed-ups and ≈ 3.5 × 10⁴ after.
* **Tier A alone:** about half of that.
* **The explicit pair and SDIRK_4_3** (3–4 stages each) are about three quarters of the bill.
* **Upper bound:** every case that diverges stops early, and at kA 0.32 most CN cases are expected
  to.

**Decide the scope after the pilot.**

## Analysis

Per case, stitch the segments in time order. Each has its own `diagnostics.csv` and VTK series, and
`load_periodic` reads one segment directory. Then:
* `energy_history(r; modes = :groups)`: the carrier against the disturbance;
* the growing modes' rates against the CN pair;
* the carrier loss reported in `CONFIGURATION.txt`.

Judge by STABILITY.md §6 (σ ≤ 0 over the run, carrier loss within the author's budget).
