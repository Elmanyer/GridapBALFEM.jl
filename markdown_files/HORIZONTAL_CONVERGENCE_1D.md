# HORIZONTAL_CONVERGENCE_1D.md — the 1-D horizontal convergence matrix

> **What this file is.** The complete result of the 1-D spatial convergence campaign (Campaign C,
> 2026-09-11/12): **18 studies, 84 solver runs, 18 OK, 0 errors, 44.1 core-hours**, sweeping the
> three Taylor-Hood horizontal pairings against all six physics models.
>
> **Why it has its own file.** `PLANNED_CAMPAIGNS.md` §2 records *that* the campaign ran and its two
> headline conclusions; `OPEN_ISSUES.md` §0b records the *one defect* it exposed. Neither is the
> place for the full matrix, and the matrix is the thing that has to be read as a whole — the
> interesting feature is a PATTERN ACROSS PAIRINGS, not any single row.
>
> **Data.** `output/local/mms/phaseB/shard_NN_v2.csv` (schema v2 — the per-level sequences are
> persisted, not just the fitted slope). Per-level diagnostics, including the achieved Newton
> residual at every level, under
> `output/local_1d/mms_convergence_campaigns/<model>/<pair>/nx<N>/diagnostics.csv`.

---

## 1. Configuration

Identical across all 18 studies; only `p_u`, `regime`, `flat_bed` and `nl_pressure` vary.

| | |
|---|---|
| vertical basis | P1LFE-2 (`M=2`, `p_vert=1`, published Yang & Liu nodes) |
| domain / mode | 1-D (`:d1`, `ny=1`, `y_wall_bc=:wall` — rule 12), static MMS (`ω=0`) |
| ladder | `nx = 4, 8, 16, 32, 64` (Q2/Q1, Q3/Q2) · `nx = 4, 8, 16, 32` (Q4/Q3 — 64 sits on the ~1e-10 algebraic floor, `PLANNED_CAMPAIGNS.md` §0) |
| manufactured amplitude | `a_eta = 0.8` |
| tolerances | `nl_tol = 1e-14`; achieved Newton residual **~1e-15 at every level of every study** |
| optimum | `p_η → p_u`, `p_u → p_u+1` (the Taylor-Hood L² pattern encoded in `run_conv_study`) |

**Model key.** 1 linear/flat/`none` · 2 linear/var/`none` · 3 nonlinear/flat/`none` ·
4 nonlinear/var/`none` · 5 nonlinear/flat/`native` · 6 nonlinear/var/`native`.

---

## 2. The matrix

### Q2/Q1

| model | regime | bed | `nl_pressure` | pair | `p_η` pairwise (nx 4→64) | `p_η` | opt | `p_u` pairwise | `p_u` | opt | `e_η`(fine) | `e_u`(fine) | min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | linear | flat | `none` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.351;2.082;2.015;2.003` | **2.003** | 3 | 1.700e-04 | 3.329e-05 | 19 |
| m2 | linear | var | `none` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.324;2.074;2.015;2.003` | **2.003** | 3 | 1.700e-04 | 3.195e-05 | 36 |
| m3 | nonlinear | flat | `none` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.463;2.117;2.023;2.004` | **2.004** | 3 | 1.701e-04 | 2.075e-05 | 111 |
| m4 | nonlinear | var | `none` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.605;2.209;2.051;2.012` | **2.012** | 3 | 1.701e-04 | 1.598e-05 | 89 |
| m5 | nonlinear | flat | `native` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.454;2.112;2.022;2.004` | **2.004** | 3 | 1.701e-04 | 2.094e-05 | 164 |
| m6 | nonlinear | var | `native` | Q2/Q1 | `1.986;1.996;1.999;2.000` | **2.000** | 2 | `2.585;2.201;2.050;2.012` | **2.012** | 3 | 1.701e-04 | 1.632e-05 | 155 |

### Q3/Q2

| model | regime | bed | `nl_pressure` | pair | `p_η` pairwise (nx 4→64) | `p_η` | opt | `p_u` pairwise | `p_u` | opt | `e_η`(fine) | `e_u`(fine) | min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | linear | flat | `none` | Q3/Q2 | `2.990;2.997;2.999;3.000` | **3.000** | 3 | `3.497;3.859;3.965;3.991` | **3.991** | 4 | 5.261e-07 | 3.225e-08 | 38 |
| m2 | linear | var | `none` | Q3/Q2 | `2.990;2.997;2.999;3.000` | **3.000** | 3 | `3.604;3.899;3.976;3.994` | **3.994** | 4 | 5.261e-07 | 2.433e-08 | 30 |
| m3 | nonlinear | flat | `none` | Q3/Q2 | `2.971;2.922;2.756;2.450` | **2.450** | 3 | `2.772;3.185;3.579;3.687` | **3.687** | 4 | 9.769e-07 | 1.317e-07 | 197 |
| m4 | nonlinear | var | `none` | Q3/Q2 | `2.971;2.922;2.756;2.449` | **2.449** | 3 | `3.074;3.558;3.792;3.708` | **3.708** | 4 | 9.774e-07 | 6.043e-08 | 133 |
| m5 | nonlinear | flat | `native` | Q3/Q2 | `2.971;2.922;2.757;2.450` | **2.450** | 3 | `2.784;3.187;3.582;3.688` | **3.688** | 4 | 9.767e-07 | 1.306e-07 | 347 |
| m6 | nonlinear | var | `native` | Q3/Q2 | `2.971;2.922;2.756;2.450` | **2.450** | 3 | `3.084;3.566;3.795;3.700` | **3.700** | 4 | 9.772e-07 | 6.157e-08 | 342 |

### Q4/Q3  (ladder stops at `nx = 32`)

| model | regime | bed | `nl_pressure` | pair | `p_η` pairwise (nx 4→64) | `p_η` | opt | `p_u` pairwise | `p_u` | opt | `e_η`(fine) | `e_u`(fine) | min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | linear | flat | `none` | Q4/Q3 | `3.994;3.999;4.000` | **4.000** | 4 | `4.332;4.661;4.596` | **4.596** | 5 | 2.292e-08 | 3.292e-09 | 36 |
| m2 | linear | var | `none` | Q4/Q3 | `3.994;3.999;4.000` | **4.000** | 4 | `4.396;4.655;4.527` | **4.527** | 5 | 2.294e-08 | 2.928e-09 | 43 |
| m3 | nonlinear | flat | `none` | Q4/Q3 | `3.995;4.001;4.002` | **4.002** | 4 | `3.434;3.728;4.298` | **4.298** | 5 | 2.285e-08 | 2.054e-08 | 90 |
| m4 | nonlinear | var | `none` | Q4/Q3 | `3.994;4.000;4.000` | **4.000** | 4 | `3.704;4.202;4.681` | **4.681** | 5 | 2.291e-08 | 7.835e-09 | 143 |
| m5 | nonlinear | flat | `native` | Q4/Q3 | `3.995;4.001;4.002` | **4.002** | 4 | `3.436;3.729;4.299` | **4.299** | 5 | 2.285e-08 | 2.044e-08 | 342 |
| m6 | nonlinear | var | `native` | Q4/Q3 | `3.994;4.000;4.000` | **4.000** | 4 | `3.706;4.204;4.682` | **4.682** | 5 | 2.291e-08 | 7.798e-09 | 335 |

⚠ **Read the pairwise SEQUENCE, not the final number or the fit** (rules 32/33). A rising sequence
with the error still dropping is *pre-asymptotic*; a flat sequence with a stalled error is
*saturated*; a **falling** sequence is a lower-order component taking over. The three appear in this
table and they mean different things.

---

## 3. The three features, pairing by pairing

### 3.1 Q2/Q1 — `p_u` is a FULL ORDER below optimal, on every model

`p_η = 2.000` exactly (optimal 2) on all six. `p_u = 2.003–2.012` against an optimal **3**.

This is not a transient and not a floor:

* the sequence **descends onto** 2.00 — `2.351 → 2.082 → 2.015 → 2.003` — so it is a converged rate,
  not a pre-asymptotic approach to 3;
* `e_u(fine) ≈ 1.6–3.3e-05` sits **five orders above** the ~1e-10 algebraic floor, so it is not
  saturation;
* it is identical for linear and nonlinear, flat and sloping bed, `:none` and `:native`. Nothing in
  the physics moves it.

It reproduces the tier-2 result, but tier 2 measured only models 1 and 3 — both flat-bed, both
`:none` — so by rule 4 it could not exercise `∇h` and never touched the `𝓝` blocks. **The finding is
now established on variable bathymetry and on the production `:native` tier.**

### 3.2 Q3/Q2 — `η` degrades for the nonlinear models; `u` is essentially optimal for all

| | `p_η` | `p_u` |
|---|---|---|
| linear (m1, m2) | `2.990;2.997;2.999;3.000` → **3.000** (opt 3) ✅ | **3.99** (opt 4) ✅ |
| nonlinear (m3–m6) | `2.971;2.922;2.756;2.450` → **2.450** (opt 3) ⚠ | **3.69–3.71** (opt 4) |

The `η` sequence **falls monotonically with refinement**, which is the one unambiguous anomaly in the
whole matrix (§4). `u` is short of its optimum by ~0.3, not by an order, and its sequence is *rising*
(`2.772 → 3.185 → 3.579 → 3.687`), i.e. still approaching 4 from below.

⚠ **Q3/Q2 linear is the only configuration in the matrix optimal in BOTH fields simultaneously**
(3.000 / 3.99). It is the natural reference point for everything else here.

### 3.3 Q4/Q3 — `η` recovers optimal; `u` is short but rising

`p_η = 4.000–4.002` (optimal 4) on all six, **including the nonlinear models** — the Q3/Q2 `η`
degradation does *not* appear. `p_u = 4.30–4.68` against an optimal 5, with sequences climbing
(`3.434 → 3.728 → 4.298`; `3.704 → 4.202 → 4.681`) and `e_u` still falling — pre-asymptotic rather
than deficient, and the ladder stops at `nx=32` precisely because 64 would sit on the algebraic
floor.

---

## 4. ⚠ THE FEATURE THAT MATTERS MOST: NO PAIRING IS DEFICIENT IN BOTH FIELDS

Collapsing the matrix to its structure:

| pairing | `p_η` | `p_u` |
|---|---|---|
| **Q2/Q1** | ✅ optimal (2.000) | ⛔ **one full order short** (2.00 vs 3) |
| **Q3/Q2** | ✅ optimal linear · ⚠ degrades nonlinear (2.45 vs 3) | ✅ ~optimal (3.69–3.99 vs 4) |
| **Q4/Q3** | ✅ optimal (4.000) | ⚠ short but rising (4.30–4.68 vs 5) |

**Every pairing is optimal in at least one field, and no pairing is deficient in both.** The
deficiency MOVES between `η` and `u` as the pairing changes, and it changes character as it moves —
a full order at Q2/Q1, a refinement-dependent decay at Q3/Q2, a pre-asymptotic shortfall at Q4/Q3.

**This is why the matrix does not localise a failure.** A genuine defect in a specific operator would
be expected to leave a *consistent* signature — the same field, the same size, the same sign, across
pairings — because the operator does not know which horizontal spaces it is being assembled on. What
is observed instead is a deficiency that is present in exactly one field at a time and that swaps
fields when `p_u` increments by one. **No single wrong term in the residual obviously produces
that.**

⚠ **So the honest reading is that WE CANNOT SAY FROM THIS MATRIX ALONE WHETHER THE SOLVER IS FAILING
AT ALL.** These may simply be the convergence properties of this mixed formulation: `η` enters
momentum undifferentiated through `∇·v` and plays the pressure role of a Stokes system (rule 2b), and
mixed methods routinely deliver different orders in the two fields, with the lowest-order pair of a
family behaving worse than the theory for higher members. The "optimum" column here is the
Taylor-Hood L² pattern (`u → p_u+1`, `η → p_u`); **whether those estimates transfer to a dispersive,
non-hydrostatic depth-integrated system has not been established anywhere in this project.** Quoting
a measured rate as "sub-optimal" presumes an optimum that is, for this system, assumed rather than
derived.

**Do not report any single row of this table as evidence that the solver is broken.**

---

## 5. The one thing in the matrix that cannot be normal

There is exactly one discriminator available, and it is the **direction of the trend under
refinement**:

* a **converged but sub-optimal** rate (Q2/Q1 `p_u` → 2.003, flat) is consistent with "this is what
  this pairing does" — it is stable, reproducible, and physics-independent;
* a **rising** rate (Q4/Q3 `p_u`, Q3/Q2 nonlinear `p_u`) is pre-asymptotic — it is not yet evidence
  of anything;
* a rate that **decays monotonically as the mesh is refined** cannot be a fixed property of a
  discretisation that is converging. It means a lower-order error component is overtaking the
  leading one, and it **gets worse on better meshes**.

Only the Q3/Q2 nonlinear `p_η` does the third thing: `2.971 → 2.922 → 2.756 → 2.450`, against linear
models holding 3.000 on the identical basis, ladder and mesh. That — and only that — is treated as a
defect, and it is the subject of **`OPEN_ISSUES.md` §0b**, where it is shown *not* to be algebraic,
not the `𝓝` blocks, not `∇h` and not quadrature (each eliminated by measurement), leaving the
nonlinear advection block.

Everything else in this matrix is recorded as **observed behaviour of the discretisation, cause
unattributed**.

---

## 6. Secondary observations worth keeping

* **`e_η(fine)` is set by the pairing, not the physics** — 1.700e-04 (Q2/Q1), 5.26e-07 (Q3/Q2
  linear), 2.29e-08 (Q4/Q3), essentially constant across all six models within each pairing. `e_u`
  by contrast varies ~2× with the physics. The surface error is a property of the discretisation;
  the velocity error is not.
* **`:native` costs nothing in rate.** m3 vs m5 and m4 vs m6 differ by ~1e-3 in both `p_η` and `p_u`
  at every pairing — independent confirmation on a convergence ladder that the `(1, 2, 4, 5)` hierarchy
  is dynamically negligible.
* **Bathymetry costs nothing in rate either**, but it does improve `e_u` by ~2× at Q2/Q1 and Q3/Q2
  (m4/m6 against m3/m5) — a level effect, not a rate effect.
* **Cost is driven by the regime, not the pairing.** A linear Q3/Q2 study runs in 30 min; the
  nonlinear one on the identical ladder takes 197–347 min. See `CAMPAIGN_COST.md` §2.

---

## 7. What would actually settle §4

Neither of these has been run; both are cheap relative to what has been spent.

1. **Extend Q4/Q3 to `nx = 64`.** If the nonlinear `p_η` degradation is a fine-mesh effect, it should
   appear there too once the fourth-order part has decayed. If `p_η` stays at 4.000 while Q3/Q2 falls
   to 2.45, the effect depends on the pairing and not only on the mesh — which would be a strong
   constraint on any explanation. ⚠ Check `e_u` against the ~1e-10 floor first; it may be
   unmeasurable at that level.
2. **Sweep the manufactured amplitude** (`a_eta = 0.8 → 0.4 → 0.2`) at Q3/Q2 nonlinear. A consistency
   error in a quadratic term scales with amplitude, so the crossover mesh should move; if the rate
   degrades identically at every amplitude, the offending term is not quadratic in the solution.

A third, and the one to do first because it costs nothing: **derive the expected orders for this
system** rather than inheriting the Taylor-Hood pattern by analogy. Until that exists, "optimal" in
the table above is a convention, and §4 cannot be resolved by more runs.
