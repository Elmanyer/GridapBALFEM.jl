# OUTPUT_NAMING_PROPOSAL.md — a standard name for every simulation output directory

> **STATUS: ACCEPTED and IMPLEMENTED (2026-09-07).** `<domain>` moved to position 2 on review.
> The generator is `output_dir_name` in `src/utilities.jl`; every driver calls it.
>
> Scope: the name of the directory a run writes into (`output_dir`), which is what makes a result
> findable, comparable and auditable months later.

---

## 0. Why the current scheme does not survive the parameter space

Today three drivers each build a name their own way:

```
output/local_1d/flume_bc_nonlinear_full_flat_A0.1_T1.6_P1LFE-2      (1-D driver)
output/local_2d/small2d_<tag>_<model>                               (2-D driver)
output/small_bcplane_nonlinear_full_flat_A0.1_T2.0_P1LFE-2          (distributed)
```

Four problems, each of which has already cost something:

1. **The prefix carries no information** — `flume_`, `small2d_`, `small_` describe which *script* ran,
   not what was *simulated*. Two identical physics cases launched from different drivers sort apart.
2. **The discretisation is invisible.** Nothing in a name says `Q2/Q1` vs `Q2/Q2`. That is exactly
   how the equal-order runs and the Taylor-Hood MMS campaign were compared for months as though they
   were the same solver (`CLAUDE.md` rule 12b). **This alone justifies the rework.**
3. **Field order is inconsistent** and some fields are silently absent (`bar` appears only when the
   bed is not flat; `M2` vs `P1LFE-2` are the same model written two ways).
4. **Ad-hoc runs escape entirely** — `ad_none_A005`, `cfl_fine_dt040`, `qq4_A010`, `lwtest_Lx90`
   are unparseable a week later.

---

## 1. The proposed grammar

```
<model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_<extra>…]
```

**`<domain>` sits at position 2** so a listing groups first by vertical basis, then by
dimensionality — the two things that decide whether two runs are comparable at all — before the
physics switches. Fields **1–6 are mandatory and in this order**. Fields 7–9 are mandatory but
carry defaults that are still written out — an absent field must never mean "the default", because
that is precisely what made the equal-order problem invisible. `<extra>` is an open, ordered tail
for anything a study varies deliberately.

### Field 1 — `<model>` · the vertical basis

`P{p_vert}LFE-{M}` → `P1LFE-2`, `P1LFE-3`, `P2LFE-2`

Already the project's canonical name (`CLAUDE.md` §0). ⚠ Retire the legacy `M2` spelling; it is the
same object written a second way and it defeats sorting.

### Field 2 — `<domain>` · geometry class

`1d` (narrow flume, `ny=1`, solid walls) | `2d` | `2dper` (y-periodic) | `2dopen` (open laterals).

⚠ **1-D is always normal-incidence** (`CLAUDE.md` rule 12), so `dir` + `1d` is a contradiction and a
name generator should refuse it — a naming scheme that can express an impossible case will
eventually label one.

### Field 3 — `<wave>` · what is generated, and how

**Generation mechanism is a prefix on the wave type, not a separate field** — `bcplane` reads better
than `plane_bc` and keeps the mechanism attached to the thing it generates.

| prefix | mechanism |
|---|---|
| *(none)* | interior Gaussian wavemaker (`:inner_res`) |
| `bc` | Dirichlet boundary generation |
| `ic` | released initial condition, no forcing |

| stem | content |
|---|---|
| `plane` | regular long-crested wave |
| `ring` | point-source ring wave |
| `bi` | bichromatic pair |
| `multi` | multichromatic, deterministic components |
| `irr` | irregular sea, JONSWAP, long-crested |
| `dir` | directional (short-crested) sea |
| `hump` | released Gaussian hump (initial condition) |
| `still` | rest state / no generation |
| `mms` | manufactured solution |

Examples: `bcplane`, `plane`, `bcirr`, `dir`, `ichump`, `bcring`.

### Field 4 — `<regime>` · `lin` | `nl`

From `regime=:linear|:nonlinear`. Abbreviated because it appears in every name.

### Field 5 — `<nlp>` · nonlinear pressure tier

`none` | `native` | `full`. **Written even when `lin`**, where it is always `none` — a fixed-width
grammar is worth more than three saved characters, and `lin_none` states that the linear model
carries no `𝓝` rather than leaving the reader to infer it.

### Field 6 — `<bed>` · sea-bed slope

| token | meaning |
|---|---|
| `flat` | `flat_bed=true`, `∇h ≡ 0` |
| `bar` | submerged bar / shoal |
| `slope` | constant or smooth slope |
| `step` | smoothed step |
| `bathy` | measured or otherwise irregular bathymetry |

⚠ **`flat` is the *model* switch, not merely a flat depth field.** `flat_bed=true` over a varying
`h` is a different model (the ∇h terms are dropped), and the name must say which was solved.

### Field 7 — `<discr>` · the discretisation **(new, and the point of the exercise)**

```
Q{p_u}Q{p_eta}[-nx{nx}[x{ny}]]
```

`Q2Q1`, `Q3Q2`, optionally `Q2Q1-nx480`. Taylor-Hood is now enforced (rule 2b), so `p_eta` is
redundant *in principle* — write it anyway. The whole lesson of 2026-09 is that a discretisation
which is not visible in the output is a discretisation nobody checks.

Include `-nx…` whenever the study varies resolution (refinement ladders, convergence work);
omit for one-off production runs where the mesh is fixed by the case.

### Fields 8–9 — `<amplitude>` and `<period>`

| forcing | amplitude | period |
|---|---|---|
| regular | `A0.1` | `T1.6` |
| irregular / directional | `Hs0.2` | `Tp2.5` |
| MMS / rest | `A0` | `T0` |

Metres and seconds, always, formatted with `%g` — **significant digits, not fixed decimals**.
`A0.1`, `A0.05`, `A0.001`. Decimal points are kept (`A0.05`, not `A005`): `A005` is ambiguous between
0.05 and 0.005, and that ambiguity is already present in the existing `ad_none_A0075` directories.
⚠ Do **not** pad to a fixed width — `%.2f` would render `A0.001` as `A0.00` and collapse the
small-amplitude tier-agreement runs onto one name.

### Field 10+ — `<extra>` · deliberate variations only

`key` + value, no separator inside the token, appended in the order below so names stay sortable:

| token | when |
|---|---|
| `dt0.04` | the time step is the variable |
| `sdirk33`, `theta`, `rk4` | integrator is not the default `SDIRK_2_2` |
| `P{n}` | duration in wave periods, when it is the variable |
| `r{N}` | MPI rank count, when comparing decompositions |
| `seed{N}` | stochastic sea realisation |
| `q{n}` | `quad_extra` |
| `ad` | AD Jacobians instead of the hand pair |

**Nothing else goes in the name.** Everything else belongs in the run manifest (§3).

---

## 2. Worked examples

| current | proposed |
|---|---|
| `flume_bc_nonlinear_full_flat_A0.1_T1.6_P1LFE-2` | `P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6` |
| `flume_bc_linear_none_flat_A0.001_T1.6_M2` | `P1LFE-2_1d_bcplane_lin_none_flat_Q2Q1_A0.001_T1.6` |
| `small_bcplane_nonlinear_full_bar_A0.1_T2.0_P1LFE-2` | `P1LFE-2_2d_bcplane_nl_full_bar_Q2Q1_A0.1_T2.0` |
| `small_directional_linear_none_flat_Hs0.2_P1LFE-2` | `P1LFE-2_2dopen_dir_lin_none_flat_Q2Q1_Hs0.2_Tp2.5` |
| `th_q3q2_sdirk33` | `P1LFE-2_1d_bcplane_nl_none_flat_Q3Q2_A0.1_T1.6_sdirk33` |
| `thref_dx0125` | `P1LFE-2_1d_bcplane_nl_none_flat_Q2Q1-nx480_A0.1_T1.6` |
| `cfl_fine_dt040` | `P1LFE-2_1d_bcplane_nl_none_flat_Q2Q1-nx480_A0.1_T1.6_dt0.04` |
| `ad_none_A005` | `P1LFE-2_1d_bcplane_nl_none_flat_Q2Q1_A0.05_T1.6_ad` |

**What the grammar buys.** The refinement ladder becomes three names differing in exactly one token,
and `ls` sorts every P1LFE-2 nonlinear flat-bed case together with the mesh as the last thing that
varies. The equal-order/Taylor-Hood confusion becomes impossible to overlook: `…_Q2Q2_…` next to
`…_Q2Q1_…` in the same listing.

---

## 3. The companion manifest — the part that actually closes the gap

A name cannot hold everything, and it should not try. **Every run directory also writes
`run_manifest.json`** with the complete configuration:

```json
{
  "name": "P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6",
  "created": "2026-09-07T14:22:51",
  "git": {"commit": "d3cf670", "branch": "main", "dirty": false},
  "model":   {"M": 2, "p_vert": 1, "c_bdy": [0.0, 0.728, 1.0], "Nsigma": 3},
  "physics": {"regime": "nonlinear", "nl_pressure": "full", "flat_bed": true, "d": 3.5},
  "discr":   {"p_u": 2, "p_eta": 1, "taylor_hood": true,
              "nx": 240, "ny": 1, "Lx": 60.0, "Ly": 0.25, "quad_extra": 0},
  "time":    {"solver": "sdirk", "tableau": "SDIRK_2_2", "dt": 0.04, "T_final": 80.0},
  "forcing": {"wave_gen": "bc", "kind": "regular", "A": 0.10, "T": 1.6,
              "relax": true, "sponge_wL": 0.0, "sponge_wR": 15.0},
  "exec":    {"mpi": false, "ranks": 1, "host": "workstation"}
}
```

**The name is for humans scanning a directory listing; the manifest is for reconstructing a run.**
Had this existed, the pairing discrepancy would have been a one-line `jq` query across every output
directory rather than a three-day investigation.

⚠ **`git.commit` and `git.dirty` are the load-bearing fields.** A result whose code state is unknown
cannot be reproduced, and `dirty: true` is the honest admission that it may not be reproducible at all.

---

## 4. Rules for the generator

1. **One function, `output_dir_name(...)`, in `src/utilities.jl`** — used by every driver, sequential
   and distributed. Three drivers building names three ways is the present defect.
2. **Fields 1–9 are never omitted.** No conditional tokens; an absent field must not encode a default.
3. **`BALFEM_OUTDIR` still overrides**, for scratch and one-off diagnostics. The standard is the
   default, not a cage.
4. **Refuse impossible combinations** (`dir` + `1d`; `lin` + `nlp≠none`) rather than name them.
5. **Never silently overwrite.** If the directory exists, suffix `_v2`, `_v3`, … ⚠ This is not
   cosmetic: on 2026-09-06 a re-execution wrote six finished runs over their own output and destroyed
   them. A naming scheme that cannot collide is a naming scheme that cannot do that.
6. **Fixed field order**, so `sort` groups by model → wave → regime → tier → bed.

---

## 4b. Implementation status (2026-09-07)

| piece | state |
|---|---|
| `output_dir_name` + tokenisers, `unique_output_dir` | ✅ `src/utilities.jl`, exported, **10/10** tests |
| `examples/local_1d/run_flume_1d.jl` | ✅ wired |
| `examples/local_2d/run_small_2d.jl` | ✅ wired |
| `examples/distributed_small/` (5 drivers) | ✅ wired |
| `examples/distributed/` (7 drivers) | ✅ wired |
| `run_manifest.json` (§3) | ⛔ **not implemented** |

⚠ `BALFEM_OUTDIR` still overrides everywhere, so existing launchers that set it are unaffected —
which is why the overnight campaign kept writing to its own directories through this change.

---

## 5. Open questions for you

1. **Length.** `P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6` is 47 characters. Acceptable, or
   should `<domain>` move into the parent directory (`output/1d/…`, `output/2d/…`) instead?
2. **`<discr>` placement.** Field 6 keeps the physics contiguous in 1–5. The alternative is putting
   it last, keeping today's names recognisable — at the cost of burying the thing that caused the
   instability.
3. **Existing output.** Leave as-is, or provide a rename script? Leaving it means two conventions
   coexist; renaming risks breaking analysis scripts that hardcode paths.
4. **Depth `d`.** Currently never in the name, though `kd` is what determines the physics. Add
   `d3.5`, or leave it to the manifest?
5. **Should the manifest come first?** It is independently valuable and much less invasive than
   renaming — it could land now and the naming follow later.
