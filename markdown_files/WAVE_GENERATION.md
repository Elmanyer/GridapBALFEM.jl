# WAVE_GENERATION.md — wave generation, absorption and boundary treatment

**What this file is.** How waves enter and leave the domain: the two generation mechanisms, the
Dirichlet boundary formulation and its polarizations, the sponge and relaxation zone, and the
WaveSpec coupling contract.

Related: [`ARCHITECTURE.md`](ARCHITECTURE.md) §3 (boundary conditions),
[`CONFIGURATION.md`](CONFIGURATION.md) §8 (sponge tuning).

---

## 1. The `wave_gen` selector — two mechanisms

| `wave_gen` | mechanism |
|---|---|
| `:inner_res` | interior Gaussian source in the continuity equation — **line** ⇒ plane wave, **point** ⇒ ring wave |
| `:bc_gen` | Dirichlet boundary generation |
| `:auto` (default) | infers: `wave_bc === nothing` ⇒ `:inner_res`, else `:bc_gen` |

For `:bc_gen` the boundary source is dispatched on the **type** of `wave_bc` — a parametrised
regular plane wave from `A_wave`/`T_wave`/`wave_dir`, a caller-supplied `WaveInput`, or a WaveSpec
`AiryState`. All three feed the same Dirichlet machinery and differ only in how the `WaveInput`
component table is populated. `resolve_wave_gen` validates the choice; `wave_dir` sets the
propagation angle vs +x.

> **The two mechanisms differ by design and by a factor ≈3.7 in delivered amplitude.**
>
> | mechanism | measured amplitude |
> |---|---|
> | interior Gaussian line source | **3.90 × A** |
> | Dirichlet boundary generation | **1.05 × A** |
>
> The interior source is correct at 3.9×: it radiates in *both* directions (the factor 2 in
> `S(x,t)`), and the two trains superpose at the source, which is where the probe maximum sits. The
> Dirichlet boundary *prescribes* the amplitude, to 5 %.
>
> **An interior-source run's `A_wave` is NOT the wave amplitude in the domain. A BC run's is.**

---

## 2. Dirichlet boundary generation (`src/waveinput.jl`)

### 2.1 What is prescribed

The unknowns are `(η, u_1..u_Nσ)` with `u_j` the horizontal velocity **at σ-node `σ_j`**. A
Dirichlet inflow therefore prescribes on the generation boundary:

* `η(x,t)` — a sum of components `A_c cos ψ_c`, `ψ_c = k_c(x cosθ_c + y sinθ_c) − ω_c t + φ_c`;
* `u_j(x,t)` — the **nodal trace of the incident velocity profile at `σ_j`**, split into x/y by the
  component direction `(cosθ_c, sinθ_c)`.

**A boundary wave must be a consistent solution — surface *and* velocity — so that it radiates
cleanly.** That is why every `:bc_gen` source is a `WaveInput`.

`𝖴y` is prescribed **only for directional seas** (`θ_c ≠ 0` present); a long-crested,
normal-incidence sea keeps the y-solid-wall BC (`v ≡ 0` is exactly consistent). **Directional runs
require `y_wall_bc=:open` plus lateral sponges** — a `v ≠ 0` inflow trace is incompatible with
corner wall tags.

### 2.2 Two vertical polarizations

For each component `(ω_c, θ_c)` the nodal amplitude vector `m_j` comes from one of:

* **`:airy`** — sample linear Airy theory at the still-water σ-levels `z_j = (σ_j − 1)d`:
  `u_j = A_c ω_c cosh(k_c d σ_j)/sinh(k_c d)`, `k_c` from `ω² = gk tanh(kd)`. Deep-water guard at
  `k_c d > 20`: `cosh(kdσ)/sinh(kd) → exp(kd(σ−1))`.
* **`:model` (default)** — the **discrete BALFE-M plane-wave eigenmode**, which makes the boundary
  data an *exact* solution of the linearised discrete system:

  ```
  u⃗_j = (g k/ω) [ (M − k²d²B)⁻¹ Φ ]_j · η · k̂          (polarization)
  ω²  = g d k² Φᵀ (M − k²d²B)⁻¹ Φ                        (model dispersion)
  ```

  Per component: solve `k_c` from the **model** dispersion (scalar Newton seeded with the Airy `k`;
  `B ≤ 0` ⇒ `M_eff > 0` always), then `m = (gk/ω) M_eff⁻¹Φ`.

Within the applicable-`kd` band the two differ by <2 % (the model's own dispersion error), but
`:model` launches a clean single rightward discrete mode with minimal spurious radiation. It
prescribes the exact discrete transport `Φ·Uamp = Aω/(kd)`.

**Consistency identity, used as a unit test:** continuity closure requires
`(dk/ω)(gk/ω) Φᵀ M_eff⁻¹ Φ = 1` — exactly the model dispersion relation.

### 2.3 Start-up

* **Hann ramp (default)** — multiply the whole BC by
  `r(t) = t < T_ramp ? 0.5(1 − cos(πt/T_ramp)) : 1`, default `T_ramp = 2·T_p`. This avoids the
  impulsive start (zero interior against a finite boundary state) that would excite spurious
  transients. A cold start `u0 = 0` is then compatible, since the BC is zero at `t=0`.
* **Hot start (`ic_from_bc=true`)** — interpolate the incident field as the initial condition
  (the η and `𝖴x`/`𝖴y` closures at `t=0` are valid in the whole domain, not just at the boundary).
  Requires `T_ramp=0`; gives zero start-up transient in the linear regime.

> **BC closures must be ForwardDiff-safe in `t`** — Gridap computes the Dirichlet time derivative by
> AD over `t`. Keep the bodies type-generic; **no `t::Float64` annotations inside**.

### 2.4 Reflection and the relaxation zone

A clamped Dirichlet boundary is perfectly reflective for the *deviation* field — anything the domain
radiates back at it. Mitigation, in order:

1. **Far-side sponge** — mandatory in practice, width ≥ `λ_p`.
2. **Optional generation/absorption relaxation zone** (`relax_bc`, default off): in a zone of width
   `relax_width` adjacent to the inflow, add

   ```
   continuity: + ∫ q μ_g(x) (η − η_inc(x,t))
   momentum:   + ∫ μ_g(x) ( Wx·M(𝖴x − 𝖴x_inc) + Wy·M(𝖴y − 𝖴y_inc) )
   ```

   with the same quadratic profile as the sponge (`μ_g` maximal at the boundary, → 0 at the inner
   edge). This is classical relaxation-zone generation: it absorbs outgoing deviations while
   enforcing the incident state. Linear in `u`, so the Jacobian contribution is exact and trivial.

   Measured: in-zone amplitude error **≤ 0.1 %**, phase **0.9°**, reflection 0.67 %, and absorption
   **145×** its control.

Driver kwargs on both drivers: `wave_bc`, `bc_side`, `bc_profile`, `T_ramp`, `ic_from_bc`,
`relax_bc`, `relax_width`. Transient trial spaces work sequentially and distributed.

---

## 3. The sponge — and why it must damp η

```
sponge: + ∫ μ ( q·η + (Wx⋅(𝗠⋅Ux)) + (Wy⋅(𝗠⋅Uy)) )
```

quadratic profile, up to 4 boundaries in 2-D, corner regions clamped at `mu_max` (not `2·mu_max`).

> **🔴 The open-boundary spurious mode is η-DOMINATED, so a velocity-only sponge cannot absorb it —
> at any `mu_max`.** At a free (`x_wall_bc=false`) outflow the model supports a boundary-localised
> mode carrying large surface displacement with little velocity. A velocity-only sponge
> (`∫ μ (W⋅𝗠U)`) is structurally blind to it. **The `+∫ μ q η` continuity term — giving
> `ηt = … − μη` — is what absorbs it, and it is not optional.** Closed-basin results are
> bit-identical (μ ≡ 0 there).

**How to recognise the mode in a run:** bin `max|η|` by `x` from the VTK output. The mode

* peaks **at the last boundary node** and decays exponentially inward;
* is spatially **disconnected** from the incident field — 60× its immediate upstream neighbour while
  the wave front is still mid-domain;
* grows with an e-folding time comparable to `Tp`;
* has `|u|/|η| ≈ 0.4`, against ≈ `ω/kd` (≈0.9) for a genuine wave and ≈1.9 in the incident train.

Once it exceeds the sea state, `eta_max` in the log stops reporting the physics and reports the
mode. `test/local/test_boundary_modes_1d.jl` reproduces it deliberately as a negative control:
`dmp/int` climbs 0.02 → 2.7 → 37.8 → **64.8** while the interior wave stays flat at 3.5e-3, in
**6 s of simulated time**.

---

## 4. WaveSpec coupling contract

* `WaveSpec.AiryState` supplies per-bin `(ω_i, θ_j)`, amplitudes `A_ij = get_amplitudes(state)`
  (energy-preserving normalisation built in) and **seeded** random phases
  `φ_ij = get_random_phases(state)` — deterministic given `state.seed`, hence reproducible across
  sessions **and MPI ranks**.
* The solver **snapshots these into plain arrays at construction** (a `WaveInput` struct). No
  WaveSpec types cross into the residual or the FE spaces; the dependency is constructor-level only.
* **`g` and `k` consistency:** WaveSpec uses `g = 9.80665` and solves `k` with `state.h`. The
  converter **discards WaveSpec's `k`** and re-solves each `k_c` with the *solver's* `g` and the
  boundary depth `d` (`:airy`) or the model dispersion (`:model`). It warns if `state.h ≠ d` — a
  mismatch means the spectrum was sampled at another depth.
* Requires constant depth `d` along the generation boundary (checked by the driver).

**Validation on record:** `spectral_fidelity.jl` at 60 `Tp` — Hs transfer **1.023**, dispersion
14/14 within 5 % of the model wavenumber, incident amplitudes 9/14 within 10 %;
`bc_irregular_sea.jl` — near-inflow Welch Hs ratio **0.979**. Gates: `test_waveinput` 30/30,
`test_bc_generation` 11/11, `test_bc_spectrum` 8/8 (Goda–Suzuki), distributed 4/4.

---

## 5. Open follow-ups

* `:bottom` / `:top` generation sides.
* Second-order (bound-wave) corrections for irregular BC generation.
* Sheared-current focusing needs an ambient-current term — a modelling extension, not a gap.
