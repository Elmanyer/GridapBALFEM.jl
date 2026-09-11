# MODEL.md — the BALFE-M model and its discrete residual

**What this file is.** The mathematical reference for the solver: how the 3-D equations become
`Nσ+2` horizontal PDEs, what each vertical σ-tensor is, and how every term of the global residual is
written in Gridap operators. It is the specification `src/problem.jl` is audited against.

**Authority.** The LaTeX project `BALFEM_models/` is the source of mathematical truth; this file is
its solver-facing condensation plus the code mapping. If you change the maths, change both.
The term-by-term classification used as the residual's *specification* lives in
`BALFEM_models/NumericalImplementation/GridapImplementation.tex`, §`subsec: term classification`.

Related: [`ARCHITECTURE.md`](ARCHITECTURE.md) (how the code is laid out),
[`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) (what is proven about it).

---

## 1. The model in one paragraph

**BALFE-M** — *Basis-Agnostic Layer-integrated Finite Element*, `M` vertical elements — generalises
Yang & Liu (2024, *JFM* 999 A32) LFE-M to an **arbitrary vertical FE basis**; piecewise-linear
elements (`p=1`) recover their model as the lowest member. It is a depth-integrated,
**non-hydrostatic** free-surface wave model. The water column `σ = (z+h)/H ∈ [0,1]` is discretised
with `M` vertical elements of order `p`, so

```
u_h(x, σ, t) = Σ_j u_j(x, t) φ_j(σ),      Nσ = num_free_dofs(U_phi) = M·p + 1
```

The vertical velocity `w` and the non-hydrostatic pressure `p_nh` are eliminated **analytically** —
no pressure Poisson solve — into small, precomputable vertical σ-tensors. What remains is `Nσ+2`
coupled 2-D PDEs in `(H, u_1…u_Nσ)`, `H = d + η`, solved on a horizontal FE mesh by Gridap.

The structure is **(small dense vertical algebra) ⊗ (large sparse horizontal FE)**: cheap dense
vertical work built once, wrapped by the large sparse horizontal solve.

**Naming, kept strictly distinct** (defined in the LaTeX at `\label{par: nomenclature}`):

| name | meaning |
|---|---|
| **BALFE-`M`** | the family this project derives, arbitrary vertical basis |
| **P`p`LFE-`M`** | a concrete member implemented/run/tabulated (`Pp` = the basis). Every instantiation in this repo is `p_vert=1` |
| **LFE-`M`** | Yang & Liu (2024)'s published piecewise-linear models — used only when citing or comparing |

Tables comparing our numbers against theirs must not label both sides the same way.

---

## 2. Vertical σ-tensors (built once, `src/vertical.jl`)

On the σ-mesh (`M` elements, order `p`, optimised nodes `c_bdy`) build `φ_j` (basis), `φ_j'`,
`varphi_j = ∫_0^σ φ_j` (solved as a BVP `dφ/dσ = φ_j`, `φ(0)=0`, degree `p+1`), and
`Φ_j = ∫_0^1 φ_j`. Then

```
M2[i,j]   = ∫ φ_i φ_j                        → M^V    vertical mass
M3[i,j,k] = ∫ φ_i φ_j φ_k                    → 𝓜^V    horizontal advection
G3[i,j,k] = ∫ (σΦ_k − varphi_k) φ_j' φ_i     → 𝓖^V    vertical advection
A2[i,j]   = ∫ φ_i θ_j                        → A^V    linear pressure, ∇h coupling
K2[i,j]   = ∫ θ_j (varphi_i − σφ_i)          → K^V    linear pressure, ∇H coupling
A3[i,j,k] = ∫ Θ_kj φ_i                       → 𝓐^V    nonlinear pressure, ∇h
K3[i,j,k] = ∫ Θ_kj (varphi_i − σφ_i)         → 𝓚^V    nonlinear pressure, ∇H
P[i,j,c]  = ∫ θ_j[c] · φ_i_int               → P^V    LEADING pressure (P[:,:,3] = −B)
Pcal[i,k,j,c] = ∫ Θ_kj[c] · φ_i_int          → 𝓟^V    nonlinear leading pressure
B[i,j]    = −∫ varphi_i varphi_j                      diagnostic (dispersion carrier)
```

with the shape vectors

```
θ_j  = [ φ_j , σφ_j , varphi_j ]                                          (3 components)
Θ_kj = [ σΦ_k φ_j , Φ_k varphi_j , φ_j φ_k , σφ_j φ_k , varphi_j φ_k ,
         σΦ_j φ_k' − varphi_j φ_k' ,
         σΦ_j φ_k + σ²Φ_j φ_k' − varphi_j φ_k − σ varphi_j φ_k' ,
         σΦ_j φ_k − varphi_j φ_k ]                                        (8 components)
```

`assemble_vertical_tensors` builds all of these once; `build_problem` reshapes them into constant
`TensorValue` / `ThirdOrderTensorValue` objects (peel the component index, store as `[i,k,j]`).

**Index order rule.** Store `𝓜^V, 𝓖^V, 𝓐^V, 𝓚^V` as `[i,k,j]` =
`[output/test layer, u_k or ∇·(Hu_k) layer, ∇u_j or u_j layer]`, so contraction over the trailing
two indices gives the mode-`i` momentum contribution directly.

**Optimised vertical node positions** (Yang & Liu Table 1; grid optimisation `Δσ_top ≈ 2.94/kd_max`):

| model | `c_bdy` | applicable `kd` |
|---|---|---|
| P1LFE-2 | `[0, 0.728, 1]` | ≈ 10.9 |
| P1LFE-3 | `[0, 0.726, 0.925, 1]` | ≈ 39.2 |
| P1LFE-4 | `[0, 0.745, 0.923, 0.977, 1]` | ≈ 127.9 |

**Notation bridge.** LaTeX `M^V, 𝓜^V, 𝓖^V, A^V, K^V, 𝓐^V, 𝓚^V, Φ, φ_j` ↔ code
`Mmat, Mcal, Gcal, A, K, Acal, Kcal, Phi/D/C, w_j = −φ_j_int`.

---

## 3. The single scalar residual

Gridap's `MultiFieldFESpace` solver consumes **one scalar** = total virtual work `∫_Ω R·v`. With
`U = [u_1…u_Nσ]` stacked, `V = [v_1…v_Nσ]`, `q` the continuity test:

```
Global = ∫_Ω [ q ∂H/∂t − ∇q·(H U)·Φ                                       (mass)
    + ( H·M^V·U̇                                                          (acceleration)
      + H·F_M(U) + F_G(H,U)                                              (advection)
      + gH∇η·Φ                                                           (gravity, η = H−h)
      − H[ ∇h·(L:A^V + N⫶𝓐^V) + ∇H·(L:K^V + N⫶𝓚^V) ] ) · V              (slope pressure)
    − H²[ (L:P^V) + (N⫶𝓟^V) ] · (∇·V) ] dΩ                              (LEADING pressure)
```

with `F_M(U)_i = Σ_kj 𝓜_ikj (u_k·∇u_j)`, `F_G_i = Σ_kj 𝓖_ikj (∇·[Hu_k]) u_j`,
`L_j = [−u̇_j·∇h, u̇_j·∇H, −∇·(Hu̇_j)]` (3 components), `N_kj` the 8 components of §5.

**Sign convention.** All three pressure blocks are positive integrals **subtracted** uniformly:

```
Global = R_mass + R_Acc + R_Adv + R_Grav − R_P − R_lin − R_nonlin ,
R_P = +∫ H²[(L:P^V) + (N⫶𝓟^V)]·(∇·V)
```

`B_stored = −B̃ ≤ 0`, and the explicit `(−1)` factors in the `R_P` and slope-pressure terms are
**load-bearing** — do not "simplify" them away.

> **`R_P` is the entire frequency dispersion of the model.** It is the only `O(η)` non-hydrostatic
> term on a flat bed. Only the *boundary* part of the integration by parts that produces it
> vanishes; the volume part remains and must be assembled. Without `R_P` the model degenerates to
> non-dispersive shallow water. This is also why **`p_u ≥ 2` is mandatory**: `Q1` elements
> cannot represent the second-order pressure coupling and zero the term.

Wave generation adds `−∫ q S(x,t)` to mass; the sponge adds `+∫ μ q η` to mass and
`+∫ μ (M^V U)·V` to momentum.

---

## 4. The residual in Gridap operators (`src/problem.jl`)

With `∂x(f) = e_x⋅∇(f)`, `∂y(f) = e_y⋅∇(f)`:

```julia
DU  = ∂x(Ux)+∂y(Uy);   DW = ∂x(Wx)+∂y(Wy)          # per-layer div of TRIAL / TEST, VectorValue{Nσ}
UgH = ∂x(H)*Ux+∂y(H)*Uy;   Ugh = ∂x(d)*Ux+∂y(d)*Uy;   S = H*DU + UgH
ū   = Operation(VectorValue)(𝚽⋅Ux, 𝚽⋅Uy)

mass  : ∫ q*ηt − H*(∇(q)⋅ū) − q*src(x,t)             # source enters with a MINUS
acc   : ∫ H*( (Wx⋅(𝗠⋅Uxt)) + (Wy⋅(𝗠⋅Uyt)) )
grav  : − ∫ (g/2)*(H*H − d*d)*(𝚽⋅DW)  − ∫ g*η*(𝚽⋅(Wx*∂x(d)+Wy*∂y(d)))
        # linear branch: same identity with H→h. The −η∇h half is IDENTICAL in both regimes
        # and is therefore assembled ONCE, outside the branch.
disp  : − ∫ H²*(sP⋅DW),   sP = P[1]⋅L1 + P[2]⋅L2 + P[3]⋅L3
        L1 = −(∂x(d)*Uxt + ∂y(d)*Uyt)
        L2 =   ∂x(H)*Uxt + ∂y(H)*Uyt
        L3 = −( H*(∂x(Uxt)+∂y(Uyt)) + L2 )
        # P_full=false → sP = P[3]⋅L3 only;  linear → −∫ d²*((P[3]⋅L3lin)⋅DW)
adv   : ∫ H*( double_contraction(𝗠3,TMx)⋅Wx + double_contraction(𝗠3,TMy)⋅Wy )
          + ( double_contraction(𝗚3,TGx)⋅Wx + double_contraction(𝗚3,TGy)⋅Wy )
        TMx = Ux⊗∂x(Ux)+Uy⊗∂y(Ux);  TMy = Ux⊗∂x(Uy)+Uy⊗∂y(Uy);  TGx = S⊗Ux;  TGy = S⊗Uy
linP  : − ∫ H*( ∂x(d)*πAx + ∂y(d)*πAy + ∂x(H)*πKx + ∂y(H)*πKy )
        πAx = (Wx⋅A[1])⋅L1 + (Wx⋅A[2])⋅L2 + (Wx⋅A[3])⋅L3     (πAy, πKx, πKy analogous)
nlP   : components {3,6,7,8} native; {1,2,4,5} split (§5); + the 𝓟-part of R_P (uses Pcal)
sponge: + ∫ μ*( q*η + (Wx⋅(𝗠⋅Ux)) + (Wy⋅(𝗠⋅Uy)) )
```

**Gravity is assembled in integrated-by-parts energy form** with the still-water baseline
subtracted, `−(g/2)(H²−d²)(𝚽·DW)`. That makes the discrete rest state exactly force-free, so an
undisturbed surface at an open wall stays at rest instead of drifting. The exact identity behind it
is

```
∇((H²−h²)/2) = H∇H − h∇h = H∇η + η∇h    ⇒    H∇η = ∇((H²−h²)/2) − η∇h
```

— **two** pieces. Keeping only the first is a real defect (it was one; see `VERIFIED_SCOPE.md` §5).

**The linear-pressure `L` is built as three separate `VectorValue{Nσ}` fields** `L1,L2,L3`, with
`L:A^V = A[1]·L1 + A[2]·L2 + A[3]·L3`. Keeping the components separate (rather than one
`VectorValue{3}`-of-vectors) is what keeps the term well-typed and first-order.

---

## 5. The nonlinear pressure `𝓝` (`src/nlpressure.jl`)

Every nonlinear-pressure component `c` appears in **three** residual blocks with different
prefactors:

| block | prefactor | integral shape |
|---|---|---|
| bed-slope `𝓐` | `H ∂_αh` | `−∫ H ∂_αh (W_α ⋅ (𝓐3[c] ⊡ N^(c)))` |
| surface-slope `𝓚` | `H ∂_αH` | `−∫ H ∂_αH (W_α ⋅ (𝓚3[c] ⊡ N^(c)))` |
| leading `𝓟` (part of `R_P`) | `H²`, tested by `∇·V` | `−∫ H² (𝗣3[c] ⊡ N^(c)) ⋅ D_W` |

`C⁰` spaces admit only **first** derivatives of the unknowns, which splits the components into
classes:

* **Native, `c ∈ {3,6,7,8}`** (`nl_pressure=:native`) — first-order everywhere. `{6,7,8}` are outer
  products of the building blocks `𝖺, 𝖻, 𝖲` (the `1/H` cancels one prefactor `H`); `c=3`'s only
  second derivative is the **analytic** bed Hessian. Assembled directly in all three blocks,
  sequential and distributed.
* **`c ∈ {1,2,4,5}`** (`nl_pressure=:full`) carry second derivatives of the unknowns and are handled
  by half:
  * **∇h half → exact integration by parts onto the test** (§5.1). Uses the analytic bed Hessian;
    machine-verified at 4.0e-15.
  * **∇H half + 𝓟 part → frozen L² projections.** The `∂²η` factor is irreducible, so per step
    project `𝖲` and `𝖻` onto the velocity FE space (SPD mass solves — direct factorisation
    sequentially, CG + Jacobi distributed) and use `∂_a(π𝖲), ∂_a(π𝖻)` **lagged one step**.

All `𝓝` blocks are `O(A²–A³)` and treated **quasi-Newton** (they enter the residual, not the
Jacobian).

> ⚠ **The frozen projections mean `:full` implements a different operator from the exact one.**
> That is a deliberate, now-quantified approximation — see `VERIFIED_SCOPE.md` §4. Do not treat its
> flat MMS rate as a bug.

### 5.1 Exact-IBP algebra for the ∇h half

With `Ψ_c := H ∂_αh (W_α ⋅ 𝓐3[c]) ∈ TensorValue{Nσ,Nσ}` (test contracted into the FIRST index),
`(Ψ·u)_k = Σ_j Ψ_{kj}u_j`, `(u·Ψ)_j = Σ_k u_kΨ_{kj}`, and `M := Σ_a (∂_a𝖲)⊗U_a`
(so `N¹ = −M`, `N² = M + 𝖲⊗D`):

```
−∫Ψ₁⊙N¹ = +∫ Σ_a (Ψ₁·U_a)⋅∂_a𝖲   →IBP→  −∫ 𝖲⋅[∂_x(Ψ₁·U_x) + ∂_y(Ψ₁·U_y)] + ∮
−∫Ψ₂⊙N² = −∫ Σ_a (Ψ₂·U_a)⋅∂_a𝖲 − ∫Ψ₂⊙(𝖲⊗D)      (first piece IBP, sign flipped; second native)
−∫Ψ₄⊙N⁴ = −∫ Σ_a (U_a·Ψ₄)⋅∂_a𝖻   →IBP→  +∫ 𝖻⋅[∂_x(U_x·Ψ₄) + ∂_y(U_y·Ψ₄)] − ∮
−∫Ψ₅⊙N⁵ = +∫ Σ_a (U_a·Ψ₅)⋅∂_a𝖲   →IBP→  −∫ 𝖲⋅[∂_x(U_x·Ψ₅) + ∂_y(U_y·Ψ₅)] + ∮
−∫Ψ₃⊙N³ = +∫ Σ_a (U_a·Ψ₃)⋅∂_a𝖺   (DIRECT — bed Hessian analytic, no IBP)
```

**Slot bookkeeping:** `c=1,2` carry the differentiated divergence in the **k** slot (`Ψ·U`),
`c=3,4,5` in the **j** slot (`U·Ψ`). That k/j asymmetry is why the machine-precision gate uses a
deliberately asymmetric state. The boundary integrals `∮ (g⋅q) n_a` vanish on solid walls
(`u·n = 0`) and behind sponges, and are dropped consistently with every other flux.

---

## 6. Physics selection — three orthogonal controls

`resolve_physics` (`src/problem.jl`) maps three high-level controls onto the seven internal booleans
stored on `BALFEMProblem` (`linearised, advection, lin_pressure, P_full, nl_pressure68,
nl_pressure_full, flat_bed`):

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core, no advection / full nonlinear core + advection |
| `nl_pressure` | `:none` \| `:native` \| `:full` | 𝓝 off / `{3,6,7,8}` / `+{1,2,4,5}` |
| `flat_bed` | `Bool` | `true` ⇔ **`∇h ≡ 0`**; `false` = variable bathymetry |

Pressure content is intrinsic to the model: `P_full = advection`,
`lin_pressure = advection ∨ ¬flat_bed`. `resolve_physics` rejects the footgun
`regime=:linear` with `nl_pressure ≠ :none`; the drivers warn on a `flat_bed` ↔ bathymetry mismatch
(`check_flat_bed_consistency`).

`build_problem` takes the three controls. `build_problem_raw` is the low-level escape hatch taking
the seven booleans, for combinations the high-level interface deliberately does not expose.

### What `flat_bed` selects

**Single point of control:** `dhx, dhy = flat_bed ? (0,0) : (∂ₓh, ∂ᵧh)` in `global_residual`,
`jacobian_u`, `jacobian_u_t`. Then `∇H → ∇η`, `L¹ = 0`, `a = u·∇h = 0`, and every `∇h`-prefixed
block vanishes automatically. Two extra touches: zero the analytic bed Hessian inside
`nlp_native_contrib` (kills `N³`), and skip `nlp_gradh_contrib` when `flat_bed`.

| residual piece | `∇h` dependence | `flat_bed=true` |
|---|---|---|
| leading `P·L` component `L¹ = −u̇·∇h` | explicit `∇h` | → 0 |
| leading `P·L` components `L², L³` | `∇H = ∇η` | survive |
| linear slope `∇h(L·A)` | explicit `∇h` | → 0 |
| linear slope `∇H(L·K)` | `∇H = ∇η` | survives |
| native `𝓝` {6}, {3} (`a = u·∇h`, bed Hessian) | `∇h`, `∇²h` | → 0 |
| native `𝓝` {7,8} | `∇H` / none | survive |
| `{1,2,4,5}` bed-slope IBP half | all `∇h` | skipped entirely |
| `{1,2,4,5}` surface-slope / leading frozen halves | `∇H = ∇η` | survive |
| advection, mass, gravity, dispersion `B` | none | survive |
| gravity's `−η∇h` IBP half | explicit `∇h` | → 0 |

---

## 7. The assembly invariant

> **Every classification row must have exactly ONE consumer, and its guard must be the CONJUNCTION
> of its three activation conditions.**

A guard testing only the bed condition (`flat_bed`) or only the amplitude condition (`regime`) is a
defect whenever the same physics has a second representation in the other regime — which, for the
`𝓐/𝓚` slope-pressure packages, it always does: the linear regime carries the linearised form
`h∇h·𝓛ˡⁱⁿ·(A+K)`, the nonlinear regime the un-expanded `H[∇h(𝓛·A) + ∇H(𝓛·K)]`, and the second
*contains* the first. **They are alternatives, never addends.**

**Corollary: a flat-bed regression can never test `∇h` code.** New bed-slope terms need a
sloping-bed gate from the start.

---

## 8. Jacobians

Hand-written (`jacobian_u`, `jacobian_u_t`). Splitting `∂R/∂u̇` (effective mass = acceleration +
`R_P` dispersion + the `𝓐/𝓚` package) from `∂R/∂u` (spatial block) lets the integrator form its
per-stage system `J = ∂R/∂u + (1/aΔt)∂R/∂u̇` directly.

Coverage differs by regime, **deliberately**:

* **Linear branch: EXACT.** The residual is affine in `(u,u̇)` and every assembled row has its exact
  derivative. Standing consequence and gate: **Newton must converge in ONE iteration per implicit
  stage**, at any amplitude, on any bathymetry (`test/test_linear_newton_gate.jl`).
* **`∂R/∂u̇`: EXACT IN BOTH REGIMES.** Every `u̇`-dependent term is differentiated exactly; the `𝓝`
  blocks carry no `u̇`-dependence. Verified `0.000e+00` against AD in all eight models.
* **Nonlinear `∂R/∂u`: QUASI-NEWTON BY CHOICE.** Advection is differentiated in full; the leading-
  and slope-pressure packages contribute no η-derivative and `𝓝` is absent. Measured gap `2.9e-2`
  (`:none`) to `8.3e-1` (`:full`), **vanishing at order 1.11–1.16 in state amplitude**. Benign
  because it costs Newton iterations, never accuracy — Newton drives the *residual* to zero.

> **The distinction that matters: an omission is benign only if it is HIGHER ORDER IN AMPLITUDE.**
> A block whose prefactor does not scale with the solution (e.g. `H·∇h`) is an `O(1)` error in the
> effective mass matrix; Newton then converges to the fixed point of the **wrong map** and no
> iteration budget rescues it. **Never assume an omission is higher-order — measure it**, with
> `test_jacobians_ad.jl`'s amplitude-scaling gate.
>
> Do not "complete" the quasi-Newton omissions without re-measuring every nonlinear reference value.

---

## 9. Field reconstruction (`src/reconstruct.jl`)

The eliminated vertical kinematics are rebuilt for output:

* `w` is the exact modal vertical-velocity FE. Per σ-level `ℓ` the modal sums collapse to three
  constant-vector contractions, `w_ℓ = −(a_ℓ⋅𝖺) + (b_ℓ⋅𝖻) − (c_ℓ⋅𝖲)` with `a_ℓ = φ(σ_ℓ)`,
  `b_ℓ = σ_ℓφ(σ_ℓ)`, `c_ℓ = φ_int(σ_ℓ)`, reusing the residual's own building blocks — so the code is
  identical sequential and distributed.
* Total pressure `p = ρgH(1−σ)` (hydrostatic) `− ρ Σ_j div(u̇_j) d² Π³_j` (non-hydrostatic, linear
  flat-bed, `u̇` by backward finite difference), where `Π³_j = ∫_σ^1 φ_j_int`. `Π³_j` vanishes at
  `σ=1`, so the reconstructed non-hydrostatic pressure satisfies the free-surface condition
  `p_nh(1) = 0` **exactly** — a free sanity check on any output.
