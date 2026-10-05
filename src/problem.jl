# ==============================================================
#  problem.jl — the BALFEMProblem bundle + the loop-free residual + hand Jacobians
#
#  This file assembles the single scalar residual of §8 of the derivation — the
#  total virtual work ∫_Ω R·v that Gridap's MultiField solver drives to zero at
#  each stage. `BALFEMProblem` gathers the precomputed vertical tensors (as Gridap
#  constants) and the physics flags; `global_residual` evaluates the weak form;
#  `jacobian_u`/`jacobian_u_t` provide the exact spatial and effective-mass
#  Jacobians so Newton needs no finite differencing or automatic differentiation.
#
#  Sign/form conventions (load-bearing — the assembly relies on them):
#    * gravity uses the integrated-by-parts energy form −∫(g/2)(H²−d²)(𝚽⋅DW);
#      subtracting the still-water baseline (H²−d²) makes the discrete rest state
#      exactly force-free, so an undisturbed surface at an open wall stays at rest;
#    * dispersion R_P: −∫ d²(𝗕⋅DUt)⋅DW (linearised) / −∫ H²(𝗕⋅∇·(Hu̇))⋅DW
#      (nonlinear) / the full P¹L¹+P²L²+P³L³ slope decomposition (P_full=true);
#      𝗕 ≤ 0 (the stored dispersion tensor) and the explicit (−1) factors give the
#      term its correct sign — this term IS the frequency dispersion of the model;
#    * the wavemaker enters continuity with a minus, −∫ q·S(x,t);
#    * the pressure slope packages sit on the momentum right-hand side (subtracted).
# ==============================================================

"""
    BALFEMProblem

Coefficient bundle carried through the time loop and consumed by the residual.
It holds the vertical tensors as Gridap constants (built once by
`assemble_vertical_tensors`, reshaped in `build_problem`), the still-water depth
`d(x,y)`, the source/sponge/relaxation profiles, and the physics flags that
switch individual residual terms on or off. Third-order tensors use the index
order `[i,k,j]` = [test layer, u_k layer, u_j layer], so contraction over the
trailing two indices directly yields the mode-i momentum contribution.
"""
struct BALFEMProblem{PV,MV,BV,PT,AT,KT,M3T,G3T,A3T,K3T,P3T}
    g            :: Float64
    h_bathy       :: Function          # still-water depth d(x,y)
    Nσ           :: Int
    Φ            :: PV                # VectorValue{Nσ}   depth weights
    Mv           :: MV                # TensorValue       vertical mass
    Bv           :: BV                # TensorValue       dispersion B ≤ 0 (= −P[:,:,3])
    P            :: PT                # NTuple{3,TensorValue}  leading pressure P^V
    Av           :: AT                # NTuple{3,TensorValue}  linear pressure (∇h)
    Kv           :: KT                # NTuple{3,TensorValue}  linear pressure (∇H)
    M3           :: M3T               # ThirdOrderTensorValue  advection 𝓜
    G3           :: G3T               # ThirdOrderTensorValue  advection 𝓖
    A3           :: A3T               # NTuple{8,ThirdOrderTensorValue}  NL pressure 𝓐
    K3           :: K3T               # NTuple{8,ThirdOrderTensorValue}  NL pressure 𝓚
    P3           :: P3T               # NTuple{8,ThirdOrderTensorValue}  NL leading 𝓟
    linearised   :: Bool              # linear regime: drop H-weights, use d²B dispersion
    advection    :: Bool              # nonlinear advection block
    lin_pressure :: Bool              # A/K linear slope-pressure package (needs a sloped bed)
    P_full       :: Bool              # keep all three slope components P¹L¹+P²L²+P³L³ in R_P
                                      #   (false = the P³L³ dispersion carrier alone)
    nl_pressure  :: Bool              # 𝓝, ALL EIGHT components on, or none. {3,6,7,8} direct; the 𝓐
                                      #   (∇h) half of {1,2,4,5} by exact IBP; the 𝓚/𝓟 halves (Class III)
                                      #   by the BROKEN formulation (src/broken.jl).
    flat_bed     :: Bool              # flat sea-bed assumption ∇h ≡ 0: drops every term carrying a
                                      #   factor ∇h (bed-slope 𝓐 packages, L¹=−u̇·∇h, N{3,6}, the
                                      #   bed-slope IBP half). ∇H = ∇h+∇η → ∇η, so surface-slope
                                      #   (∇η) and dispersion terms are kept. false = variable bathymetry.
    WK3          :: ThirdOrderTensorValue  # reduced Class-III weight for the 𝓚 (surface-slope) block:
    WP3          :: ThirdOrderTensorValue  # …and for the 𝓟 (leading-pressure) block.
                                      #   W[i,k,j] = −T¹[i,k,j] + T²[i,k,j] − T⁵[i,j,k] collapses the
                                      #   THREE ∇𝖲-carrying components {1,2,5} onto ONE contraction.
                                      #   EXACT (verified 4.4e-16) — see NEW_TREATMENT.md §A.2 and
                                      #   `alg_class3_weight`. The 𝓐 (∇h) block is NOT reduced: it is
                                      #   handled by exact IBP and vanishes on a flat bed (§A.5).
    mu_sponge    :: Function
    wm_src       :: Function
    relax_bc     :: Bool              # generation/absorption relaxation zone (Dirichlet inflow)
    relax_mu     :: Function          # zone profile μ_g(x) (quadratic, max at the boundary)
    relax_tg     :: Any               # incident_fields NamedTuple (eta, ux, uy) or nothing
    mms_src      :: Any               # analytic MMS forcing, or `nothing` (default).
                                      #   NamedTuple (Seta, Sx, Sy) of (x,t) callables; Sx/Sy return
                                      #   VectorValue{Nσ}. Subtracted from the residual as
                                      #   F = ∫(q Sη + Wx⋅Sx + Wy⋅Sy), making u* the exact solution
                                      #   of the forced problem. Independent of u ⇒ NO Jacobian
                                      #   contribution. See markdown_files/VERIFIED_SCOPE.md (the analytic-MMS plan itself is gone).
    skel         :: Base.RefValue{Any} # skeleton (interior-facet) context, or `nothing`. ALWAYS present
                                      #   when `nl_pressure` (the broken Class-III layer lives on it; built
                                      #   by `build_problem_raw`); otherwise present only if a stabilisation
                                      #   is attached (`attach_skeleton!`, γ > 0). `nothing` ⇒ the residual
                                      #   has no skeleton term at all. src/broken.jl.
end

"""
    resolve_physics(; regime=:nonlinear, nl_pressure=false, flat_bed=false) → NamedTuple

Translate the high-level physics selection into the internal boolean flags
(`linearised, advection, lin_pressure, P_full, nl_pressure, flat_bed`). This is the single
place the flag couplings are defined and validated:

  * `regime`      — `:linear` (⇒ `linearised`, no advection) or `:nonlinear`
                    (⇒ full nonlinear core with advection);
  * `nl_pressure` — `Bool`: the nonlinear pressure operator 𝓝, all eight components, on or off;
  * `flat_bed`    — the sea-bed geometry: `false` = variable bathymetry (∇h≠0, full
                    model), `true` = flat bed (∇h≡0, every ∇h-term dropped; ∇η-terms kept).

A `Symbol` (the v1 tiers `:none / :native / :full`) is REFUSED rather than translated: a
silently converted `:native` would run the full operator under a name promising a partial one.

Nonlinear pressure is meaningful only in the `:nonlinear` regime, so
`regime=:linear` with `nl_pressure=true` is rejected.
"""
function resolve_physics(; regime::Symbol = :nonlinear,
                             nl_pressure = false,
                             flat_bed::Bool = false)
    regime in (:linear, :nonlinear) ||
        error("resolve_physics: regime must be :linear or :nonlinear (got :$regime)")
    nl_pressure isa Bool ||
        error("resolve_physics: nl_pressure must be a Bool (got $(repr(nl_pressure))). v2 has no " *
              "tiers: true = all eight 𝓝 components, false = none. The v1 symbols :none/:native/:full " *
              "are gone (markdown_files/V2_SOLVER_PLAN.md).")
    regime == :linear && nl_pressure &&
        error("resolve_physics: nl_pressure=true requires regime=:nonlinear " *
              "(a linear model carries no nonlinear pressure)")
    advection = regime == :nonlinear
    return (linearised   = regime == :linear,
            advection    = advection,
            lin_pressure = advection || !flat_bed,   # nonlinear: always; linear: variable-bed only
            P_full       = advection,                # the nonlinear leading pressure is always complete
            nl_pressure  = nl_pressure,
            flat_bed     = flat_bed)
end

"""
    build_problem(vert; g, h_bathy, regime=:nonlinear, nl_pressure=false, flat_bed=false,
                        model=nothing, quad_degree=nothing, mu_sponge, wm_src, relax_*) → BALFEMProblem

Assemble the problem bundle from the high-level physics selection (see
[`resolve_physics`](@ref) for the `regime`/`nl_pressure`/`flat_bed` semantics).
`flat_bed=true` solves the chosen model over a flat sea bed (∇h≡0); `false` over
variable bathymetry. For fine-grained control of the individual boolean flags —
e.g. `lin_pressure` without `P_full` — call [`build_problem_raw`](@ref) directly.

⚠ `nl_pressure=true` REQUIRES `model` and `quad_degree`: the Class-III terms are assembled
by the broken formulation, whose skeleton layer lives on the mesh's interior facets. The
skeleton context is built here, so a full-pressure problem can never exist without it.
"""
function build_problem(vert;
        g            :: Float64  = g,
        h_bathy      :: Function = (x -> 3.5),
        regime       :: Symbol   = :nonlinear,
        nl_pressure              = false,
        flat_bed     :: Bool     = false,
        model                    = nothing,
        quad_degree              = nothing,
        mu_sponge    :: Function = (x -> 0.0),
        wm_src       :: Function = ((x, t) -> 0.0),
        relax_bc     :: Bool     = false,
        relax_mu     :: Function = (x -> 0.0),
        relax_tg                 = nothing,
        mms_src                  = nothing)
    phys = resolve_physics(; regime=regime, nl_pressure=nl_pressure,
                             flat_bed=flat_bed)
    return build_problem_raw(vert; g=g, h_bathy=h_bathy,
        linearised=phys.linearised, advection=phys.advection,
        lin_pressure=phys.lin_pressure, P_full=phys.P_full,
        nl_pressure=phys.nl_pressure, flat_bed=phys.flat_bed,
        model=model, quad_degree=quad_degree,
        mu_sponge=mu_sponge, wm_src=wm_src,
        relax_bc=relax_bc, relax_mu=relax_mu, relax_tg=relax_tg,
        mms_src=mms_src)
end

"""
    build_problem_raw(vert; g, h_bathy, <6 boolean flags>, model, quad_degree, mu_sponge,
                      wm_src, relax_*) → BALFEMProblem

Low-level constructor taking the individual physics booleans directly
(`linearised, advection, lin_pressure, P_full, nl_pressure, flat_bed`).
Use [`build_problem`](@ref) for the ordinary high-level interface; this raw form
exists for the cases that need a flag combination the high-level interface
deliberately does not expose (e.g. `lin_pressure` without `P_full`, used by the
oracle-equivalence test). `flat_bed=false` (default) assembles the bed-slope (∇h)
terms; `true` drops them. Reshapes the `assemble_vertical_tensors` NamedTuple
into constant Gridap tensors and bundles the flags. With `nl_pressure=true` it also builds
the skeleton context (`model`, `quad_degree` required).
"""
function build_problem_raw(vert;
        g            :: Float64  = g,
        h_bathy      :: Function = (x -> 3.5),
        linearised   :: Bool     = false,
        advection    :: Bool     = true,
        lin_pressure :: Bool     = false,
        P_full       :: Bool     = false,
        nl_pressure  :: Bool     = false,
        flat_bed     :: Bool     = false,
        model                    = nothing,
        quad_degree              = nothing,
        mu_sponge    :: Function = (x -> 0.0),
        wm_src       :: Function = ((x, t) -> 0.0),
        relax_bc     :: Bool     = false,
        relax_mu     :: Function = (x -> 0.0),
        relax_tg                 = nothing,
        mms_src                  = nothing)
    Φ  = alg_to_vec(vert.Phi)
    Mv = alg_to_tensor2(vert.Mmat)
    Bv = alg_to_tensor2(vert.B)
    P  = ntuple(c -> alg_to_tensor2(vert.P[:, :, c]), 3)
    Av = ntuple(c -> alg_to_tensor2(vert.A[:, :, c]), 3)
    Kv = ntuple(c -> alg_to_tensor2(vert.K[:, :, c]), 3)
    M3 = alg_to_tensor3(vert.Mcal)
    G3 = alg_to_tensor3(vert.Gcal)
    A3 = ntuple(c -> alg_to_tensor3(vert.Acal[:, :, :, c]), 8)
    K3 = ntuple(c -> alg_to_tensor3(vert.Kcal[:, :, :, c]), 8)
    P3 = ntuple(c -> alg_to_tensor3(vert.Pcal[:, :, :, c]), 8)
    # Reduced Class-III weights for the 𝓚 and 𝓟 blocks. Built from the raw arrays before
    # tensorisation so the (k,j) transpose of component 5 is explicit.
    WK3 = alg_to_tensor3(alg_class3_weight(vert.Kcal[:, :, :, 1],
                                           vert.Kcal[:, :, :, 2],
                                           vert.Kcal[:, :, :, 5]))
    WP3 = alg_to_tensor3(alg_class3_weight(vert.Pcal[:, :, :, 1],
                                           vert.Pcal[:, :, :, 2],
                                           vert.Pcal[:, :, :, 5]))
    relax_bc && relax_tg === nothing &&
        error("build_problem: relax_bc=true requires relax_tg (incident_fields NamedTuple)")
    nl_pressure && linearised &&
        error("build_problem_raw: nl_pressure=true with linearised=true — a linear model carries no 𝓝")
    nl_pressure && (model === nothing || quad_degree === nothing) &&
        error("build_problem: nl_pressure=true needs `model` and `quad_degree` — the Class-III " *
              "terms are assembled by the broken formulation, whose skeleton layer lives on the " *
              "mesh's interior facets (src/broken.jl)")
    prob = BALFEMProblem(g, h_bathy, vert.N_dof, Φ, Mv, Bv, P, Av, Kv, M3, G3,
                         A3, K3, P3, linearised, advection, lin_pressure,
                         P_full, nl_pressure, flat_bed, WK3, WP3,
                         mu_sponge, wm_src,
                         relax_bc, relax_mu, relax_tg, mms_src, Ref{Any}(nothing))
    nl_pressure && (prob.skel[] = skeleton_nt(build_skeleton_ctx(model; degree = quad_degree)))
    return prob
end

"""
    global_residual(t, u, v, prob, trian, dΩh)

Single scalar Gridap residual, stacked layout. `u` is a TransientCellField
(`∂t(u)` available); `u[1]=η`, `u[2]=𝖴x`, `u[3]=𝖴y`. No per-layer loops.
"""
function global_residual(t::Real, u, v, prob::BALFEMProblem, trian, dΩh)
    ut = ∂t(u)
    η,  Ux,  Uy  = u[1], u[2], u[3]
    ηt, Uxt, Uyt = ut[1], ut[2], ut[3]
    q,  Wx,  Wy  = v[1], v[2], v[3]

    g   = prob.g
    lin = prob.linearised
    d_cf   = CellField(prob.h_bathy, trian)
    src_cf = CellField(x -> prob.wm_src(x, t), trian)
    mu_cf  = CellField(prob.mu_sponge, trian)
    H = d_cf + η

    # derived scalar / layer-vector fields.
    # Bed slope ∇h — identically zero under the flat-bed assumption. Zeroing it here
    # is the single point of control for `flat_bed`: it drops L¹=−u̇·∇h, the ∇h
    # prefactor of every bed-slope (𝓐) block, and the a=u·∇h field uniformly, while
    # ∇H = ∇h+∇η collapses to ∇η so the surface-slope (∇η) terms survive.
    dhx = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)
    dhy = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
    dHx = dhx + alg_dx(η);     dHy = dhy + alg_dy(η)           # ∇H = ∇h + ∇η
    DW  = alg_dx(Wx) + alg_dy(Wy)                              # test-divergence vector
    DUt = alg_dx(Uxt) + alg_dy(Uyt)                            # layer div(u̇)
    ub  = alg_vec2(alg_dot(prob.Φ, Ux), alg_dot(prob.Φ, Uy))   # depth-averaged velocity

    # ---- mass continuity + wavemaker source ----------------------------------
    #  LINEARISED regime: still-water depth h(x,y) in the flux, ∇·(h ū) — ∇·(η ū) is O(ε²),
    #  dropped like the H-weighting of acceleration and gravity (LinearModel.tex). Keeping H
    #  here is invisible at wave amplitude but stalls the analytic MMS rate.
    h_cf = lin ? d_cf : H
    r = ∫( q*ηt - h_cf*(∇(q) ⋅ ub) - q*src_cf ) * dΩh

    # ---- acceleration ----------------------------------------------------------
    accx = alg_mul(prob.Mv, Uxt); accy = alg_mul(prob.Mv, Uyt)
    #  LINEAR: h-WEIGHTED (Σⱼ h Mᵢⱼ u̇ⱼ, LinearModel.tex). An h-divided form would be a different
    #  model over variable bathymetry.
    r = lin ? r + ∫( d_cf*((Wx ⋅ accx) + (Wy ⋅ accy)) ) * dΩh :
              r + ∫( H*(Wx ⋅ accx) + H*(Wy ⋅ accy) ) * dΩh

    # ---- gravity (integrated-by-parts energy form; rest-state baseline removed) -
    #  The IBP of the gravity term produces TWO pieces in BOTH regimes, and the
    #  ∇h one is easy to lose because it has no counterpart in the strong form
    #  (`tab: term classification` flags exactly this). Nonlinear identity, EXACT:
    #
    #      ∇((H²−h²)/2) = H∇H − h∇h = H∇η + η∇h    ⇒    H∇η = ∇((H²−h²)/2) − η∇h
    #
    #  so   ∫ g Φᵢ (H∇η)·vᵢ  =  −∫ (g/2)(H²−h²) Φ·(∇·v)  −  ∫ g η Φ·(v·∇h).
    #  The linear case is the same identity with H→h: h∇η = ∇(hη) − η∇h.
    #  The SECOND piece is therefore IDENTICAL in the two regimes — only the first
    #  differs (h η vs (H²−h²)/2) — so it is assembled once, outside the branch.
    #  ⚠ It is needed in BOTH regimes: without it the nonlinear momentum carries an O(∇h)
    #  error that only a variable-bed nonlinear MMS can see.
    PhiDW = alg_dot(prob.Φ, DW)
    r = lin ? r + ∫( (-g)*η*d_cf*PhiDW ) * dΩh :
              r + ∫( (-0.5*g)*(H*H - d_cf*d_cf)*PhiDW ) * dΩh
    #  Shared ∇h half of the IBP. dhx/dhy are already zeroed under flat_bed, so
    #  this vanishes identically there and costs a flat-bed run nothing.
    r = r + ∫( (-g)*η*alg_dot(prob.Φ, dhx*Wx + dhy*Wy) ) * dΩh

    # ---- leading pressure R_P (dispersion — MANDATORY) -------------------------
    if lin
        #  Full linearised leading pressure ∇(h²Σⱼ𝓛ⱼ·Pᵢⱼ), IBP'd onto the test
        #  divergence, with the LINEARISED 𝓛ⱼ = [−u̇ⱼ·∇h, u̇ⱼ·∇h, −∇·(h u̇ⱼ)]
        #  (LinearModel.tex `eq: linearised linear pressure operator`). All three components
        #  collapse to P³L³ on a flat bed.
        LgT = dhx*Uxt + dhy*Uyt                       # u̇ⱼ·∇h, stacked
        Ll1 = (-1.0)*LgT
        Ll2 = LgT
        Ll3 = (-1.0)*(d_cf*DUt + LgT)                 # −∇·(h u̇ⱼ)
        #  L² = −L¹ EXACTLY in the linearised model, so P¹L¹+P²L² = (P¹−P²)L¹. Collapsing
        #  3 contractions to 2 (and 6 to 2 for A/K below) keeps the operator tree shallow.
        Pm  = prob.P[1] - prob.P[2]
        sPl = alg_mul(Pm, Ll1) + alg_mul(prob.P[3], Ll3)
        r = r + ∫( (-1.0)*(d_cf*d_cf)*(sPl ⋅ DW) ) * dΩh
        #  Bed-slope pressure package  −h ∇h Σⱼ𝓛ⱼ·(Aᵢⱼ+Kᵢⱼ). Carries an explicit ∇h,
        #  so it exists only over a non-flat bed — `lin_pressure` is exactly that
        #  condition (resolve_physics: advection || !flat_bed).
        if prob.lin_pressure
            AKm = (prob.Av[1] + prob.Kv[1]) - (prob.Av[2] + prob.Kv[2])
            AK3 =  prob.Av[3] + prob.Kv[3]
            sAK = alg_mul(AKm, Ll1) + alg_mul(AK3, Ll3)
            r = r + ∫( (-1.0)*d_cf*dhx*(sAK ⋅ Wx) ) * dΩh +
                    ∫( (-1.0)*d_cf*dhy*(sAK ⋅ Wy) ) * dΩh
        end
    else
        UgHt = dHx*Uxt + dHy*Uyt                               # u̇ⱼ·∇H stacked
        if prob.P_full
            L1 = (-1.0)*(dhx*Uxt + dhy*Uyt)
            L2 = UgHt
            L3 = (-1.0)*(H*DUt + UgHt)
            sP = alg_mul(prob.P[1], L1) + alg_mul(prob.P[2], L2) + alg_mul(prob.P[3], L3)
            r  = r + ∫( (-1.0)*(H*H)*(sP ⋅ DW) ) * dΩh
        else
            r = r + ∫( (-1.0)*(H*H)*((alg_mul(prob.Bv, H*DUt + UgHt)) ⋅ DW) ) * dΩh
        end
    end

    # ---- sponge (damps velocity AND the free surface η, same μ profile) ---------
    #  +∫ μ q η in continuity gives ηt = … − μη → exponential decay of the surface
    #  in the layer (Newtonian relaxation); essential to absorb the η-dominated
    #  open-boundary mode that a velocity-only sponge leaves under-damped.
    r = r + ∫( mu_cf*( q*η
                     + (Wx ⋅ alg_mul(prob.Mv, Ux))
                     + (Wy ⋅ alg_mul(prob.Mv, Uy)) ) ) * dΩh

    # ---- generation/absorption relaxation zone (Dirichlet inflow) ---------------
    #  Relax the state toward the incident wave in a zone adjacent to the
    #  generation boundary: absorbs outgoing deviations, enforces the incident
    #  field (classical relaxation-zone practice). Linear in u.
    if prob.relax_bc
        tg     = prob.relax_tg
        mug_cf = CellField(prob.relax_mu, trian)
        eta_i  = CellField(x -> tg.eta(x, t), trian)
        ux_i   = CellField(x -> tg.ux(x, t),  trian)
        uy_i   = CellField(x -> tg.uy(x, t),  trian)
        r = r + ∫( mug_cf*q*(η - eta_i)
                 + mug_cf*((Wx ⋅ alg_mul(prob.Mv, Ux - ux_i))
                         + (Wy ⋅ alg_mul(prob.Mv, Uy - uy_i))) ) * dΩh
    end

    # ---- nonlinear advection (𝓜/𝓖 block) ---------------------------------------
    if prob.advection
        DU  = alg_dx(Ux) + alg_dy(Uy)
        UgH = dHx*Ux + dHy*Uy
        S   = H*DU + UgH                                       # ∇·(H u_j) stacked
        TMx = alg_outer(Ux, alg_dx(Ux)) + alg_outer(Uy, alg_dy(Ux))
        TMy = alg_outer(Ux, alg_dx(Uy)) + alg_outer(Uy, alg_dy(Uy))
        advx = alg_dc3(prob.M3, TMx);  advy = alg_dc3(prob.M3, TMy)
        gvx  = alg_dc3(prob.G3, alg_outer(S, Ux))
        gvy  = alg_dc3(prob.G3, alg_outer(S, Uy))
        r = r + ∫( H*(advx ⋅ Wx) + H*(advy ⋅ Wy) + (gvx ⋅ Wx) + (gvy ⋅ Wy) ) * dΩh
    end

    # ---- linear non-hydrostatic pressure (A/K slope package) --------------------
    #  The NONLINEAR representation of the 𝓐/𝓚 slope package, H[∇h(𝓛·A) + ∇H(𝓛·K)], which
    #  CONTAINS the linearised `sAK` block above. ⚠ The two are ALTERNATIVES, never addends, so
    #  the guard is the CONJUNCTION with `!lin` (`lin_pressure` alone is also true in the linear
    #  variable-bed case). Assembly invariant: CLAUDE.md rule 4, MODEL.md §7.
    if !lin && prob.lin_pressure
        UgHt = dHx*Uxt + dHy*Uyt
        L1 = (-1.0)*(dhx*Uxt + dhy*Uyt)
        L2 = UgHt
        L3 = (-1.0)*(H*DUt + UgHt)
        LA = alg_mul(prob.Av[1], L1) + alg_mul(prob.Av[2], L2) + alg_mul(prob.Av[3], L3)
        LK = alg_mul(prob.Kv[1], L1) + alg_mul(prob.Kv[2], L2) + alg_mul(prob.Kv[3], L3)
        r = r + ∫( (-1.0)*H*( dhx*(Wx ⋅ LA) + dHx*(Wx ⋅ LK)
                            + dhy*(Wy ⋅ LA) + dHy*(Wy ⋅ LK) ) ) * dΩh
    end

    # ---- nonlinear pressure 𝓝 — all eight components (nlpressure.jl, broken.jl) -----
    #  {3,6,7,8}: first-order everywhere, assembled directly in all three blocks.
    #  {1,2,4,5}, bed-slope 𝓐 half: EXACT integration by parts onto the test (pure ∇h;
    #             skipped on a flat bed).
    #  {1,2,4,5}, surface-slope 𝓚 and leading 𝓟 halves (Class III): the BROKEN formulation —
    #             the distributional gradients of 𝖲 and 𝖻, i.e. cellwise Hessians + the skeleton
    #             layer. Both the ∇𝖲 and the ∇𝖻 arm, always.
    if prob.nl_pressure
        DU  = alg_dx(Ux) + alg_dy(Uy)
        UgH = dHx*Ux + dHy*Uy
        Ugh = dhx*Ux + dhy*Uy
        S   = H*DU + UgH
        r = r + nlp_direct_contrib(prob, d_cf, η, H, dhx, dhy, dHx, dHy,
                                   Ux, Uy, Wx, Wy, DW, Ugh, UgH, S, DU, dΩh)
        prob.flat_bed || (r = r + nlp_gradh_contrib(prob, d_cf, η, H, dhx, dhy,
                                  Ux, Uy, Wx, Wy, Ugh, UgH, S, DU, dΩh))
        r = r + broken_class3_residual(prob, prob.skel[], d_cf, η, H, dHx, dHy,
                                       Ux, Uy, DU, S, UgH, Wx, Wy, DW, dΩh)
    end

    # ---- skeleton stabilisation (src/broken.jl): :jumpgrad or :ghostvolume, any model ------
    #  Sign-definite damping of the inter-element mismatch; consistent (weakly). Absent unless a
    #  γ > 0 is attached (`attach_skeleton!`).
    has_cip(prob) && (r = r + stab_contrib(prob, prob.skel[], η, Ux, Uy, q, Wx, Wy))

    # ---- analytic MMS forcing (verification only; `nothing` in every physical run) --
    #  Subtract F(t;q,vᵢ) = ∫(q Sη + Σᵢ Sᵢ·vᵢ) so that the manufactured field u* is the
    #  EXACT solution of R − F = 0. The forcing is derived from the governing equations
    #  in closed form (src/mms.jl) and never calls this file — that independence is the
    #  whole point: it is what lets the measured convergence RATE detect a residual that
    #  is self-consistently wrong. `t` here is the value the integrator passes (RK stage
    #  time, or t_{n+θ}), which is exactly the "forcing exact in time" requirement.
    #  No Jacobian contribution: F is independent of u and u̇.
    if prob.mms_src !== nothing
        S   = prob.mms_src
        Sη  = CellField(x -> S.Seta(x, t), trian)
        Sxf = CellField(x -> S.Sx(x, t),   trian)
        Syf = CellField(x -> S.Sy(x, t),   trian)
        r = r - ∫( q*Sη + (Wx ⋅ Sxf) + (Wy ⋅ Syf) ) * dΩh
    end

    return r
end

# ----------------------------------------------------------
#  Hand Jacobians. Splitting ∂R/∂u̇ (effective mass) from ∂R/∂u (spatial) lets
#  the time integrator form its per-stage system J = ∂R/∂u + (1/aΔt)∂R/∂u̇
#  directly.
#
#  COVERAGE — differs by regime, deliberately (MODEL.md §6–7):
#
#   * LINEAR branch (`prob.linearised`): the Jacobians are EXACT. The residual is
#     affine in (u,u̇) there, so every assembled row has its exact derivative here —
#     C1/M1/M10/M14 in ∂R/∂u̇, C2,C3/M3/sponge/relax in ∂R/∂u. Consequence, and the
#     standing gate that protects it: Newton MUST converge in ONE iteration for any
#     linear configuration, at any amplitude, on any bathymetry
#     (test/test_linear_newton_gate.jl). Any excess iteration is proof of a
#     residual↔Jacobian inconsistency.
#
#   * NONLINEAR branch: ∂R/∂u̇ is EXACT; ∂R/∂u is QUASI-NEWTON by choice.
#
#     ∂R/∂u̇ (jacobian_u_t): every u̇-dependent term is differentiated exactly — mass, the
#     H-weighted acceleration, the leading pressure R_P and the 𝓐/𝓚 slope-pressure package.
#     The 𝓝 blocks carry no u̇-dependence.
#
#     ∂R/∂u (jacobian_u): still quasi-Newton, deliberately. Advection is differentiated
#     in full (so Newton stays quadratic on the dominant nonlinearity), and so are the
#     Class-III 𝓚/𝓟 blocks (broken_class3_jacobian); the leading- and slope-pressure
#     packages contribute no η-derivative (their dependence through H and ∇H is frozen) and
#     the other 𝓝 blocks ({3,6,7,8}, the ∇h IBP half) add to the residual but not here.
#
#     An omission is benign only if it is HIGHER ORDER IN AMPLITUDE (it costs iterations,
#     never accuracy); one with an O(1) prefactor (e.g. H·∇h) makes Newton converge to the
#     wrong map. Measure it with test_jacobians_ad.jl (rule 5); never assume it.
#
#     Do not extend the remaining ∂R/∂u omissions without re-measuring every nonlinear
#     reference value.
# ----------------------------------------------------------

"∂R/∂u̇ — effective mass operator (acceleration + R_P dispersion)."
function jacobian_u_t(t::Real, u, dut, v, prob::BALFEMProblem, trian, dΩh)
    η = u[1]
    dηt, dUxt, dUyt = dut[1], dut[2], dut[3]
    q, Wx, Wy = v[1], v[2], v[3]
    lin  = prob.linearised
    d_cf = CellField(prob.h_bathy, trian)
    H = d_cf + η

    r = ∫( q*dηt ) * dΩh
    accx = alg_mul(prob.Mv, dUxt); accy = alg_mul(prob.Mv, dUyt)
    #  h-weighted in the linear regime, matching global_residual (MMS_VARBED_PLAN §0.A)
    r = lin ? r + ∫( d_cf*((Wx ⋅ accx) + (Wy ⋅ accy)) ) * dΩh :
              r + ∫( H*(Wx ⋅ accx) + H*(Wy ⋅ accy) ) * dΩh

    DW   = alg_dx(Wx) + alg_dy(Wy)
    dDUt = alg_dx(dUxt) + alg_dy(dUyt)
    if lin
        #  Exact derivative of the full linearised leading pressure AND the bed-slope
        #  package; both are linear in u̇, so both belong here (MMS_VARBED_PLAN §0.A).
        dhxl = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)
        dhyl = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
        dLgT = dhxl*dUxt + dhyl*dUyt
        dLl1 = (-1.0)*dLgT
        dLl2 = dLgT
        dLl3 = (-1.0)*(d_cf*dDUt + dLgT)
        dPm  = prob.P[1] - prob.P[2]
        dsPl = alg_mul(dPm, dLl1) + alg_mul(prob.P[3], dLl3)
        r = r + ∫( (-1.0)*(d_cf*d_cf)*(dsPl ⋅ DW) ) * dΩh
        if prob.lin_pressure
            dAKm = (prob.Av[1] + prob.Kv[1]) - (prob.Av[2] + prob.Kv[2])
            dAK3 =  prob.Av[3] + prob.Kv[3]
            dsAK = alg_mul(dAKm, dLl1) + alg_mul(dAK3, dLl3)
            r = r + ∫( (-1.0)*d_cf*dhxl*(dsAK ⋅ Wx) ) * dΩh +
                    ∫( (-1.0)*d_cf*dhyl*(dsAK ⋅ Wy) ) * dΩh
        end
    else
        dhx = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)   # ∇h ≡ 0 on a flat bed
        dhy = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
        dHx = dhx + alg_dx(η); dHy = dhy + alg_dy(η)
        dUgHt = dHx*dUxt + dHy*dUyt
        #  𝓛 differentiated w.r.t. u̇. Needed by BOTH the leading-pressure block and
        #  the 𝓐/𝓚 slope package below, so it is computed once here rather than
        #  inside `if P_full` — `build_problem_raw` can set lin_pressure WITHOUT
        #  P_full (the oracle-equivalence split), and that combination needs these.
        dL1 = (-1.0)*(dhx*dUxt + dhy*dUyt)
        dL2 = dUgHt
        dL3 = (-1.0)*(H*dDUt + dUgHt)
        if prob.P_full
            sP  = alg_mul(prob.P[1], dL1) + alg_mul(prob.P[2], dL2) + alg_mul(prob.P[3], dL3)
            r   = r + ∫( (-1.0)*(H*H)*(sP ⋅ DW) ) * dΩh
        else
            r = r + ∫( (-1.0)*(H*H)*((alg_mul(prob.Bv, H*dDUt + dUgHt)) ⋅ DW) ) * dΩh
        end
        #  ---- 𝓐/𝓚 slope-pressure package ----------------------------------------
        #  MANDATORY here: its prefactor H·∇h does not scale with the solution, so leaving it
        #  out is an O(1) error in the effective mass and Newton converges to the wrong map
        #  (rule 5). It is EXACTLY LINEAR IN u̇, so its derivative is the residual expression
        #  with u̇ → du̇, in the same 6-contraction form (L2 ≠ −L1 in the nonlinear branch).
        if prob.lin_pressure
            dLA = alg_mul(prob.Av[1], dL1) + alg_mul(prob.Av[2], dL2) + alg_mul(prob.Av[3], dL3)
            dLK = alg_mul(prob.Kv[1], dL1) + alg_mul(prob.Kv[2], dL2) + alg_mul(prob.Kv[3], dL3)
            r = r + ∫( (-1.0)*H*( dhx*(Wx ⋅ dLA) + dHx*(Wx ⋅ dLK)
                                + dhy*(Wy ⋅ dLA) + dHy*(Wy ⋅ dLK) ) ) * dΩh
        end
    end
    return r
end

"∂R/∂u — continuity + gravity + sponge + (nonlinear Acc η-term) + FULL advection derivative."
function jacobian_u(t::Real, u, du, v, prob::BALFEMProblem, trian, dΩh)
    η,  Ux,  Uy  = u[1], u[2], u[3]
    dη, dUx, dUy = du[1], du[2], du[3]
    q,  Wx,  Wy  = v[1], v[2], v[3]
    g = prob.g; lin = prob.linearised
    d_cf  = CellField(prob.h_bathy, trian)
    mu_cf = CellField(prob.mu_sponge, trian)
    H = d_cf + η

    ub  = alg_vec2(alg_dot(prob.Φ, Ux),  alg_dot(prob.Φ, Uy))
    dub = alg_vec2(alg_dot(prob.Φ, dUx), alg_dot(prob.Φ, dUy))

    # continuity — the ∂/∂η term exists only in the nonlinear flux ∇·(H ū); the
    # linearised flux ∇·(h ū) carries no η dependence at all (see global_residual).
    r = lin ? ∫( (-1.0)*( d_cf*(∇(q) ⋅ dub) ) ) * dΩh :
              ∫( (-1.0)*( (∇(q) ⋅ ub)*dη + H*(∇(q) ⋅ dub) ) ) * dΩh

    # gravity
    DW = alg_dx(Wx) + alg_dy(Wy)
    PhiDW = alg_dot(prob.Φ, DW)
    #  Exact derivative of the IBP gravity form. Both regimes carry the shared ∇h
    #  half (see the derivation at the gravity block in global_residual); only the
    #  ∇·v half differs, h η vs (H²−h²)/2, whose η-derivatives are h and H.
    dhxg = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)
    dhyg = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
    r = lin ? r + ∫( (-g)*dη*d_cf*PhiDW ) * dΩh :
              r + ∫( (-g)*H*dη*PhiDW ) * dΩh
    #  ∂/∂η of  −∫ g η Φ·(v·∇h)  — the shared ∇h half. Present in BOTH regimes,
    #  matching the residual; it was previously in the linear branch only.
    r = r + ∫( (-g)*dη*alg_dot(prob.Φ, dhxg*Wx + dhyg*Wy) ) * dΩh

    # nonlinear-Acc η derivative (needs current u̇)
    if !lin
        ut = ∂t(u)
        r = r + ∫( dη*(Wx ⋅ alg_mul(prob.Mv, ut[2])) + dη*(Wy ⋅ alg_mul(prob.Mv, ut[3])) ) * dΩh
    end

    # sponge (velocity + free-surface η damping; exact linear derivative)
    r = r + ∫( mu_cf*( q*dη
                     + (Wx ⋅ alg_mul(prob.Mv, dUx))
                     + (Wy ⋅ alg_mul(prob.Mv, dUy)) ) ) * dΩh

    # relaxation zone (linear in u — exact derivative)
    if prob.relax_bc
        mug_cf = CellField(prob.relax_mu, trian)
        r = r + ∫( mug_cf*q*dη
                 + mug_cf*((Wx ⋅ alg_mul(prob.Mv, dUx))
                         + (Wy ⋅ alg_mul(prob.Mv, dUy))) ) * dΩh
    end

    # full advection derivative
    if prob.advection
        dhx = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)   # ∇h ≡ 0 on a flat bed
        dhy = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
        dHx = dhx + alg_dx(η); dHy = dhy + alg_dy(η)
        DU  = alg_dx(Ux) + alg_dy(Uy)
        UgH = dHx*Ux + dHy*Uy
        S   = H*DU + UgH
        dDU = alg_dx(dUx) + alg_dy(dUy)
        dS  = dη*DU + H*dDU + (alg_dx(dη)*Ux + alg_dy(dη)*Uy) + (dHx*dUx + dHy*dUy)

        TMx  = alg_outer(Ux, alg_dx(Ux)) + alg_outer(Uy, alg_dy(Ux))
        TMy  = alg_outer(Ux, alg_dx(Uy)) + alg_outer(Uy, alg_dy(Uy))
        dTMx = alg_outer(dUx, alg_dx(Ux)) + alg_outer(Ux, alg_dx(dUx)) +
               alg_outer(dUy, alg_dy(Ux)) + alg_outer(Uy, alg_dy(dUx))
        dTMy = alg_outer(dUx, alg_dx(Uy)) + alg_outer(Ux, alg_dx(dUy)) +
               alg_outer(dUy, alg_dy(Uy)) + alg_outer(Uy, alg_dy(dUy))
        dTGx = alg_outer(dS, Ux) + alg_outer(S, dUx)
        dTGy = alg_outer(dS, Uy) + alg_outer(S, dUy)

        r = r + ∫( dη*((alg_dc3(prob.M3, TMx)) ⋅ Wx) + dη*((alg_dc3(prob.M3, TMy)) ⋅ Wy)
                 + H*((alg_dc3(prob.M3, dTMx)) ⋅ Wx) + H*((alg_dc3(prob.M3, dTMy)) ⋅ Wy)
                 + ((alg_dc3(prob.G3, dTGx)) ⋅ Wx) + ((alg_dc3(prob.G3, dTGy)) ⋅ Wy) ) * dΩh
    end

    # Skeleton stabilisation: LINEAR in (η,𝖴) ⇒ its exact derivative is the same form on (dη,d𝖴).
    has_cip(prob) && (r = r + stab_contrib(prob, prob.skel[], dη, dUx, dUy, q, Wx, Wy))

    # Class-III 𝓚/𝓟 blocks (broken: volume + skeleton layer): EXACT linearisation. They depend on
    # the iterate through cellwise Hessians and facet jumps that scale like 1/h; leaving them out
    # stalled Newton at 32 cells/λ (GHOST_PENALTY_PLAN.md §5.5). The other 𝓝 blocks ({3,6,7,8} and
    # the ∇h IBP half) stay quasi-Newton — O(A²) and benign (rule 5). src/broken.jl.
    if prob.nl_pressure
        r = r + broken_class3_jacobian(prob, prob.skel[], d_cf, η, Ux, Uy, dη, dUx, dUy, Wx, Wy, dΩh)
    end

    return r
end

"TransientFEOperator with the hand Jacobians (default; fast)."
function build_ode_operator(prob::BALFEMProblem, U, V, trian, dΩh)
    r  = (t, u, v)      -> global_residual(t, u, v, prob, trian, dΩh)
    j  = (t, u, du, v)  -> jacobian_u(t, u, du, v, prob, trian, dΩh)
    jt = (t, u, dut, v) -> jacobian_u_t(t, u, dut, v, prob, trian, dΩh)
    return TransientFEOperator(r, j, jt, U, V)
end

"Operator variant that lets Gridap build the Jacobians by automatic
differentiation of `global_residual` instead of using the hand Jacobians —
useful for cross-checking. The three-field stacked residual keeps the AD
expression tree small enough to be practical."
function build_ode_operator_ad(prob::BALFEMProblem, U, V, trian, dΩh)
    r = (t, u, v) -> global_residual(t, u, v, prob, trian, dΩh)
    return TransientFEOperator(r, U, V)
end
