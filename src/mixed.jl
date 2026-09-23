# ==============================================================
#  mixed.jl — THE MIXED (projection-free) FORMULATION of the Class-III blocks
#
#  WHY THIS EXISTS. It is a DIAGNOSTIC, not a production path. The `:full` tier is
#  unstable on a flat bed (OPEN_ISSUES.md §0c) and the leading suspect is the frozen
#  L²-projection treatment of the Class-III components {1,2,4,5}. This module removes
#  the projection entirely so that suspicion can be TESTED rather than argued:
#
#      stable here   ⇒ strong evidence the projection treatment is the culprit
#      unstable here ⇒ the projection is exonerated; look elsewhere
#
#  ⚠ THE INFERENCE IS ASYMMETRIC. A mixed formulation carries its OWN inf-sup condition
#  — adding 𝖦 creates a saddle-point structure whose stability depends on how the 𝖦
#  space pairs with the velocity space. A STABLE result is therefore clean evidence; an
#  UNSTABLE one is ambiguous between "the projection was not the cause" and "this 𝖦
#  space is a bad choice". 𝖦 and 𝖥 are taken in the VELOCITY space here, and that
#  choice must be quoted with any negative result.
#
#  ⚠ AND IT IS *CONSISTENT*, NOT *EXACT*. 𝖦 is still a finite-element approximation of
#  ∇𝖲, living in its own space. What is removed is the ad-hoc recovery (project, then
#  differentiate) and the one-step lag, replaced by a variationally consistent,
#  fully-coupled definition that enters the Jacobian. Do not claim "no approximation".
#
#  THE FORMULATION. Two auxiliary vector unknowns per horizontal direction:
#
#      𝖦x ≈ ∂ₓ𝖲 ,  𝖦y ≈ ∂ᵧ𝖲       with  𝖲 = ∇·(H𝗎)
#      𝖥x ≈ ∂ₓ𝖻 ,  𝖥y ≈ ∂ᵧ𝖻       with  𝖻 = 𝗎·∇H
#
#  defined weakly WITH THE INTEGRATION BY PARTS BUILT IN, so no derivative ever lands
#  on 𝖲 or 𝖻 (this is the whole point — every unknown carries at most one derivative):
#
#      ∫ 𝖦ₐ·Ψ  +  ∫ 𝖲·(∂ₐΨ)  =  0        ∀Ψ
#      ∫ 𝖥ₐ·Ψ  +  ∫ 𝖻·(∂ₐΨ)  =  0        ∀Ψ
#
#  ⚠ THE BOUNDARY TERM ∮ 𝖲 Ψ nₐ IS ASSEMBLED, NOT DROPPED — AND THAT IS LOAD-BEARING.
#  The bed-slope IBP (`nlp_gradh_contrib`) may drop its boundary term because 𝗎·n = 0
#  enters every one of its integrands. HERE IT MAY NOT: 𝖲 = ∇·(H𝗎) is generally NON-ZERO
#  on a solid wall even when 𝗎·n = 0, so the same argument does not transfer.
#  MEASURED (test_mixed_formulation.jl G4, against an analytic ∂ₓ𝖲): dropping it leaves
#  the interior recovery essentially exact at 2.7e-06 relative, while the BOUNDARY cells
#  are wrong by a relative 22.4 — a factor of twenty-two, not a perturbation.
#  ⚠ THAT IS EXACTLY WHERE THE INSTABILITY UNDER INVESTIGATION LIVES. The `:full` mode is
#  velocity-led AT THE GENERATION BOUNDARY (CLAUDE.md §5.2b: x = 0.50 m, the inflow). A
#  diagnostic that is locally wrong by 22x precisely there could not have distinguished
#  "the projection is the culprit" from "my boundary treatment is the culprit".
#
#  MultiField:  [η, 𝖴x, 𝖴y]  →  [η, 𝖴x, 𝖴y, 𝖦x, 𝖦y, 𝖥x, 𝖥y]     (3 → 7 fields)
#  The physics indices are UNCHANGED, so `global_residual` needs no reindexing.
#
#  ⚠ 𝖥 (≈ ∇𝖻, component 4) IS OPTIONAL, AND BY DEFAULT IT IS GONE. The layout follows
#  `prob.c3_mask`:
#      c3_mask = (true, true)   → 7 fields  [η, 𝖴x, 𝖴y, 𝖦x, 𝖦y, 𝖥x, 𝖥y]
#      c3_mask = (true, false)  → 5 fields  [η, 𝖴x, 𝖴y, 𝖦x, 𝖦y]        ← the diagnostic
#  MEASURED JUSTIFICATION, not a guess:
#    * STATICALLY (test_class3_split.jl, assembled residual, flat bed, Q3/Q2):
#        |∇𝖲 arm| = 1.2837e-02   |∇𝖻 arm| = 1.1271e-05   — a ratio of ~1140x.
#    * DYNAMICALLY (run/local/run_1dc3iso_*.sh, the 60 m flume that dies at t = 12.60 s):
#        ∇𝖲 ALONE reproduces the failure — t = 12.00 s, u_max 0.61 → 3.17 in one output
#        interval, x_at_max pinned at 0.25 m (the inflow), Newton → cap.
#        ∇𝖻 ALONE ran the full 25 periods to t = 40 s with Newton flat at 4.
#  So component 4 is neither large nor destabilising, and carrying it here would double
#  the auxiliary unknowns and add a second untested space pairing for nothing.
#  ⚠ CONSEQUENCE FOR THE COMPARISON: the 5-field path is the projection-free version of
#  the `gs` ARM, so its control is the PROJECTED gs run (t = 12.00 s), not plain `:full`.
#
#  JACOBIANS: AD ONLY. `jacobian_u`/`jacobian_u_t` are hand-derived for the 3-field
#  layout and do not know about 𝖦/𝖥. `build_ode_operator_mixed` therefore builds the
#  AD operator, which also removes any quasi-Newton ambiguity from the result.
# ==============================================================

"""
    global_residual_mixed(t, u, v, prob, trian, dΩh)

Residual for the 7-field mixed layout. Delegates the physics rows to
`global_residual` with `aux` supplied — so the momentum/continuity equations are the
SAME code the 3-field path runs — and adds the four auxiliary weak definitions.

⚠ The auxiliary rows carry NO time derivative: they are algebraic constraints, making
the system an index-1 DAE. `∂R/∂u̇` has zero rows there; the stage system stays
non-singular because `∂R_aux/∂aux` is the auxiliary mass matrix.
"""
function global_residual_mixed(t::Real, u, v, prob::BALFEMProblem, trian, dΩh, dΓ, nΓ)
    with_F = prob.c3_mask[2]
    aux = with_F ? (Gx = u[4], Gy = u[5], Fx = u[6], Fy = u[7]) :
                   (Gx = u[4], Gy = u[5])
    r   = global_residual(t, u, v, prob, trian, dΩh; aux = aux)

    η, Ux, Uy = u[1], u[2], u[3]
    Pgx, Pgy = v[4], v[5]

    d_cf = CellField(prob.h_bathy, trian)
    H    = d_cf + η
    # flat_bed ⇔ ∇h ≡ 0 — the SINGLE control point, matching `global_residual` exactly.
    dhx  = prob.flat_bed ? 0.0 * d_cf : alg_dx(d_cf)
    dhy  = prob.flat_bed ? 0.0 * d_cf : alg_dy(d_cf)
    dHx  = dhx + alg_dx(η)
    dHy  = dhy + alg_dy(η)
    DU   = alg_dx(Ux) + alg_dy(Uy)
    bf   = dHx * Ux + dHy * Uy          # 𝖻 = 𝗎·∇H
    S    = H * DU + bf                  # 𝖲 = ∇·(H𝗎) = H∇·𝗎 + 𝖻

    nx = Operation(n -> n ⋅ Ex)(nΓ)
    ny = Operation(n -> n ⋅ Ey)(nΓ)

    # 𝖦 ≈ ∇𝖲 :  ∫ 𝖦ₐ·Ψ + ∫ 𝖲·(∂ₐΨ) − ∮ 𝖲 Ψ nₐ = 0
    r = r + ∫( (u[4] ⋅ Pgx) + (S ⋅ alg_dx(Pgx)) ) * dΩh - ∫( (S ⋅ Pgx) * nx ) * dΓ
    r = r + ∫( (u[5] ⋅ Pgy) + (S ⋅ alg_dy(Pgy)) ) * dΩh - ∫( (S ⋅ Pgy) * ny ) * dΓ

    # 𝖥 ≈ ∇𝖻 : only when component 4 is active (see the header — it usually is not)
    if with_F
        Pfx, Pfy = v[6], v[7]
        r = r + ∫( (u[6] ⋅ Pfx) + (bf ⋅ alg_dx(Pfx)) ) * dΩh - ∫( (bf ⋅ Pfx) * nx ) * dΓ
        r = r + ∫( (u[7] ⋅ Pfy) + (bf ⋅ alg_dy(Pfy)) ) * dΩh - ∫( (bf ⋅ Pfy) * ny ) * dΓ
    end
    return r
end

"How many auxiliary VectorValue{Nσ} fields the mixed layout needs for this problem."
mixed_n_aux(prob::BALFEMProblem) = prob.c3_mask[2] ? 4 : 2

"""
    mixed_delta_S(prob, u, du, trian) -> (δ𝖲, δ𝖻)

Directional derivative of `𝖲 = ∇·(H𝗎)` and `𝖻 = 𝗎·∇H` in the direction `du = (dη,d𝖴x,d𝖴y)`.

    δ𝖲 = ∇·(dη·𝗎) + ∇·(H·d𝗎) = dη(∇·𝗎) + ∇dη·𝗎 + H(∇·d𝗎) + ∇H·d𝗎
    δ𝖻 = ∇dη·𝗎 + ∇H·d𝗎

⚠ EVERY TERM IS A PRODUCT OF FIRST DERIVATIVES, so this introduces no admissibility problem:
`𝖲` was always first-order — it is `∇𝖲` that is not, which is why `𝖦` exists in the first place.
"""
function mixed_delta_S(prob::BALFEMProblem, u, du, trian)
    η,  Ux,  Uy  = u[1],  u[2],  u[3]
    dη, dUx, dUy = du[1], du[2], du[3]
    d_cf = CellField(prob.h_bathy, trian)
    H    = d_cf + η
    dhx  = prob.flat_bed ? 0.0 * d_cf : alg_dx(d_cf)
    dhy  = prob.flat_bed ? 0.0 * d_cf : alg_dy(d_cf)
    dHx  = dhx + alg_dx(η);   dHy = dhy + alg_dy(η)
    DU   = alg_dx(Ux) + alg_dy(Uy)
    dDU  = alg_dx(dUx) + alg_dy(dUy)
    δb   = alg_dx(dη)*Ux + alg_dy(dη)*Uy + dHx*dUx + dHy*dUy
    δS   = dη*DU + H*dDU + δb
    return δS, δb
end

"""
    mixed_coupling_jacobian(prob, u, du, v, trian, dΩh, dΓ, nΓ; coupling=true)

The auxiliary rows and their couplings in the mixed Jacobian (NEW_TREATMENT.md Part G):

    M  = ∂R_𝖦/∂𝖦     the auxiliary mass matrix  (always)
    C  = ∂R_𝖦/∂(η,𝗎)  the 𝖲-dependence of the auxiliary rows   — `coupling=true`
    B  = ∂R_phys/∂𝖦   the Class-III blocks' 𝖦-dependence        — `coupling=true`

⚠ `C` CARRIES THE BOUNDARY TERM, because the residual does. An auxiliary Jacobian without it
would be inconsistent with its own residual at the boundary — reintroducing exactly the defect
that measured a relative error of 22.4 in the boundary cells (§F.1).

⚠ `B` is trivial only because `GU = Σₐ 𝖦ₐ⊗𝖴ₐ` is LINEAR in `𝖦`: `δGU = Σₐ d𝖦ₐ⊗𝖴ₐ`, with `SD` and
`N4` untouched since neither depends on `𝖦`. It is pushed through the SAME reduced contributors the
residual uses, so the two cannot drift apart.

With `coupling=false` this degenerates to the original block-diagonal form, kept so the two can be
compared directly.
"""
function mixed_coupling_jacobian(prob::BALFEMProblem, u, du, v, trian, dΩh, dΓ, nΓ;
                                 coupling::Bool = true)
    with_F = prob.c3_mask[2]
    Pgx, Pgy = v[4], v[5]
    # ---- M : the auxiliary mass block ---------------------------------------
    r = ∫( (du[4] ⋅ Pgx) + (du[5] ⋅ Pgy) ) * dΩh
    with_F && (r = r + ∫( (du[6] ⋅ v[6]) + (du[7] ⋅ v[7]) ) * dΩh)
    coupling || return r

    nx = Operation(n -> n ⋅ Ex)(nΓ)
    ny = Operation(n -> n ⋅ Ey)(nΓ)
    δS, δb = mixed_delta_S(prob, u, du, trian)

    # ---- C : ∂R_𝖦/∂(η,𝗎), volume + boundary, matching the residual exactly ----
    r = r + ∫( δS ⋅ alg_dx(Pgx) ) * dΩh - ∫( (δS ⋅ Pgx) * nx ) * dΓ
    r = r + ∫( δS ⋅ alg_dy(Pgy) ) * dΩh - ∫( (δS ⋅ Pgy) * ny ) * dΓ
    if with_F
        r = r + ∫( δb ⋅ alg_dx(v[6]) ) * dΩh - ∫( (δb ⋅ v[6]) * nx ) * dΓ
        r = r + ∫( δb ⋅ alg_dy(v[7]) ) * dΩh - ∫( (δb ⋅ v[7]) * ny ) * dΓ
    end

    # ---- B : ∂R_phys/∂𝖦, through the SAME reduced contributors ---------------
    if any(prob.c3_mask)
        η, Ux, Uy = u[1], u[2], u[3]
        q, Wx, Wy = v[1], v[2], v[3]
        DW   = alg_dx(Wx) + alg_dy(Wy)
        d_cf = CellField(prob.h_bathy, trian)
        H    = d_cf + η
        dhx  = prob.flat_bed ? 0.0 * d_cf : alg_dx(d_cf)
        dhy  = prob.flat_bed ? 0.0 * d_cf : alg_dy(d_cf)
        dHx  = dhx + alg_dx(η);  dHy = dhy + alg_dy(η)
        #  δGU from 𝖦 alone; SD and N4 do not depend on 𝖦, so they are ZERO here.
        δGU  = alg_outer(du[4], Ux) + alg_outer(du[5], Uy)
        Z    = 0.0 * δGU
        δN4  = with_F ? (alg_outer(Ux, du[6]) + alg_outer(Uy, du[7])) : Z
        r = r + nlp_gradH_reduced_contrib(prob, H, dHx, dHy, Wx, Wy, δGU, Z, δN4, dΩh)
        r = r + nlp_P_reduced_contrib(prob, H, DW, δGU, Z, δN4, dΩh)
    end
    return r
end

"""
    build_ode_operator_mixed(prob, U, V, trian, dΩh)

`TransientFEOperator` for the mixed layout, with **AD Jacobians** (the hand Jacobians
are 3-field only). Slower per assembly and that is accepted: this path exists to answer
a stability question, not to run production.
"""
function build_ode_operator_mixed(prob::BALFEMProblem, U, V, trian, dΩh; bdeg::Int = 10,
                                  use_ad::Bool = false, coupling::Bool = true)
    #  ⚠ The FE layout and `c3_mask` must agree, or u[4]/u[5] silently index the wrong
    #  field. `build_fe_spaces(...; n_aux = mixed_n_aux(prob))` is the only correct call.
    nfields = length(U.spaces)
    nfields == 3 + mixed_n_aux(prob) ||
        error("build_ode_operator_mixed: c3_mask=$(prob.c3_mask) needs " *
              "$(3 + mixed_n_aux(prob)) fields, got $nfields — build the FE spaces with " *
              "n_aux = mixed_n_aux(prob)")
    Γ  = BoundaryTriangulation(get_background_model(trian))
    dΓ = Measure(Γ, bdeg)
    nΓ = get_normal_vector(Γ)
    res(t, u, v) = global_residual_mixed(t, u, v, prob, trian, dΩh, dΓ, nΓ)
    use_ad && return TransientFEOperator(res, U, V)          # exact, and far too slow
    #  QUASI-NEWTON (the default): hand Jacobians on the physics rows — they index only
    #  u[1..3], so they are valid unchanged on the wider MultiField — plus the auxiliary
    #  mass block. `∂R/∂u̇` needs nothing extra: the auxiliary rows carry no time derivative,
    #  so its rows there are correctly zero.
    jac(t, u, du, v)    = jacobian_u(t, u, du, v, prob, trian, dΩh) +
                          mixed_coupling_jacobian(prob, u, du, v, trian, dΩh, dΓ, nΓ;
                                                  coupling = coupling)
    jac_t(t, u, dut, v) = jacobian_u_t(t, u, dut, v, prob, trian, dΩh)
    return TransientFEOperator(res, jac, jac_t, U, V)
end

"""
    make_initial_conditions_mixed(U, Nσ; eta0_func=nothing, ux0_func=nothing, uy0_func=nothing)

Zero (or prescribed) start for the 7-field layout. ⚠ From REST every field is zero and
the auxiliary constraints are satisfied exactly. From a non-rest start they are NOT —
𝖦 and 𝖥 would have to be solved for consistently first, and this helper does not do
that, so a non-rest mixed start is inconsistent at t₀ by construction.
"""
function make_initial_conditions_mixed(U, Nσ::Int; n_aux::Int = 2, eta0_func = nothing,
                                       ux0_func = nothing, uy0_func = nothing)
    zvv  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    eta0 = isnothing(eta0_func) ? (x -> 0.0) : eta0_func
    ux0  = isnothing(ux0_func)  ? (x -> zvv) : ux0_func
    uy0  = isnothing(uy0_func)  ? (x -> zvv) : uy0_func
    zf   = x -> zvv
    return interpolate_everywhere([eta0, ux0, uy0, fill(zf, n_aux)...], U)
end

# ---------------------------------------------------------------
#  Known weaknesses of this diagnostic — read before quoting a result
# ---------------------------------------------------------------
#  1. INF-SUP. 𝖦, 𝖥 are in the velocity space. Untested pairing; a negative result is
#     ambiguous (see the header).
#  2. ✅ THE BOUNDARY TERM IS NOW ASSEMBLED (was dropped in the first draft; G4 measured
#     a relative error of 22.4 in the boundary cells without it, against 2.7e-06 in the
#     interior). Re-check G4 after ANY change to the auxiliary rows.
#  3. COST. 7 fields against 3, with AD Jacobians — expect a large factor per step. 1-D
#     only.
#  4. THE LIVE-KNOB CONTROL IS MANDATORY. Compare every mixed `:full` run against
#     `:native` at the same settings. If they agree, the Class-III blocks are not
#     contributing and any stability is meaningless — that is exactly how the in-loop
#     attempt failed (NEW_TREATMENT.md §E.3), and rule 38d exists for it.
