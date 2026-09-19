# ==============================================================
#  nlpressure.jl — FULL nonlinear pressure blocks (𝓐 / 𝓚 / 𝓟 families)
#
#  Completes the model physics: all eight 𝓝_kj components in the three
#  residual blocks of main.tex §8 (bed-slope 𝓐, surface-slope 𝓚, leading 𝓟).
#
#  Admissibility classes (main.tex ground rules) and treatments:
#    NATIVE  c ∈ {3,6,7,8}  first-order everywhere (c=3's only 2nd derivative
#            is the ANALYTIC bed Hessian) → direct, all blocks, all paths.
#    ∇h half c ∈ {1,2,4,5}  EXACT integration by parts onto the test function
#            (Ψ ∝ ∂_αh smooth; ∂Ψ gives first test derivatives + bed Hessian)
#            → nlp_gradh_contrib, residual-only, serial + distributed.
#    ∇H half + 𝓟 part, c ∈ {1,2,4,5}  irreducible (∂²η resp. ∂²(test)) →
#            FROZEN L²-PROJECTIONS: per step project 𝖲 and 𝖻 onto the velocity
#            FE space (π𝖲, π𝖻; SPD mass solves) and use ∂_a(π𝖲), ∂_a(π𝖻) lagged
#            one step. O(dt) on an O(A³) term. Serial (direct factorisation) and
#            distributed (CG + Jacobi mass solve) — see `build_nlp_ctx`.
#
#  All blocks are O(A²–A³) and treated quasi-Newton (they add to the residual
#  but not to the Jacobian — their contribution to convergence is negligible at
#  these amplitudes). Slot bookkeeping (a⊗b)[k,j]=a[k]b[j]:
#    N¹,N² carry the differentiated divergence in the k slot (Ψ·U_a),
#    N³,N⁴,N⁵ carry the velocity in the k slot (U_a·Ψ).
#  Verified at machine precision by test_nlpressure.jl gate G1.
# ==============================================================

"Contract the test layer-vector into the FIRST index of a constant 3-tensor:
(W ⋅ 𝓣)[k,j] = Σᵢ W[i]𝓣[i,k,j] → TensorValue{Nσ,Nσ}-CellField (native).
Linear in W ⇒ ∂_a(W ⋅ 𝓣) = (∂_aW) ⋅ 𝓣 — used to keep ∇ off composed test
expressions (∇ of an Operation-of-basis is not implemented for block arrays)."
alg_cont1(T::ThirdOrderTensorValue, W) = Operation(w -> w ⋅ T)(W)

"Bed-Hessian components ∂_a∂_b h of the ANALYTIC bathymetry (exact, via ∇∇)."
function alg_bed_hessian(d_cf)
    h2  = ∇∇(d_cf)
    hxx = Operation(T -> T[1,1])(h2)
    hxy = Operation(T -> T[1,2])(h2)
    hyy = Operation(T -> T[2,2])(h2)
    return hxx, hxy, hyy
end

"""
    nlp_native_contrib(prob, d_cf, η, H, dhx, dhy, dHx, dHy, Ux, Uy, Wx, Wy, DW,
                       af, bf, S, DU, dΩh)

First-order (native) nonlinear-pressure components c ∈ {3,6,7,8} in ALL THREE
blocks. `af = u·∇h`, `bf = u·∇H`, `S = ∇·(Hu)`, `DU = ∇·u` (stacked). The 1/H
of c ∈ {6,7,8} is cancelled analytically against one prefactor H (M_c ≡ H·N_c).
∂_a(𝖺) is expanded by hand so ∇ acts only on raw fields:
    ∂_a𝖺 = (∂_a∂_x h)𝖴x + (∂_xh)∂_a𝖴x + (∂_a∂_y h)𝖴y + (∂_yh)∂_a𝖴y.
Every term is subtracted (momentum RHS).
"""
function nlp_native_contrib(prob::BALFEMProblem, d_cf, η, H, dhx, dhy, dHx, dHy,
                            Ux, Uy, Wx, Wy, DW, af, bf, S, DU, dΩh)
    # M_c = H·N_c for c = 6,7,8 (first slot carries k = the u_k-type factor)
    M6 = (-1.0)*alg_outer(af, S)
    M7 = alg_outer(bf, S)
    M8 = (-1.0)*alg_outer(S, S)
    # N³ = −[U_a ⊗ ∂_a(𝖺)] with the hand-expanded ∂𝖺 (bed Hessian analytic).
    # On a flat bed both ∇h (dhx,dhy → 0, passed in) and the bed Hessian vanish, so N³ → 0.
    hxx, hxy, hyy = prob.flat_bed ? (0.0*d_cf, 0.0*d_cf, 0.0*d_cf) : alg_bed_hessian(d_cf)
    dxaf = hxx*Ux + dhx*alg_dx(Ux) + hxy*Uy + dhy*alg_dx(Uy)
    dyaf = hxy*Ux + dhx*alg_dy(Ux) + hyy*Uy + dhy*alg_dy(Uy)
    N3 = (-1.0)*(alg_outer(Ux, dxaf) + alg_outer(Uy, dyaf))

    # slope blocks (𝓐 with ∂_αh, 𝓚 with ∂_αH); H cancelled for 6–8, kept for 3
    NA68 = alg_dc3(prob.A3[6], M6) + alg_dc3(prob.A3[7], M7) + alg_dc3(prob.A3[8], M8)
    NK68 = alg_dc3(prob.K3[6], M6) + alg_dc3(prob.K3[7], M7) + alg_dc3(prob.K3[8], M8)
    NA3  = alg_dc3(prob.A3[3], N3)
    NK3  = alg_dc3(prob.K3[3], N3)
    r = ∫( (-1.0)*( dhx*(Wx ⋅ NA68) + dHx*(Wx ⋅ NK68)
                  + dhy*(Wy ⋅ NA68) + dHy*(Wy ⋅ NK68) ) ) * dΩh
    r = r + ∫( (-1.0)*H*( dhx*(Wx ⋅ NA3) + dHx*(Wx ⋅ NK3)
                        + dhy*(Wy ⋅ NA3) + dHy*(Wy ⋅ NK3) ) ) * dΩh

    # leading 𝓟 part: −∫H²(𝗣3[c]⊡N_c)·D_W ; with M_c = H·N_c → −∫H(𝗣3⊡M_c)·D_W
    NP68 = alg_dc3(prob.P3[6], M6) + alg_dc3(prob.P3[7], M7) + alg_dc3(prob.P3[8], M8)
    NP3  = alg_dc3(prob.P3[3], N3)
    r = r + ∫( (-1.0)*H*(NP68 ⋅ DW) ) * dΩh
    r = r + ∫( (-1.0)*(H*H)*(NP3 ⋅ DW) ) * dΩh
    return r
end

"""
    nlp_gradh_contrib(prob, d_cf, η, H, dhx, dhy, Ux, Uy, Wx, Wy, af, bf, S, DU, dΩh)

Bed-slope (𝓐, ∇h-prefactored) half of components c ∈ {1,2,4,5} via the EXACT
integration by parts onto the test function (main.tex §8):

    −∫Ψ₁⊙N¹ →  −∫ 𝖲⋅[∂x(Ψ₁⋅Ux)+∂y(Ψ₁⋅Uy)]
    −∫Ψ₂⊙N² →  +∫ 𝖲⋅[∂x(Ψ₂⋅Ux)+∂y(Ψ₂⋅Uy)] − ∫Ψ₂⊙(𝖲⊗DU)
    −∫Ψ₄⊙N⁴ →  +∫ 𝖻⋅[∂x(Ux⋅Ψ₄)+∂y(Uy⋅Ψ₄)]
    −∫Ψ₅⊙N⁵ →  −∫ 𝖲⋅[∂x(Ux⋅Ψ₅)+∂y(Uy⋅Ψ₅)]

with Ψ_c = (H∂_αh)(W_α ⋅ 𝓐3[c]), summed over α ∈ {x,y}. The IBP'd derivatives
are expanded BY HAND (∇ of an Operation-of-basis is not block-implemented),
using linearity of the contraction in the test and Σ_a Ψ·∂_aU_a = Ψ·DU:

    ∂x(Ψ⋅Ux)+∂y(Ψ⋅Uy) = (∂xΨ)⋅Ux + (∂yΨ)⋅Uy + Ψ⋅DU        (k-slot, c=1,2)
    ∂x(Ux⋅Ψ)+∂y(Uy⋅Ψ) = Ux⋅(∂xΨ) + Uy⋅(∂yΨ) + DU⋅Ψ        (j-slot, c=4,5)
    ∂_aΨ_c = [(∂_aH)∂_αh + H ∂_a∂_αh]·(W_α⋅𝓐3[c]) + (H∂_αh)·((∂_aW_α)⋅𝓐3[c])

so ∇ acts only on raw bases/fields; the surviving second derivative is the
ANALYTIC bed Hessian. Boundary integrals vanish on solid walls (u·n = 0
enters every q) / behind sponges — dropped. Vanishes identically on a flat
bed. Residual-only (no per-step state) → works serial AND distributed.
"""
function nlp_gradh_contrib(prob::BALFEMProblem, d_cf, η, H, dhx, dhy,
                           Ux, Uy, Wx, Wy, af, bf, S, DU, dΩh)
    hxx, hxy, hyy = alg_bed_hessian(d_cf)
    ddx = dhx + alg_dx(η)                       # ∂_x H
    ddy = dhy + alg_dy(η)                       # ∂_y H
    SD  = alg_outer(S, DU)                      # native remainder of N²
    r = nothing
    for (sα, hαx, hαy, Wα) in ((dhx, hxx, hxy, Wx), (dhy, hxy, hyy, Wy))
        Hs = H * sα
        gx = ddx*sα + H*hαx                     # ∂_x(H ∂_αh)
        gy = ddy*sα + H*hαy                     # ∂_y(H ∂_αh)
        dWαx = alg_dx(Wα); dWαy = alg_dy(Wα)

        # per component: C = W_α⋅𝓣, Cx/Cy = (∂W_α)⋅𝓣 (linearity in the test)
        C1  = alg_cont1(prob.A3[1], Wα)
        C1x = alg_cont1(prob.A3[1], dWαx); C1y = alg_cont1(prob.A3[1], dWαy)
        C2  = alg_cont1(prob.A3[2], Wα)
        C2x = alg_cont1(prob.A3[2], dWαx); C2y = alg_cont1(prob.A3[2], dWαy)
        C4  = alg_cont1(prob.A3[4], Wα)
        C4x = alg_cont1(prob.A3[4], dWαx); C4y = alg_cont1(prob.A3[4], dWαy)
        C5  = alg_cont1(prob.A3[5], Wα)
        C5x = alg_cont1(prob.A3[5], dWαx); C5y = alg_cont1(prob.A3[5], dWαy)

        # div_q, k-slot (Ψ·U):  Σ_a ∂_a(Ψ⋅U_a)
        dq1 = ((gx*C1 + Hs*C1x) ⋅ Ux) + ((gy*C1 + Hs*C1y) ⋅ Uy) + Hs*(C1 ⋅ DU)
        dq2 = ((gx*C2 + Hs*C2x) ⋅ Ux) + ((gy*C2 + Hs*C2y) ⋅ Uy) + Hs*(C2 ⋅ DU)
        # div_q, j-slot (U·Ψ):  Σ_a ∂_a(U_a⋅Ψ)
        dq4 = (Ux ⋅ (gx*C4 + Hs*C4x)) + (Uy ⋅ (gy*C4 + Hs*C4y)) + Hs*(DU ⋅ C4)
        dq5 = (Ux ⋅ (gx*C5 + Hs*C5x)) + (Uy ⋅ (gy*C5 + Hs*C5y)) + Hs*(DU ⋅ C5)

        t = ∫( (-1.0)*(S ⋅ dq1)
             + (S ⋅ dq2) - Hs*(C2 ⊙ SD)
             + (bf ⋅ dq4)
             + (-1.0)*(S ⋅ dq5) ) * dΩh
        r = r === nothing ? t : r + t
    end
    return r
end

"""
    nlp_frozen_N(Ux, Uy, S, DU, piS, pib) -> (N1, N2, N4, N5)

Components c ∈ {1,2,4,5} as TensorValue{Nσ,Nσ} fields with the irreducible
second derivatives replaced by derivatives of the FROZEN projections
`π𝖲, π𝖻` (previous-step FEFunctions on the velocity space):
    N¹ = −[(∂ₓπ𝖲)⊗Ux + (∂ᵧπ𝖲)⊗Uy]        N² = −N¹ + 𝖲⊗DU
    N⁴ = +[Ux⊗(∂ₓπ𝖻) + Uy⊗(∂ᵧπ𝖻)]        N⁵ = −[Ux⊗(∂ₓπ𝖲) + Uy⊗(∂ᵧπ𝖲)]
"""
function nlp_frozen_N(Ux, Uy, S, DU, piS, pib)
    dSx = alg_dx(piS); dSy = alg_dy(piS)
    dbx = alg_dx(pib); dby = alg_dy(pib)
    N1 = (-1.0)*(alg_outer(dSx, Ux) + alg_outer(dSy, Uy))
    N2 = (-1.0)*N1 + alg_outer(S, DU)
    N4 = alg_outer(Ux, dbx) + alg_outer(Uy, dby)
    N5 = (-1.0)*(alg_outer(Ux, dSx) + alg_outer(Uy, dSy))
    return N1, N2, N4, N5
end

"""
    nlp_gradH_frozen_contrib(prob, H, dHx, dHy, Wx, Wy, N1, N2, N4, N5, dΩh)

Surface-slope (𝓚, ∇H-prefactored) half of c ∈ {1,2,4,5} using the frozen
projections (irreducible ∂²η — IBP does not help here).
"""
function nlp_gradH_frozen_contrib(prob::BALFEMProblem, H, dHx, dHy,
                                  Wx, Wy, N1, N2, N4, N5, dΩh)
    NK = alg_dc3(prob.K3[1], N1) + alg_dc3(prob.K3[2], N2) +
         alg_dc3(prob.K3[4], N4) + alg_dc3(prob.K3[5], N5)
    return ∫( (-1.0)*H*( dHx*(Wx ⋅ NK) + dHy*(Wy ⋅ NK) ) ) * dΩh
end

"""
    nlp_P_frozen_contrib(prob, H, DW, N1, N2, N4, N5, dΩh)

Leading-pressure (𝓟) part of c ∈ {1,2,4,5} using the frozen projections
(IBP unusable — it would need second TEST derivatives through D_W).
"""
function nlp_P_frozen_contrib(prob::BALFEMProblem, H, DW, N1, N2, N4, N5, dΩh)
    NP = alg_dc3(prob.P3[1], N1) + alg_dc3(prob.P3[2], N2) +
         alg_dc3(prob.P3[4], N4) + alg_dc3(prob.P3[5], N5)
    return ∫( (-1.0)*(H*H)*(NP ⋅ DW) ) * dΩh
end

# ----------------------------------------------------------
#  Frozen-projection machinery (sequential AND distributed time loop)
# ----------------------------------------------------------

"""
    build_nlp_ctx(model, p_u, Nσ, trian, dΩh; distributed=false,
                  cg_rtol=1e-10, cg_maxiter=500)

Projection context for `nl_pressure_full`: an UNCONSTRAINED VectorValue{Nσ}
FE space (same reffe as the velocities) and its mass matrix, solved ONCE per
step for two right-hand sides (`π𝖲, π𝖻`). The mass matrix is SPD and
well-conditioned (unlike the advection-dominated Jacobian), so:
  * `distributed=false` (sequential): direct `lu` factorisation, ONE-OFF.
  * `distributed=true`: `CGSolver(JacobiLinearSolver())` from GridapSolvers —
    the same Jacobi-preconditioned Krylov family as the main distributed
    Newton solve (`build_ode_solver_distributed`), since base `lu` has no
    method for a partitioned `PSparseMatrix`. `numerical_setup` is built ONCE
    from the assembled mass matrix and reused every step (`solve!`).

RHS/solution vectors are allocated from the matrix itself (`allocate_in_range`/
`allocate_in_domain`) and the RHS is assembled IN-PLACE into that buffer
(`assemble_vector!`) — required distributed: an `assemble_vector(f, V)` result
and `allocate_in_domain(A)` are only isomorphic, not the SAME `PRange` object
(distinct ghost/assembly-cache layout), and `mul!`/CG inside `solve!` asserts
exact partition equality. `allocate_in_range`/`allocate_in_domain` are the
matrix-derived vector types that are *guaranteed* compatible with `A` on both
paths.
"""
function build_nlp_ctx(model, p_u::Int, Nσ::Int, trian, dΩh;
                       distributed::Bool = false,
                       cg_rtol::Float64 = 1e-10, cg_maxiter::Int = 500)
    reffe = ReferenceFE(lagrangian, VectorValue{Nσ,Float64}, p_u)
    Vp = FESpace(model, reffe; conformity=:H1)
    Up = TrialFESpace(Vp)
    a(u, v) = ∫( u ⋅ v ) * dΩh
    Mmass = assemble_matrix(a, Up, Vp)
    if distributed
        ls = CGSolver(JacobiLinearSolver(); rtol=cg_rtol, atol=1e-14, maxiter=cg_maxiter)
        ns = numerical_setup(symbolic_setup(ls, Mmass), Mmass)
        solvefun = (x, r) -> solve!(x, ns, r)
    else
        Mlu = lu(Mmass)
        solvefun = (x, r) -> ldiv!(x, Mlu, r)
    end
    return (Vp=Vp, Up=Up, dΩh=dΩh, trian=trian, Mmass=Mmass, solve=solvefun)
end

"""
    update_nlp_state!(prob, ctx, u_n)

After an accepted step: project `𝖲 = ∇·(Hu)` and `𝖻 = u·∇H` (from `u_n`) onto
the projection space and store `(π𝖲, π𝖻)` on `prob.nlp_state[]` for the next
step's residual (project-then-differentiate, one-step lag). Uses `ctx.solve`
(direct `lu` sequential / `CGSolver` distributed — set by `build_nlp_ctx`);
RHS/solution vectors are matrix-derived (see `build_nlp_ctx` docstring) so
they are exactly compatible with `ctx.Mmass`'s partition on both paths.
"""
function update_nlp_state!(prob::BALFEMProblem, ctx, u_n)
    η  = u_n[1];  Ux = u_n[2];  Uy = u_n[3]
    d_cf = CellField(prob.h_bathy, ctx.trian)
    H  = d_cf + η
    dHx = alg_dx(d_cf) + alg_dx(η);  dHy = alg_dy(d_cf) + alg_dy(η)
    DU = alg_dx(Ux) + alg_dy(Uy)
    b  = dHx*Ux + dHy*Uy
    S  = H*DU + b

    rS = allocate_in_range(ctx.Mmass); fill!(rS, zero(eltype(rS)))
    rb = allocate_in_range(ctx.Mmass); fill!(rb, zero(eltype(rb)))
    assemble_vector!(v -> ∫( S ⋅ v ) * ctx.dΩh, rS, ctx.Vp)
    assemble_vector!(v -> ∫( b ⋅ v ) * ctx.dΩh, rb, ctx.Vp)

    piS_vec = allocate_in_domain(ctx.Mmass); fill!(piS_vec, zero(eltype(piS_vec)))
    pib_vec = allocate_in_domain(ctx.Mmass); fill!(pib_vec, zero(eltype(pib_vec)))
    ctx.solve(piS_vec, rS)
    ctx.solve(pib_vec, rb)
    piS = FEFunction(ctx.Up, piS_vec)
    pib = FEFunction(ctx.Up, pib_vec)
    prob.nlp_state[] = (piS=piS, pib=pib)
    return nothing
end

# ==============================================================
#  REDUCED Class-III assembly + IN-LOOP (static-condensation) projections
#  Added on branch `new-classIII-treatment`.
#  Full derivation, verification and rationale: markdown_files/NEW_TREATMENT.md
# ==============================================================

"""
    nlp_class3_reduced_fields(Ux, Uy, S, DU, piS, pib) -> (GU, SD, N4)

The three TensorValue{Nσ,Nσ} objects the REDUCED Class-III assembly needs:

    GU = Σₐ ∂ₐ(π𝖲) ⊗ Uₐ     the single ∇𝖲-carrying object for components {1,2,5}
    SD = 𝖲 ⊗ DU              the ADMISSIBLE remainder of 𝓝² (first-order product)
    N4 = Σₐ Uₐ ⊗ ∂ₐ(π𝖻)      component 4, which needs ∇𝖻 and cannot join the others

Replaces the four-object `nlp_frozen_N` for the two PROJECTED blocks. `nlp_frozen_N`
is deliberately KEPT as the reference implementation that
`test_class3_residual_parity.jl` compares against — it is not dead code.
"""
function nlp_class3_reduced_fields(Ux, Uy, S, DU, piS, pib)
    dSx = alg_dx(piS); dSy = alg_dy(piS)
    dbx = alg_dx(pib); dby = alg_dy(pib)
    GU = alg_outer(dSx, Ux) + alg_outer(dSy, Uy)
    SD = alg_outer(S, DU)
    N4 = alg_outer(Ux, dbx) + alg_outer(Uy, dby)
    return GU, SD, N4
end

"""
    nlp_gradH_reduced_contrib(prob, H, dHx, dHy, Wx, Wy, GU, SD, N4, dΩh)

Surface-slope (𝓚) half of c ∈ {1,2,4,5}, reduced: `W ⊙ GU + 𝓚³⁽²⁾ ⊙ SD + 𝓚³⁽⁴⁾ ⊙ N4`
instead of four separate contractions. EXACT — see NEW_TREATMENT.md §A.2.
"""
function nlp_gradH_reduced_contrib(prob::BALFEMProblem, H, dHx, dHy,
                                   Wx, Wy, GU, SD, N4, dΩh)
    NK = alg_dc3(prob.WK3, GU) + alg_dc3(prob.K3[2], SD) + alg_dc3(prob.K3[4], N4)
    return ∫( (-1.0)*H*( dHx*(Wx ⋅ NK) + dHy*(Wy ⋅ NK) ) ) * dΩh
end

"""
    nlp_P_reduced_contrib(prob, H, DW, GU, SD, N4, dΩh)

Leading-pressure (𝓟) part of c ∈ {1,2,4,5}, reduced. EXACT — NEW_TREATMENT.md §A.2.
"""
function nlp_P_reduced_contrib(prob::BALFEMProblem, H, DW, GU, SD, N4, dΩh)
    NP = alg_dc3(prob.WP3, GU) + alg_dc3(prob.P3[2], SD) + alg_dc3(prob.P3[4], N4)
    return ∫( (-1.0)*(H*H)*(NP ⋅ DW) ) * dΩh
end

"""
    nlp_enable_inloop!(prob, ctx)

Switch `prob` into IN-LOOP (static-condensation) mode: the Class-III `L²` projections
are refreshed from the CURRENT Newton iterate inside `global_residual` instead of being
frozen from the previous accepted step. Pass `nothing` to return to lagged mode.
"""
function nlp_enable_inloop!(prob::BALFEMProblem, ctx)
    prob.nlp_ctx[] = ctx
    return prob
end

"""
    nlp_plain_iterate(u) -> Bool

⚠ THE AD GUARD, AND IT IS LOAD-BEARING.

`global_residual` is called in two very different situations:
  * **residual assembly** — `u` carries genuine `Float64` free values (the Newton iterate).
    This is the ONLY case in which a mass solve makes sense, and the only case that occurs
    at all on the default hand-Jacobian path (`TransientFEOperator(r,j,jt,U,V)` never
    differentiates `r`).
  * **AD Jacobian assembly** (`use_ad=true` ⇒ `TransientFEOperator(r,U,V)`) — Gridap
    differentiates `r`, so the cell data carries `ForwardDiff.Dual` numbers. Assembling a
    `Float64` RHS from those would throw, or silently truncate.

Returning `false` there makes the refresh a no-op, so the AD Jacobian is taken with `π`
held at the value the preceding residual evaluation computed — i.e. the SAME quasi-Newton
treatment these blocks already receive (rule 17b: the Jacobian sets the path, not the root).

Anything unexpected also returns `false`: failing safe here degrades to the legacy lagged
behaviour rather than throwing inside an assembly loop.
"""
function nlp_plain_iterate(u)
    try
        for k in 1:3
            v = get_free_dof_values(u[k])
            eltype(v) === Float64 || return false
        end
        return true
    catch
        return false
    end
end

"""
    refresh_nlp_state!(prob, ctx, S, b) -> Float64

Project `𝖲` and `𝖻` — as evaluated at the CURRENT Newton iterate — onto the projection
space and store them on `prob.nlp_state[]`. Returns the `∞`-norm change in `π𝖲`, which
`test_nlp_inloop.jl` G1 uses to assert the knob is live (rule 38d) and G2 to assert the
fixed point is reached.

This is `update_nlp_state!` without the one-step lag. Same mass matrix, same factorisation,
same cost per solve — only the state it is evaluated at differs. See NEW_TREATMENT.md §B.1
for why removing the lag costs no unknowns: the projection is static condensation of the
mixed formulation, and the freezing was only ever a decoupling convenience.
"""
"Diagnostic counter: how many times the in-loop refresh has actually fired.
Reset it before a run and read it after to check WHEN the refresh happens
(expect ~ Newton iterations x stages per step, not once per step)."
const NLP_REFRESH_COUNT = Ref(0)

function refresh_nlp_state!(prob::BALFEMProblem, ctx, S, b)
    NLP_REFRESH_COUNT[] += 1
    rS = allocate_in_range(ctx.Mmass); fill!(rS, zero(eltype(rS)))
    rb = allocate_in_range(ctx.Mmass); fill!(rb, zero(eltype(rb)))
    assemble_vector!(v -> ∫( S ⋅ v ) * ctx.dΩh, rS, ctx.Vp)
    assemble_vector!(v -> ∫( b ⋅ v ) * ctx.dΩh, rb, ctx.Vp)

    piS_vec = allocate_in_domain(ctx.Mmass); fill!(piS_vec, zero(eltype(piS_vec)))
    pib_vec = allocate_in_domain(ctx.Mmass); fill!(pib_vec, zero(eltype(pib_vec)))
    ctx.solve(piS_vec, rS)
    ctx.solve(pib_vec, rb)

    st  = prob.nlp_state[]
    chg = st === nothing ? Inf :
          maximum(abs, get_free_dof_values(st.piS) .- piS_vec)
    prob.nlp_state[] = (piS = FEFunction(ctx.Up, piS_vec),
                        pib = FEFunction(ctx.Up, pib_vec))
    return chg
end
