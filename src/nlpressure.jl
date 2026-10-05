# ==============================================================
#  nlpressure.jl — the nonlinear pressure operator 𝓝 (𝓐 / 𝓚 / 𝓟 families)
#
#  All eight 𝓝_kj components in the three residual blocks of the derivation (bed-slope 𝓐,
#  surface-slope 𝓚, leading 𝓟), assembled whenever `prob.nl_pressure` (all eight, or none).
#
#  Treatments, by the regularity each term demands on C⁰ spaces:
#    c ∈ {3,6,7,8}           first-order everywhere (c=3's only 2nd derivative is the
#                            ANALYTIC bed Hessian) → direct, all three blocks
#                            (`nlp_direct_contrib`).
#    c ∈ {1,2,4,5}, 𝓐 half   EXACT integration by parts onto the test function
#                            (Ψ ∝ ∂_αh smooth; ∂Ψ gives first test derivatives + bed Hessian)
#                            → `nlp_gradh_contrib`.
#    c ∈ {1,2,4,5}, 𝓚, 𝓟     irreducible second derivatives of the unknowns (Class III) → the
#                            BROKEN formulation (src/broken.jl): the reduced contractions below
#                            (`_c3_sum`, `nlp_gradH_reduced_contrib`, `nlp_P_reduced_contrib`)
#                            fed with the cellwise Hessians, plus the skeleton layer.
#
#  Only `broken_class3_jacobian` linearises 𝓝 exactly; the {3,6,7,8} and 𝓐 blocks stay
#  quasi-Newton (O(A²), benign — rule 5). Slot bookkeeping
#  (a⊗b)[k,j]=a[k]b[j]: N¹,N² carry the differentiated divergence in the k slot (Ψ·U_a),
#  N³,N⁴,N⁵ carry the velocity in the k slot (U_a·Ψ).
# ==============================================================

"Contract the test layer-vector into the FIRST index of a constant 3-tensor:
(W ⋅ 𝓣)[k,j] = Σᵢ W[i]𝓣[i,k,j] → TensorValue{Nσ,Nσ}-CellField.
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
    nlp_direct_contrib(prob, d_cf, η, H, dhx, dhy, dHx, dHy, Ux, Uy, Wx, Wy, DW,
                       af, bf, S, DU, dΩh)

First-order nonlinear-pressure components c ∈ {3,6,7,8}, assembled directly in ALL THREE
blocks. `af = u·∇h`, `bf = u·∇H`, `S = ∇·(Hu)`, `DU = ∇·u` (stacked). The 1/H
of c ∈ {6,7,8} is cancelled analytically against one prefactor H (M_c ≡ H·N_c).
∂_a(𝖺) is expanded by hand so ∇ acts only on raw fields:
    ∂_a𝖺 = (∂_a∂_x h)𝖴x + (∂_xh)∂_a𝖴x + (∂_a∂_y h)𝖴y + (∂_yh)∂_a𝖴y.
Every term is subtracted (momentum RHS).
"""
function nlp_direct_contrib(prob::BALFEMProblem, d_cf, η, H, dhx, dhy, dHx, dHy,
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
    SD  = alg_outer(S, DU)                      # first-order remainder of N²
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

# ==============================================================
#  REDUCED Class-III contractions (NEW_TREATMENT.md §A.2): components {1,2,5} collapse EXACTLY
#  onto one contraction in ∇𝖲 (verified 4.4e-16). Fed by src/broken.jl with the cellwise parts
#  of the distributional gradients.
# ==============================================================

"""
    _c3_sum(prob, W, T2, T4, GU, SD, N4)

The Class-III contraction for ONE residual block — both irreducible objects the algebraic
reduction leaves (NEW_TREATMENT.md §A.2):

  * the ∇𝖲 family — `W ⊙ GU` (components {1,2,5} collapsed) **plus** `T² ⊙ SD`, the
    first-order remainder `s_k ∇·u_j` that 𝓝² contributes;
  * ∇𝖻 — `T⁴ ⊙ N4`, component 4.

Both are always assembled.
"""
_c3_sum(prob::BALFEMProblem, W, T2, T4, GU, SD, N4) =
    alg_dc3(W, GU) + alg_dc3(T2, SD) + alg_dc3(T4, N4)

"""
    nlp_gradH_reduced_contrib(prob, H, dHx, dHy, Wx, Wy, GU, SD, N4, dΩh)

Surface-slope (𝓚) half of c ∈ {1,2,4,5}, reduced: `W ⊙ GU + 𝓚³⁽²⁾ ⊙ SD + 𝓚³⁽⁴⁾ ⊙ N4`
instead of four separate contractions. EXACT — see NEW_TREATMENT.md §A.2.
"""
function nlp_gradH_reduced_contrib(prob::BALFEMProblem, H, dHx, dHy,
                                   Wx, Wy, GU, SD, N4, dΩh)
    NK = _c3_sum(prob, prob.WK3, prob.K3[2], prob.K3[4], GU, SD, N4)
    return ∫( (-1.0)*H*( dHx*(Wx ⋅ NK) + dHy*(Wy ⋅ NK) ) ) * dΩh
end

"""
    nlp_P_reduced_contrib(prob, H, DW, GU, SD, N4, dΩh)

Leading-pressure (𝓟) part of c ∈ {1,2,4,5}, reduced. EXACT — NEW_TREATMENT.md §A.2.
"""
function nlp_P_reduced_contrib(prob::BALFEMProblem, H, DW, GU, SD, N4, dΩh)
    NP = _c3_sum(prob, prob.WP3, prob.P3[2], prob.P3[4], GU, SD, N4)
    return ∫( (-1.0)*(H*H)*(NP ⋅ DW) ) * dΩh
end
