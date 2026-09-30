# ==============================================================
#  broken.jl — the BROKEN (skeleton) Class-III formulation and the C⁰-IP penalty
#
#  Mathematics: LaTeX §"Broken Weak Formulation: Term-by-Term Audit"
#  (latex_docs/BALFEM_models/NumericalImplementation/BrokenAudit.tex).
#  Implementation record and campaign: markdown_files/BROKEN_FORMULATION_PLAN.md.
#
#  TWO ORTHOGONAL OPTIONS, both carried by `prob.skel[]` (nothing ⇒ bit-identical to main):
#
#  (1) `broken = true` — a third treatment of the Class-III 𝓚 (surface-slope) and 𝓟
#      (leading-pressure) blocks, beside the PROJECTED (frozen L² projections) and the MIXED
#      (auxiliary unknowns) ones. It inserts the DISTRIBUTIONAL gradient of the broken fields
#      𝖲 = ∇·(H𝗎) and 𝖻 = 𝗎·∇H directly:
#
#          ∇𝖲 = ∇_𝒯𝖲 − Σ_F ⟦𝖲⟧_F n_F δ_F           (⟦f⟧_F ≡ f⁺ − f⁻, n_F = n⁺)
#
#      — cellwise second derivatives of the C⁰ unknowns (via ∇∇, hand-expanded) PLUS one
#      skeleton integral for the layer. Form (A) of the audit: trial Hessians only, first
#      test derivatives only. No auxiliary unknowns, no projection, no lag. The 𝓐 (∇h) block
#      is NOT touched: its global IBP is already the exact distributional operator (audit
#      §"Bed-slope block").
#
#  (2) the C⁰ interior penalty  J_h = Σ_F ∫_F γ_u τ_u h_F^s Σₐ ⟦∂ₙ𝖴ₐ⟧·Mv⟦∂ₙ𝖵ₐ⟧
#                                         + γ_η τ_η h_F^s ⟦∂ₙη⟧⟦∂ₙq⟧,
#      τ_u = d√(gd), τ_η = √(gd), d = still-water depth. LINEAR in the unknowns (still-water
#      scales), so its Jacobian is exact and it adds nothing to ∂R/∂u̇ (a damping, not a mass
#      modification). Sign-definite, consistent (zero on the exact solution and on global
#      polynomials of degree ≤ p), mass-conserving (q ≡ 1 ⇒ 0). Works with ANY tier and ANY
#      Class-III treatment.
#
#  ⚠ JUMP CONVENTION. Every jump here is written EXPLICITLY as `f.plus − f.minus` with the
#  single normal n_F = nΛ.plus, rather than through Gridap's `jump`, so that the sign of the
#  layer is fixed by this file and verified by test_broken_formulation.jl G1 against the
#  (independently verified) mixed weak gradient, to round-off.
#
#  ⚠ EXCLUDED, DELIBERATELY: the skeleton term ⟦𝒫ᵢ⟧·vᵢ of the leading-pressure integration
#  by parts. It is consistent but turns the effective mass matrix into the element-wise strong
#  form (h⁻⁴ conditioning, non-convergent; audit Table 6.1). Do not add it.
#
#  SEQUENTIAL ONLY (like the mixed path). On an x-periodic mesh the identified edge is an
#  interior facet of the skeleton.
# ==============================================================

"(a,b) component of the ELEMENT-WISE Hessian of a scalar or stacked `VectorValue{Nσ}` field.
`∇∇(𝖴)` is a `ThirdOrderTensorValue{2,2,Nσ}` indexed [a,b,j] = ∂ₐ∂_b 𝖴ⱼ (symmetric in a,b);
`∇∇(η)` a `TensorValue{2,2}`. `∇` of an Operation-composed field is not implemented (rule 6),
so every cellwise second derivative is built from these raw-field Hessians."
alg_hess(u, a::Int, b::Int) = Operation(T -> _hess_ab(T, a, b))(∇∇(u))
_hess_ab(T::ThirdOrderTensorValue{2,2,N}, a, b) where {N} =
    VectorValue{N}(ntuple(j -> T[a, b, j], Val(N)))
_hess_ab(T::TensorValue{2,2}, a, b) = T[a, b]

"Mean of a field over the two sides of a facet (test factors, discontinuous prefactors)."
_avg(f) = 0.5 * (f.plus + f.minus)

"""
    build_skeleton_ctx(model; degree) -> NamedTuple

Interior-facet (skeleton) objects, built ONCE: `Λ`, `dΛ`, the unique facet normal
`n_F = nΛ.plus` split into scalar CellFields `nx, ny`, and the facet measure `hF`.
"""
function build_skeleton_ctx(model; degree::Int)
    Λ  = SkeletonTriangulation(model)
    dΛ = Measure(Λ, degree)
    nΛ = get_normal_vector(Λ)
    n⁺ = nΛ.plus
    nx = Operation(n -> n ⋅ Ex)(n⁺)
    ny = Operation(n -> n ⋅ Ey)(n⁺)
    hF = CellField(collect(get_cell_measure(Λ)), Λ)
    return (Λ = Λ, dΛ = dΛ, nx = nx, ny = ny, hF = hF, nfacets = num_cells(Λ))
end

"""
    attach_skeleton!(prob, model; broken=false, cip_gamma_u=0.0, cip_gamma_eta=0.0,
                     cip_hexp=2.0, degree) -> prob

Select the broken Class-III treatment and/or the C⁰-IP penalty by attaching the skeleton
context to `prob.skel[]`. With `broken=false` and both γ = 0 the field is left `nothing`
and every residual path is bit-identical to the Galerkin one (test_broken_formulation G3).
"""
function attach_skeleton!(prob::BALFEMProblem, model; broken::Bool = false,
                          cip_gamma_u::Real = 0.0, cip_gamma_eta::Real = 0.0,
                          cip_hexp::Real = 2.0, degree::Int)
    cip_gamma_u ≥ 0 && cip_gamma_eta ≥ 0 ||
        error("attach_skeleton!: the C⁰-IP coefficients must be ≥ 0 (a negative penalty " *
              "is an energy SOURCE)")
    broken && !prob.nl_pressure_full &&
        error("attach_skeleton!: broken=true replaces the Class-III treatment and needs " *
              "nl_pressure=:full")
    if !broken && cip_gamma_u == 0 && cip_gamma_eta == 0
        prob.skel[] = nothing
        return prob
    end
    ctx = build_skeleton_ctx(model; degree = degree)
    g   = prob.g
    hs  = Operation(h -> h^cip_hexp)(ctx.hF)
    #  still-water scales — the penalty stays LINEAR in (η, 𝖴), so its Jacobian is exact
    τu  = CellField(x -> (d = prob.h_bathy(x); d * sqrt(g * d)), ctx.Λ)
    τη  = CellField(x -> sqrt(g * prob.h_bathy(x)), ctx.Λ)
    prob.skel[] = (ctx..., broken = broken,
                   gu = Float64(cip_gamma_u), ge = Float64(cip_gamma_eta),
                   hexp = Float64(cip_hexp),
                   cu = Float64(cip_gamma_u) * τu * hs,
                   ce = Float64(cip_gamma_eta) * τη * hs)
    return prob
end

"True when `prob` assembles the Class-III 𝓚/𝓟 blocks by the broken (distributional) path."
is_broken(prob::BALFEMProblem) = (sk = prob.skel[]; sk !== nothing && sk.broken)

"True when `prob` carries a C⁰-IP penalty."
has_cip(prob::BALFEMProblem) = (sk = prob.skel[]; sk !== nothing && (sk.gu > 0 || sk.ge > 0))

"""
    broken_class3_cell_fields(prob, d_cf, η, H, dHx, dHy, Ux, Uy, DU, S)
        -> (GU_𝒯, SD, N4_𝒯, (dxS, dyS, dxb, dyb))

The CELLWISE part of the Class-III objects, from the element-wise Hessians (audit eq.
"cellwise grad s b"):

    ∂ₐ𝖲 = ∂ₐH·D + H·∂ₐD + ∂ₐ𝖻,              ∂ₐD = ∂ₐ∂ₓ𝖴x + ∂ₐ∂ᵧ𝖴y
    ∂ₐ𝖻 = (∂ₐ∂ₓH)𝖴x + ∂ₓH·∂ₐ𝖴x + (∂ₐ∂ᵧH)𝖴y + ∂ᵧH·∂ₐ𝖴y,    ∂ₐ∂_bH = ∂ₐ∂_b h + ∂ₐ∂_b η

    GU_𝒯 = Σₐ ∂ₐ𝖲 ⊗ 𝖴ₐ,   SD = 𝖲 ⊗ D,   N4_𝒯 = Σₐ 𝖴ₐ ⊗ ∂ₐ𝖻

— the exact analogues of `nlp_class3_reduced_fields`, with ∂ₐ(π𝖲) → ∂ₐ^𝒯𝖲. The skeleton
layer that completes them is `broken_class3_skeleton_contrib`.
"""
function broken_class3_cell_fields(prob::BALFEMProblem, d_cf, η, H, dHx, dHy, Ux, Uy, DU, S)
    hxx, hxy, hyy = prob.flat_bed ? (0.0*d_cf, 0.0*d_cf, 0.0*d_cf) : alg_bed_hessian(d_cf)
    Hxx = hxx + alg_hess(η, 1, 1)
    Hxy = hxy + alg_hess(η, 1, 2)
    Hyy = hyy + alg_hess(η, 2, 2)
    dxUx = alg_dx(Ux); dyUx = alg_dy(Ux)
    dxUy = alg_dx(Uy); dyUy = alg_dy(Uy)
    dxD = alg_hess(Ux, 1, 1) + alg_hess(Uy, 1, 2)
    dyD = alg_hess(Ux, 2, 1) + alg_hess(Uy, 2, 2)
    dxb = Hxx*Ux + dHx*dxUx + Hxy*Uy + dHy*dxUy
    dyb = Hxy*Ux + dHx*dyUx + Hyy*Uy + dHy*dyUy
    dxS = dHx*DU + H*dxD + dxb
    dyS = dHy*DU + H*dyD + dyb
    GU = alg_outer(dxS, Ux) + alg_outer(dyS, Uy)
    SD = alg_outer(S, DU)
    N4 = alg_outer(Ux, dxb) + alg_outer(Uy, dyb)
    return GU, SD, N4, (dxS, dyS, dxb, dyb)
end

"Skeleton-layer counterpart of `_c3_sum`: the same weights and the same `c3_mask` bits, acting
on the layer objects. `SD` has no layer (it is undifferentiated), so it is absent here."
function _c3_layer(prob::BALFEMProblem, W, T4, LS, Lb)
    use_gs, use_gb = prob.c3_mask
    if use_gs && use_gb
        return alg_dc3(W, LS) + alg_dc3(T4, Lb)
    elseif use_gs
        return alg_dc3(W, LS)
    else
        return alg_dc3(T4, Lb)
    end
end

"""
    broken_class3_skeleton_contrib(prob, sk, H, dHx, dHy, Ux, Uy, S, bf, Wx, Wy, DW)

The skeleton layer of the Class-III 𝓚 and 𝓟 blocks (audit eq. "broken audit final"):

    + Σ_F ∫_F H Σₐ{∂ₐH}{𝖵ₐ}·( W_K ⊙ (⟦𝖲⟧⊗𝖴ₙ) + 𝓚⁽⁴⁾ ⊙ (𝖴ₙ⊗⟦𝖻⟧) )
    + Σ_F ∫_F H² {D_W}·( W_P ⊙ (⟦𝖲⟧⊗𝖴ₙ) + 𝓟⁽⁴⁾ ⊙ (𝖴ₙ⊗⟦𝖻⟧) )

with 𝖴ₙ = n_x𝖴x + n_y𝖴y (continuous), ⟦f⟧ = f⁺ − f⁻ along n_F = n⁺. The PLUS sign: the layer
enters ∇f = ∇_𝒯f − Σ⟦f⟧n_Fδ_F with a minus, and the Class-III blocks are SUBTRACTED in the
residual. Verified (sign, normal, slot order) by test_broken_formulation.jl G1.
"""
function broken_class3_skeleton_contrib(prob::BALFEMProblem, sk, H, dHx, dHy,
                                        Ux, Uy, S, bf, Wx, Wy, DW)
    nx, ny = sk.nx, sk.ny
    Hf = H.plus                                     # H is continuous
    Un = nx*Ux.plus + ny*Uy.plus                    # 𝖴ₙ, continuous
    jS = S.plus - S.minus                           # ⟦𝖲⟧_F
    jb = bf.plus - bf.minus                         # ⟦𝖻⟧_F
    LS = alg_outer(jS, Un)
    Lb = alg_outer(Un, jb)
    NK = _c3_layer(prob, prob.WK3, prob.K3[4], LS, Lb)
    NP = _c3_layer(prob, prob.WP3, prob.P3[4], LS, Lb)
    return ∫( Hf*( _avg(dHx)*(_avg(Wx) ⋅ NK) + _avg(dHy)*(_avg(Wy) ⋅ NK) )
            + (Hf*Hf)*(NP ⋅ _avg(DW)) ) * sk.dΛ
end

"Jump of the normal derivative along n_F, ⟦∂ₙf⟧ = n_F·(∇f⁺ − ∇f⁻), for a scalar or stacked
field (or test basis)."
_jdn(sk, f) = sk.nx*(alg_dx(f).plus - alg_dx(f).minus) + sk.ny*(alg_dy(f).plus - alg_dy(f).minus)

"""
    cip_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)

The C⁰-IP penalty J_h. LINEAR in (η, 𝖴): the same function is its own exact Jacobian
(call it with the increments `dη, dUx, dUy`). ⚠ Callers check `has_cip(prob)` first.
"""
function cip_contrib(prob::BALFEMProblem, sk, η, Ux, Uy, q, Wx, Wy)
    r = nothing
    if sk.gu > 0
        t = ∫( sk.cu * ( (alg_mul(prob.Mv, _jdn(sk, Ux)) ⋅ _jdn(sk, Wx))
                       + (alg_mul(prob.Mv, _jdn(sk, Uy)) ⋅ _jdn(sk, Wy)) ) ) * sk.dΛ
        r = t
    end
    if sk.ge > 0
        t = ∫( sk.ce * (_jdn(sk, η) * _jdn(sk, q)) ) * sk.dΛ
        r = r === nothing ? t : r + t
    end
    return r
end
