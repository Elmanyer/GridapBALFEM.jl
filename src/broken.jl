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
                          cip_hexp::Real = 2.0, cip_order::Int = 1,
                          stabilization::Symbol = :jumpgrad,
                          p_u::Int = 2, p_eta::Int = 2, degree::Int)
    stabilization in (:jumpgrad, :ghostvolume) ||
        error("attach_skeleton!: stabilization must be :jumpgrad or :ghostvolume (got :$stabilization)")
    stabilization === :ghostvolume && cip_order != 1 &&
        error("attach_skeleton!: cip_order applies to :jumpgrad only — :ghostvolume always covers " *
              "every order 0…p (GHOST_PENALTY_PLAN.md §1.4)")
    cip_gamma_u ≥ 0 && cip_gamma_eta ≥ 0 ||
        error("attach_skeleton!: the C⁰-IP coefficients must be ≥ 0 (a negative penalty " *
              "is an energy SOURCE)")
    1 ≤ cip_order ≤ CIP_MAX_ORDER ||
        error("attach_skeleton!: cip_order must be in 1:$CIP_MAX_ORDER (got $cip_order). Gridap " *
              "evaluates derivatives of the FE bases only up to order 2 " *
              "(BROKEN_FORMULATION_PLAN.md §5.3)")
    broken && !prob.nl_pressure_full &&
        error("attach_skeleton!: broken=true replaces the Class-III treatment and needs " *
              "nl_pressure=:full")
    if !broken && cip_gamma_u == 0 && cip_gamma_eta == 0
        prob.skel[] = nothing
        return prob
    end
    ctx = build_skeleton_ctx(model; degree = degree)
    g   = prob.g
    #  still-water scales — the penalty stays LINEAR in (η, 𝖴), so its Jacobian is exact
    τu  = CellField(x -> (d = prob.h_bathy(x); d * sqrt(g * d)), ctx.Λ)
    τη  = CellField(x -> sqrt(g * prob.h_bathy(x)), ctx.Λ)
    #  hp-CIP (Burman–Ern): order j ≤ min(cip_order, p_field), weight h_F^(s + 2(j−1)) — so the
    #  default s = 2 gives h^(2j), the scaling under which every order contributes O(1) damping
    #  on grid-scale modes and O((kh)^(2p+2)) on resolved ones. A jump of order j > p_field is
    #  identically zero (the j-th derivative of a degree-p polynomial vanishes for j > p), so the
    #  cap is a cost saving, not a change of operator.
    ou  = min(cip_order, p_u)
    oe  = min(cip_order, p_eta)
    hw  = j -> Operation(h -> h^(cip_hexp + 2*(j - 1)))(ctx.hF)
    gh  = stabilization === :ghostvolume && (cip_gamma_u > 0 || cip_gamma_eta > 0) ?
          build_ghost_ctx(model, ctx, τu, τη, Float64(cip_gamma_u), Float64(cip_gamma_eta),
                          Float64(cip_hexp), p_u, p_eta) : nothing
    prob.skel[] = (ctx..., broken = broken, stab = stabilization,
                   gu = Float64(cip_gamma_u), ge = Float64(cip_gamma_eta),
                   hexp = Float64(cip_hexp), order_u = ou, order_eta = oe,
                   cu = ntuple(j -> Float64(cip_gamma_u) * τu * hw(j), ou),
                   ce = ntuple(j -> Float64(cip_gamma_eta) * τη * hw(j), oe),
                   ghost = gh)
    return prob
end

"Highest normal-derivative order the C⁰-IP can penalise: Gridap's FE bases evaluate ∇ and ∇∇
only (Fields/FieldsInterfaces.jl:85–86; Polynomials/PolynomialInterfaces.jl)."
const CIP_MAX_ORDER = 2

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

"Jump of the SECOND normal derivative, ⟦∂ₙ²f⟧ = n_F·(∇∇f⁺ − ∇∇f⁻)·n_F, from the element-wise
Hessians (`alg_hess`); scalar or stacked field, or test basis."
function _jdn2(sk, f)
    nx, ny = sk.nx, sk.ny
    Hxx = alg_hess(f, 1, 1); Hxy = alg_hess(f, 1, 2); Hyy = alg_hess(f, 2, 2)
    return (nx*nx)*(Hxx.plus - Hxx.minus) + (2.0*nx*ny)*(Hxy.plus - Hxy.minus) +
           (ny*ny)*(Hyy.plus - Hyy.minus)
end

"⟦∂ₙʲf⟧ for j = 1, 2."
_jdnj(sk, f, j::Int) = j == 1 ? _jdn(sk, f) : _jdn2(sk, f)

"""
    cip_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)

The hp C⁰-IP penalty

    J_h = Σ_F Σ_{j ≤ order_u} ∫_F γ_u τ_u h_F^(s+2(j−1)) Σₐ ⟦∂ₙʲ𝖴ₐ⟧·Mv⟦∂ₙʲ𝖵ₐ⟧
        + Σ_F Σ_{j ≤ order_η} ∫_F γ_η τ_η h_F^(s+2(j−1)) ⟦∂ₙʲη⟧⟦∂ₙʲq⟧

LINEAR in (η, 𝖴): the same function is its own exact Jacobian (call it with the increments
`dη, dUx, dUy`). ⚠ Callers check `has_cip(prob)` first.
"""
function cip_contrib(prob::BALFEMProblem, sk, η, Ux, Uy, q, Wx, Wy)
    r = nothing
    if sk.gu > 0
        for j in 1:sk.order_u
            t = ∫( sk.cu[j] * ( (alg_mul(prob.Mv, _jdnj(sk, Ux, j)) ⋅ _jdnj(sk, Wx, j))
                              + (alg_mul(prob.Mv, _jdnj(sk, Uy, j)) ⋅ _jdnj(sk, Wy, j)) ) ) * sk.dΛ
            r = r === nothing ? t : r + t
        end
    end
    if sk.ge > 0
        for j in 1:sk.order_eta
            t = ∫( sk.ce[j] * (_jdnj(sk, η, j) * _jdnj(sk, q, j)) ) * sk.dΛ
            r = r === nothing ? t : r + t
        end
    end
    return r
end


# ==============================================================
#  :ghostvolume — the DIRECT (volume) ghost penalty (markdown_files/GHOST_PENALTY_PLAN.md)
#
#      J_F(u,v) = γ τ h^(s−3) ∫_{T⁺∪T⁻} (E u⁺ − E u⁻)·(E v⁺ − E v⁻)
#
#  E u^± = the polynomial of T^± EXTENDED over the patch (Preuß 2018; Lehrenfeld & Olshanskii
#  2019). Equivalent to the COMPLETE jump family Σ_{j,k ≤ p} G_jk ⟦∂ₙʲ·⟧⟦∂ₙᵏ·⟧ with the cross
#  terms (plan eq. ★) — every order up to p, WITHOUT any derivative (Gridap has none beyond 2).
#
#  EVALUATION (plan §1.3). The patch integral is a facet integral of a normal line integral,
#      ∫_{T⁺∪T⁻} g = Σ_{σ=±1} ∫_F ∫_0^h g(x_F + σ t e_a) dt     (e_a = the facet-normal axis),
#  and E u^± at x_F + σ t e_a is the cell's OWN polynomial composed with a constant REFERENCE
#  shift δ = σ t e_a / h_a (uniform Cartesian cells ⇒ the same δ in every cell). The skeleton
#  `.plus` / `.minus` of the shifted field then evaluate both polynomials at the SAME physical
#  point — the one inside its own cell, the other extended — and Gridap places the plus/minus
#  DOFs as for any skeleton term. Summing over BOTH signs σ covers T⁺ and T⁻ whatever the
#  plus/minus orientation, so no orientation bookkeeping is needed.
#  ⚠ Verified against the closed form (★) to round-off, including the x-periodic WRAP facet
#  (the shift is in each cell's reference coordinates) — test_broken_formulation.jl G9c.
# ==============================================================

"Constant 2-D reference-space translation ξ ↦ ξ + δ."
_shift_map(δ) = Gridap.Fields.AffineField(TensorValue(1.0, 0.0, 0.0, 1.0), Point(δ[1], δ[2]))

"Cell data of `f` composed with the translation δ (reference coordinates)."
function _shift_data(f, δ)
    d = get_data(f)
    return lazy_map(Gridap.Fields.Broadcasting(∘), d, fill(_shift_map(δ), length(d)))
end

"""
    _shift(f, δ)

`f` evaluated at reference points shifted by δ, i.e. each cell's polynomial EXTENDED by a
constant reference translation. Preserves the field's kind, so it can enter an assembled form:

  * multi-field basis component → shift its single-field basis and re-block it (Gridap's
    `similar_cell_field` is `@notimplemented` for `MultiFieldFEBasisComponent`);
  * single-field basis → `similar_cell_field` keeps the basis style;
  * any other CellField (FE function, transient component, AD-dual function) → `GenericCellField`.
"""
_shift(f::Gridap.MultiField.MultiFieldFEBasisComponent, δ) =
    Gridap.MultiField.MultiFieldFEBasisComponent(_shift(f.single_field, δ), f.fieldid, f.nfields)
_shift(f::Gridap.FESpaces.SingleFieldFEBasis, δ) =
    Gridap.CellData.similar_cell_field(f, _shift_data(f, δ), get_triangulation(f), DomainStyle(f))
function _shift(f::CellField, δ)
    DomainStyle(f) isa ReferenceDomain ||
        error("_shift: expected a ReferenceDomain CellField (the shift is in reference coordinates)")
    return GenericCellField(_shift_data(f, δ), get_triangulation(f), DomainStyle(f))
end

"""
    build_ghost_ctx(model, ctx, τu, τη, γu, γη, s, p_u, p_eta) -> NamedTuple

The geometry and quadrature of the direct ghost penalty:
  * uniform cell sizes (hx, hy) — ASSERTED, the reference shift is only a physical shift then;
  * the facet-normal axes that actually occur among the interior facets (an `ny = 1` flume has
    x-facets only), each with its indicator n_a² (exactly 1 on its facets, 0 elsewhere);
  * Gauss rules on (0, 1) for the normal line integral, (p+1) points per field — exact for the
    degree-2p integrand;
  * the coefficients γ τ h_a^(s−2) (= γ τ h^(s−3) × the physical length h_a of dt).
"""
function build_ghost_ctx(model, ctx, τu, τη, γu, γη, s, p_u, p_eta)
    trian = Triangulation(model)
    cc = get_cell_coordinates(trian)
    ext(c, k) = maximum(x -> x[k], c) - minimum(x -> x[k], c)
    hx = ext(cc[1], 1); hy = ext(cc[1], 2)
    all(c -> isapprox(ext(c, 1), hx; rtol = 1e-10) && isapprox(ext(c, 2), hy; rtol = 1e-10), cc) ||
        error("build_ghost_ctx: :ghostvolume needs a UNIFORM Cartesian mesh (the reference shift " *
              "is a physical shift only then)")
    occ(ind) = sum(∫(ind) * ctx.dΛ) > 1e-12
    axes = Any[]
    indx = ctx.nx * ctx.nx; indy = ctx.ny * ctx.ny
    occ(indx) && push!(axes, (1, hx, indx))
    occ(indy) && push!(axes, (2, hy, indy))
    rule(p) = (Q = Gridap.ReferenceFEs.Quadrature(SEGMENT, 2 * p);
               ([t[1] for t in Gridap.ReferenceFEs.get_coordinates(Q)],
                collect(Gridap.ReferenceFEs.get_weights(Q))))
    tu, wu = rule(p_u); te, we = rule(p_eta)
    cu = [γu * τu * ha^(s - 2) for (_, ha, _) in axes]
    ce = [γη * τη * ha^(s - 2) for (_, ha, _) in axes]
    return (hx = hx, hy = hy, axes = axes, tu = tu, wu = wu, te = te, we = we, cu = cu, ce = ce)
end

"""
    ghost_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)

The direct ghost penalty J_h (plan §1.1). LINEAR in (η, 𝖴): the same function is its own exact
Jacobian (call it with the increments). ⚠ Callers check `has_cip(prob)` first.
"""
function ghost_contrib(prob::BALFEMProblem, sk, η, Ux, Uy, q, Wx, Wy)
    gh = sk.ghost
    r = nothing
    for (ia, (ax, _, ind)) in enumerate(gh.axes), σ in (-1.0, 1.0)
        if sk.gu > 0
            for (t, w) in zip(gh.tu, gh.wu)
                δ = ax == 1 ? VectorValue(σ * t, 0.0) : VectorValue(0.0, σ * t)
                Uxs = _shift(Ux, δ); Uys = _shift(Uy, δ); Wxs = _shift(Wx, δ); Wys = _shift(Wy, δ)
                dUx = Uxs.plus - Uxs.minus; dUy = Uys.plus - Uys.minus
                dWx = Wxs.plus - Wxs.minus; dWy = Wys.plus - Wys.minus
                term = ∫( (w * gh.cu[ia]) * ind * ( (alg_mul(prob.Mv, dUx) ⋅ dWx)
                                                  + (alg_mul(prob.Mv, dUy) ⋅ dWy) ) ) * sk.dΛ
                r = r === nothing ? term : r + term
            end
        end
        if sk.ge > 0
            for (t, w) in zip(gh.te, gh.we)
                δ = ax == 1 ? VectorValue(σ * t, 0.0) : VectorValue(0.0, σ * t)
                ηs = _shift(η, δ); qs = _shift(q, δ)
                term = ∫( (w * gh.ce[ia]) * ind * ((ηs.plus - ηs.minus) * (qs.plus - qs.minus)) ) * sk.dΛ
                r = r === nothing ? term : r + term
            end
        end
    end
    return r
end

"The stabilisation term selected on `prob.skel[]`: :jumpgrad (hp C⁰-IP) or :ghostvolume."
stab_contrib(prob::BALFEMProblem, sk, η, Ux, Uy, q, Wx, Wy) =
    sk.stab === :ghostvolume ? ghost_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy) :
                               cip_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)


# ==============================================================
#  EXACT Jacobian of the broken Class-III blocks (volume + skeleton layer), 2026-10-02
#
#  The broken 𝓚/𝓟 Class-III residual (`global_residual` → `broken_class3_cell_fields` +
#  `nlp_gradH_reduced_contrib` + `nlp_P_reduced_contrib` + `broken_class3_skeleton_contrib`) is a
#  polynomial in (η, 𝖴) and their first and second derivatives, so its directional derivative
#  in (dη, d𝖴) follows by the product rule, term by term. Every reduced contraction (_c3_sum,
#  _c3_layer) is LINEAR in its tensor arguments, so δ(contraction) = contraction(δ arguments).
#  Verified against a central finite difference of the Class-III residual alone —
#  test_broken_formulation.jl G10.
# ==============================================================

"State fields shared by the broken Class-III residual and its Jacobian (cellwise)."
function _broken_state(prob::BALFEMProblem, d_cf, η, Ux, Uy)
    H   = d_cf + η
    dhx = prob.flat_bed ? 0.0*alg_dx(d_cf) : alg_dx(d_cf)
    dhy = prob.flat_bed ? 0.0*alg_dy(d_cf) : alg_dy(d_cf)
    dHx = dhx + alg_dx(η);  dHy = dhy + alg_dy(η)
    DU  = alg_dx(Ux) + alg_dy(Uy)
    bf  = dHx*Ux + dHy*Uy
    S   = H*DU + bf
    hxx, hxy, hyy = prob.flat_bed ? (0.0*d_cf, 0.0*d_cf, 0.0*d_cf) : alg_bed_hessian(d_cf)
    Hxx = hxx + alg_hess(η, 1, 1); Hxy = hxy + alg_hess(η, 1, 2); Hyy = hyy + alg_hess(η, 2, 2)
    dxUx = alg_dx(Ux); dyUx = alg_dy(Ux); dxUy = alg_dx(Uy); dyUy = alg_dy(Uy)
    dxD = alg_hess(Ux, 1, 1) + alg_hess(Uy, 1, 2)
    dyD = alg_hess(Ux, 2, 1) + alg_hess(Uy, 2, 2)
    dxb = Hxx*Ux + dHx*dxUx + Hxy*Uy + dHy*dxUy
    dyb = Hxy*Ux + dHx*dyUx + Hyy*Uy + dHy*dyUy
    dxS = dHx*DU + H*dxD + dxb
    dyS = dHy*DU + H*dyD + dyb
    return (H=H, dHx=dHx, dHy=dHy, DU=DU, bf=bf, S=S, Hxx=Hxx, Hxy=Hxy, Hyy=Hyy,
            dxUx=dxUx, dyUx=dyUx, dxUy=dxUy, dyUy=dyUy, dxD=dxD, dyD=dyD,
            dxb=dxb, dyb=dyb, dxS=dxS, dyS=dyS)
end

"Directional derivative of `_broken_state` in the direction (dη, d𝖴x, d𝖴y)."
function _broken_state_lin(F, dη, dUx, dUy, Ux, Uy)
    δH   = dη
    δdHx = alg_dx(dη);  δdHy = alg_dy(dη)
    δDU  = alg_dx(dUx) + alg_dy(dUy)
    δbf  = δdHx*Ux + δdHy*Uy + F.dHx*dUx + F.dHy*dUy
    δS   = δH*F.DU + F.H*δDU + δbf
    δHxx = alg_hess(dη, 1, 1); δHxy = alg_hess(dη, 1, 2); δHyy = alg_hess(dη, 2, 2)
    δdxUx = alg_dx(dUx); δdyUx = alg_dy(dUx); δdxUy = alg_dx(dUy); δdyUy = alg_dy(dUy)
    δdxD = alg_hess(dUx, 1, 1) + alg_hess(dUy, 1, 2)
    δdyD = alg_hess(dUx, 2, 1) + alg_hess(dUy, 2, 2)
    δdxb = δHxx*Ux + F.Hxx*dUx + δdHx*F.dxUx + F.dHx*δdxUx + δHxy*Uy + F.Hxy*dUy + δdHy*F.dxUy + F.dHy*δdxUy
    δdyb = δHxy*Ux + F.Hxy*dUx + δdHx*F.dyUx + F.dHx*δdyUx + δHyy*Uy + F.Hyy*dUy + δdHy*F.dyUy + F.dHy*δdyUy
    δdxS = δdHx*F.DU + F.dHx*δDU + δH*F.dxD + F.H*δdxD + δdxb
    δdyS = δdHy*F.DU + F.dHy*δDU + δH*F.dyD + F.H*δdyD + δdyb
    return (δH=δH, δdHx=δdHx, δdHy=δdHy, δDU=δDU, δbf=δbf, δS=δS, δdxb=δdxb, δdyb=δdyb,
            δdxS=δdxS, δdyS=δdyS)
end

"""
    broken_class3_jacobian(prob, sk, d_cf, η, Ux, Uy, dη, dUx, dUy, Wx, Wy, dΩh)

Exact ∂R/∂(η,𝖴)·(dη,d𝖴) of the broken Class-III 𝓚 and 𝓟 blocks, volume AND skeleton layer.
"""
function broken_class3_jacobian(prob::BALFEMProblem, sk, d_cf, η, Ux, Uy, dη, dUx, dUy, Wx, Wy, dΩh)
    F  = _broken_state(prob, d_cf, η, Ux, Uy)
    δ  = _broken_state_lin(F, dη, dUx, dUy, Ux, Uy)
    DW = alg_dx(Wx) + alg_dy(Wy)
    # ---- volume: R_K = −∫ H Σₐ ∂ₐH (Wₐ⋅NK),  R_P = −∫ H² (NP⋅D_W) -------------------------
    GU  = alg_outer(F.dxS, Ux) + alg_outer(F.dyS, Uy)
    SD  = alg_outer(F.S, F.DU)
    N4  = alg_outer(Ux, F.dxb) + alg_outer(Uy, F.dyb)
    δGU = alg_outer(δ.δdxS, Ux) + alg_outer(F.dxS, dUx) + alg_outer(δ.δdyS, Uy) + alg_outer(F.dyS, dUy)
    δSD = alg_outer(δ.δS, F.DU) + alg_outer(F.S, δ.δDU)
    δN4 = alg_outer(dUx, F.dxb) + alg_outer(Ux, δ.δdxb) + alg_outer(dUy, F.dyb) + alg_outer(Uy, δ.δdyb)
    NK  = _c3_sum(prob, prob.WK3, prob.K3[2], prob.K3[4], GU, SD, N4)
    δNK = _c3_sum(prob, prob.WK3, prob.K3[2], prob.K3[4], δGU, δSD, δN4)
    NP  = _c3_sum(prob, prob.WP3, prob.P3[2], prob.P3[4], GU, SD, N4)
    δNP = _c3_sum(prob, prob.WP3, prob.P3[2], prob.P3[4], δGU, δSD, δN4)
    r = ∫( (-1.0)*( δ.δH*( F.dHx*(Wx ⋅ NK) + F.dHy*(Wy ⋅ NK) )
                  + F.H*( δ.δdHx*(Wx ⋅ NK) + δ.δdHy*(Wy ⋅ NK) )
                  + F.H*( F.dHx*(Wx ⋅ δNK) + F.dHy*(Wy ⋅ δNK) ) )
           + (-1.0)*( (2.0*F.H*δ.δH)*(NP ⋅ DW) + (F.H*F.H)*(δNP ⋅ DW) ) ) * dΩh
    # ---- skeleton layer (broken_class3_skeleton_contrib, plan eq. "broken audit final") -----
    nx, ny = sk.nx, sk.ny
    Hf  = F.H.plus;      δHf = dη.plus
    Un  = nx*Ux.plus + ny*Uy.plus
    δUn = nx*dUx.plus + ny*dUy.plus
    jS  = F.S.plus - F.S.minus;     δjS = δ.δS.plus - δ.δS.minus
    jb  = F.bf.plus - F.bf.minus;   δjb = δ.δbf.plus - δ.δbf.minus
    LS  = alg_outer(jS, Un);  Lb = alg_outer(Un, jb)
    δLS = alg_outer(δjS, Un) + alg_outer(jS, δUn)
    δLb = alg_outer(δUn, jb) + alg_outer(Un, δjb)
    LK  = _c3_layer(prob, prob.WK3, prob.K3[4], LS, Lb);   δLK = _c3_layer(prob, prob.WK3, prob.K3[4], δLS, δLb)
    LP  = _c3_layer(prob, prob.WP3, prob.P3[4], LS, Lb);   δLP = _c3_layer(prob, prob.WP3, prob.P3[4], δLS, δLb)
    aHx = _avg(F.dHx); aHy = _avg(F.dHy); δaHx = _avg(δ.δdHx); δaHy = _avg(δ.δdHy)
    aWx = _avg(Wx); aWy = _avg(Wy); aDW = _avg(DW)
    r = r + ∫( δHf*( aHx*(aWx ⋅ LK) + aHy*(aWy ⋅ LK) )
             + Hf*( δaHx*(aWx ⋅ LK) + δaHy*(aWy ⋅ LK) )
             + Hf*( aHx*(aWx ⋅ δLK) + aHy*(aWy ⋅ δLK) )
             + (2.0*Hf*δHf)*(LP ⋅ aDW) + (Hf*Hf)*(δLP ⋅ aDW) ) * sk.dΛ
    return r
end
