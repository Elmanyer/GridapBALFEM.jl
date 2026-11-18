# ==============================================================
#  vms.jl — RESIDUAL-BASED (variational multiscale) STABILISATION  (branch mixed-vms-solver)
#
#  Plan and record: markdown_files/MIXED_VMS_PLAN.md §3. The only stabiliser of the mixed solver;
#  also available on models 3–4 (3 fields, already conforming).
#
#  STRONG RESIDUALS. R_c (scalar) and R_m,a (VectorValue{Nσ}, a = x, y): the strong form of the
#  equations `global_residual` discretises, obtained by undoing each of its integrations by parts
#  (`vms_strong_residuals`). ⚠ They are checked against the weak form, sign by sign, by
#  test_mixed_vms.jl G3: ∫ v·R_strong equals the assembled weak residual up to the inter-element
#  jumps of a C⁰ field, which vanish under refinement.
#
#  SUBSCALES. Quasi-static, algebraic:  η′ = −τ_η R_c ,  𝗎′ₐ = −T (H M)⁻¹ R_m,a  with
#
#      T M⁻¹ = τ_w e eᵀ + τ_a (M⁻¹ − e eᵀ)          e = (1,…,1): the depth-uniform velocity
#
#  P = e Φᵀ is the M-orthogonal projector onto e (M e = Φ, Φ·e = 1), so the first part acts on the
#  depth-INTEGRATED momentum residual (e·R = Σᵢ Rᵢ) with the wave time scale, the second on the
#  VERTICAL-STRUCTURE residual with the advective one. Both matrices are symmetric PSD.
#
#      τ_w = C_w h_e/(2c₀)     τ_η = C_η h_e/(2c₀)     c₀ = √(g d),  h_e = h/p_u
#      τ_a = C_a [ k_Δt + (2/h_e)² (u₀² + 𝗎·M𝗎) ]^(−1/2)      u₀ = 10⁻³ c₀,  k_Δt = (2/Δt)² or 0
#
#  TERMS (derivative on the test function; all ADDED to the residual):
#
#      ∫ Σₐ (N_b ∂_b Wₐ)·Zₐ ,  Zₐ = τ_w e (e·R_m,a) + τ_a (M⁻¹ − e eᵀ) R_m,a      SUPG, momentum
#      ∫ τ_w Σₐ ∂ₐq (e·R_m,a)                                                      PSPG, continuity
#      ∫ τ_η (ū·∇q) R_c                                                            SUPG, continuity
#      ∫ τ_η g H (Φ·∇·W) R_c                                                       LSIC, momentum
#
#  N_b[i,j] = Σ_k 𝓜[i,k,j] u_{b,k}: the layer-coupled advection matrix (symmetric, 𝓜 being fully
#  symmetric). NOT included, by choice: the adjoint of the vertical-advection tensor 𝓖 (not
#  symmetric), zero-order adjoint terms, cross/Reynolds subscale stresses, ∂ₜ of the subscales.
#
#  RESIDUAL AND JACOBIANS FROM ONE FORMULA (src/jets.jl): `vms_contrib(…)` is the residual,
#  `vms_contrib(…; du)` the exact ∂/∂u (τ_a's dependence on 𝗎 included), `vms_contrib(…; dut)` the
#  exact ∂/∂u̇.
#
#  LIMITS, STATED. R_P needs cellwise second derivatives of u̇ and η (use Q3/Q2: at Q2/Q1 the η
#  Hessian vanishes in 1-D). With 𝓝, ∇𝖲 → 𝖦 and ∇𝖻 → 𝖥 everywhere, so the only second derivative
#  of 𝗎 left is ∂ₐ(∇·𝗎) in ∂ₐ(H² NP). On a variable bed ∂ₐ(H² NP₃) needs the THIRD derivative of the
#  bathymetry; if it cannot be evaluated the stabilisation is refused, never truncated.
#
#  ⚠ MEASURED DEFECT OF THE ASGS FORM ABOVE (2026-11-17, MIXED_VMS_PLAN.md §9.3): R_m,a carries the
#  time-derivative part (acceleration + the dispersive −d²B̃∇∇·𝗎̇), so momentum SUPG adds an O(C_a),
#  non-symmetric, second-derivative block to ∂R/∂u̇ and the effective mass turns INDEFINITE once the
#  state moves. Exact-Newton CN runs then fail. The OSS option below removes it by construction.
#
#  ORTHOGONAL SUBGRID SCALES (`oss = true`; Codina 2000, 2002). The subscale is the part of the
#  residual ORTHOGONAL to the finite-element space:  𝗎′ₐ = −T (HM)⁻¹ P⊥(Rˢ_m,a),  η′ = −τ_η P⊥(Rˢ_c),
#  with P⊥ = I − P_h (L² projection) and Rˢ the STATIC residual: R with every u̇ term removed
#  (acceleration, the whole 𝓛(u̇) pressure package, η_t). Removing them is the OSS argument, not a
#  truncation: P⊥ annihilates what lies in the FE space (P⊥(∂ₜu_h) = 0 for a constant mass), and what
#  is left of the time part (the H and dispersive weights) is O(h^p) after projection, the order to
#  which OSS is consistent anyway. The four terms keep their shape with R → Rˢ − P_h(Rˢ).
#
#  The projections are UNKNOWNS (monolithic, exact Newton; the pattern of the mixed 𝖦, 𝖥):
#
#      Πₐ = P_h(Rˢ_m,a) ∈ VectorValue{Nσ}, velocity order     δ = P_h(Rˢ_c) ∈ scalar, surface order
#      rows   ∫ (Πₐ − Rˢ_m,a)·χₐ = 0 ,   ∫ (δ − Rˢ_c) ϑ = 0       (no time derivative: index-1 DAE)
#
#  appended after every other field (index `oss_i0` = 4 + mixed_n_aux). CONSEQUENCES: ∂R/∂u̇ is the
#  Galerkin mass exactly (the stabilisation carries no u̇); a hot start needs consistent projections
#  (`oss_consistent_ic`); the Jacobian of the terms and rows is exact (jets).
# ==============================================================

"True when `prob` carries a residual-based stabilisation record."
has_vms(prob::BALFEMProblem) = prob.vms[] !== nothing

"True when the attached stabilisation uses orthogonal subgrid scales (projection fields present)."
has_oss(prob::BALFEMProblem) = has_vms(prob) && prob.vms[].oss

"Number of OSS projection fields (Πx, Πy, δ) in the layout, 0 without OSS."
const OSS_N_FIELDS = 3
oss_n_fields(prob::BALFEMProblem) = has_oss(prob) ? OSS_N_FIELDS : 0

"""
    extra_field_zeros(prob) -> Vector{Any}

Zero initial functions for every field after [η, 𝖴x, 𝖴y], in layout order: the mixed 𝖦x, 𝖦y, 𝖥x,
𝖥y, then the OSS Πx, Πy (VectorValue{Nσ}) and δ (scalar). For `interpolate_everywhere`.
"""
function extra_field_zeros(prob::BALFEMProblem)
    zvv = VectorValue(ntuple(_ -> 0.0, prob.Nσ)...)
    zf  = x -> zvv
    out = Any[fill(zf, mixed_n_aux(prob))...]
    has_oss(prob) && append!(out, Any[zf, zf, x -> 0.0])
    return out
end

"Cell size √|K| as a cellwise-constant CellField; per rank on a distributed triangulation."
_cell_size(trian::Triangulation) = CellField(sqrt.(collect(get_cell_measure(trian))), trian)
_cell_size(trian::GridapDistributed.DistributedTriangulation) =
    CellField(map(t -> sqrt.(collect(get_cell_measure(t))), GridapDistributed.local_views(trian)), trian)

"ThirdOrderTensorValue → plain array (index-preserving)."
_t3_array(T::ThirdOrderTensorValue{N,N,N}) where {N} = [T[i, k, j] for i in 1:N, k in 1:N, j in 1:N]

"""
    _bed_third(prob, trian) -> (hxxx, hxxy, hxyy, hyyy) as CellFields, or `nothing`

Third derivatives of the analytic bathymetry, by differentiating the components of its Hessian.
`nothing` when the bathymetry function cannot be differentiated three times.
"""
function _bed_third(prob::BALFEMProblem, trian)
    h    = prob.h_bathy
    hess = ∇∇(h)
    comp(a, b) = x -> hess(x)[a, b]
    g11 = ∇(comp(1, 1)); g12 = ∇(comp(1, 2)); g22 = ∇(comp(2, 2))
    try
        x0 = Point(0.123, 0.456)
        v  = (g11(x0)[1], g11(x0)[2], g12(x0)[2], g22(x0)[2])
        all(isfinite, v) || return nothing
    catch
        return nothing
    end
    return (CellField(x -> g11(x)[1], trian), CellField(x -> g11(x)[2], trian),
            CellField(x -> g12(x)[2], trian), CellField(x -> g22(x)[2], trian))
end

"""
    vms_ctx(prob, trian; p_u, c_adv=1.0, c_wave=1.0, c_eta=1.0, dt=nothing, oss=false) -> NamedTuple

The record the stabilisation reads: the constant vertical matrices `e eᵀ`, `M⁻¹ − e eᵀ`, the
bed-slope Class-III weight `W_A`, the element size `h_e = √|K|/p_u`, the still-water speed and the
coefficients. `dt::Real` switches the `(2/Δt)²` term of `τ_a` on (off by default: with it the
semi-discrete operator depends on Δt). `oss = true`: orthogonal subgrid scales (file header); the
FE spaces must carry the three projection fields after every other field
(`build_fe_spaces(…; n_oss = 3)`).
"""
function vms_ctx(prob::BALFEMProblem, trian; p_u::Int, c_adv::Real = 1.0, c_wave::Real = 1.0,
                 c_eta::Real = 1.0, dt = nothing, oss::Bool = false)
    prob.linearised &&
        error("vms: the stabilisation is defined for the nonlinear regime only (models 3–6). The " *
              "linear models are stable and carry no advection to stabilise.")
    prob.nl_pressure && prob.formulation === :primal &&
        error("vms: nl_pressure=true needs formulation=:mixed. The strong residual of the primal 𝓝 " *
              "form carries third derivatives of the unknowns (markdown_files/MIXED_VMS_PLAN.md §1).")
    (prob.P_full && prob.lin_pressure && prob.advection) ||
        error("vms: needs the complete nonlinear core (P_full, lin_pressure, advection): build the " *
              "problem with build_problem, not a partial build_problem_raw flag set")
    (c_adv ≥ 0 && c_wave ≥ 0 && c_eta ≥ 0) ||
        error("vms: the coefficients must be ≥ 0 (a negative τ is an energy source)")
    N    = prob.Nσ
    Mm   = [prob.Mv[i, j] for i in 1:N, j in 1:N]
    E    = ones(N, N)
    Ebar = alg_to_tensor2(E)
    Qv   = alg_to_tensor2(inv(Mm) - E)
    evec = alg_to_vec(ones(N))
    WA3  = alg_to_tensor3(alg_class3_weight(_t3_array(prob.A3[1]), _t3_array(prob.A3[2]),
                                            _t3_array(prob.A3[5])))
    third = nothing
    if prob.nl_pressure && !prob.flat_bed
        third = _bed_third(prob, trian)
        third === nothing &&
            error("vms: model 6 (𝓝 over a variable bed) needs the third derivative of the bathymetry " *
                  "for ∂ₐ(H² NP₃), and this bathymetry function cannot be differentiated three times. " *
                  "Refused: the term is never dropped (MIXED_VMS_PLAN.md §3.4).")
    end
    g    = prob.g
    he   = (1.0 / p_u) * _cell_size(trian)
    c0   = CellField(x -> sqrt(g * prob.h_bathy(x)), trian)
    kdt  = dt === nothing ? 0.0 : (2.0 / Float64(dt))^2
    return (Ebar = Ebar, Qv = Qv, evec = evec, WA3 = WA3, third = third, he = he, c0 = c0,
            c_adv = Float64(c_adv), c_wave = Float64(c_wave), c_eta = Float64(c_eta),
            kdt = kdt, dt = dt, p_u = p_u, oss = oss, oss_i0 = 4 + mixed_n_aux(prob))
end

"""
    attach_vms!(prob, trian; p_u, c_adv=1.0, c_wave=1.0, c_eta=1.0, dt=nothing, oss=false) -> prob

Attach the residual-based stabilisation. From then on `global_residual`, `jacobian_u` and
`jacobian_u_t` add its terms. All three coefficients zero ⇒ the record is removed (refused with
`oss = true`: the projection fields would be left without rows, a singular system).
"""
function attach_vms!(prob::BALFEMProblem, trian; p_u::Int, c_adv::Real = 1.0, c_wave::Real = 1.0,
                     c_eta::Real = 1.0, dt = nothing, oss::Bool = false)
    if c_adv == 0 && c_wave == 0 && c_eta == 0
        oss && error("attach_vms!: oss=true with all coefficients 0 would leave the projection " *
                     "fields without rows; build the spaces without them instead")
        prob.vms[] = nothing
        return prob
    end
    prob.vms[] = vms_ctx(prob, trian; p_u, c_adv, c_wave, c_eta, dt, oss)
    return prob
end

"""
    vms_fields(prob, trian, u; du=nothing, dut=nothing, static=false) -> NamedTuple of jets

Everything the strong residuals read: the core fields, the time derivatives, the cellwise
Hessians of η, 𝗎 and u̇, the bed Hessian and, with 𝓝, the auxiliary fields and their gradients.
`du` seeds the tangents of the unknowns, `dut` those of u̇. `static = true` sets every u̇ jet to
the exact zero, so the strong residuals become the STATIC residuals Rˢ of OSS (every u̇ term of
R is linear and homogeneous in u̇, so it is pruned, not evaluated at zero).
"""
function vms_fields(prob::BALFEMProblem, trian, u; du = nothing, dut = nothing, static::Bool = false)
    C  = jet_core_fields(prob, trian, u; du = du)
    p(i)  = _part(du, i)
    if static
        ηt = Uxt = Uyt = dxUxt = dyUxt = dxUyt = dyUyt = dxDUt = dyDUt = jzero()
    else
        ut = ∂t(u)
        pt(i) = _part(dut, i)
        ηt    = _js(ut[1], pt(1))
        Uxt   = _js(ut[2], pt(2));       Uyt   = _js(ut[3], pt(3))
        dxUxt = _jdx(ut[2], pt(2));      dyUxt = _jdy(ut[2], pt(2))
        dxUyt = _jdx(ut[3], pt(3));      dyUyt = _jdy(ut[3], pt(3))
        #  ∂ₐ(∇·u̇) from the element-wise Hessians of u̇
        dxDUt = _jh(ut[2], pt(2), 1, 1) + _jh(ut[3], pt(3), 1, 2)
        dyDUt = _jh(ut[2], pt(2), 1, 2) + _jh(ut[3], pt(3), 2, 2)
    end
    #  Hessian of H = bed Hessian (analytic) + element-wise Hessian of η
    if prob.flat_bed
        hxx = jzero(); hxy = jzero(); hyy = jzero()
    else
        bxx, bxy, byy = alg_bed_hessian(C.d_cf)
        hxx = Jet(bxx); hxy = Jet(bxy); hyy = Jet(byy)
    end
    Hxx = hxx + _jh(u[1], p(1), 1, 1)
    Hxy = hxy + _jh(u[1], p(1), 1, 2)
    Hyy = hyy + _jh(u[1], p(1), 2, 2)
    F = (; ηt, Uxt, Uyt, dxUxt, dyUxt, dxUyt, dyUyt, dxDUt, dyDUt, hxx, hxy, hyy, Hxx, Hxy, Hyy)
    if prob.nl_pressure
        Gx = _js(u[4], p(4));  Gy = _js(u[5], p(5))
        Fx = _js(u[6], p(6));  Fy = _js(u[7], p(7))
        A = (; Gx, Gy, Fx, Fy,
               dxGx = _jdx(u[4], p(4)), dyGx = _jdy(u[4], p(4)),
               dxGy = _jdx(u[5], p(5)), dyGy = _jdy(u[5], p(5)),
               dxFx = _jdx(u[6], p(6)), dyFx = _jdy(u[6], p(6)),
               dxFy = _jdx(u[7], p(7)), dyFy = _jdy(u[7], p(7)),
               Uxxx = _jh(u[2], p(2), 1, 1), Uxxy = _jh(u[2], p(2), 1, 2), Uxyy = _jh(u[2], p(2), 2, 2),
               Uyxx = _jh(u[3], p(3), 1, 1), Uyxy = _jh(u[3], p(3), 1, 2), Uyyy = _jh(u[3], p(3), 2, 2))
        return merge(C, F, A)
    end
    return merge(C, F)
end

"""
    vms_strong_residuals(prob, vc, t, trian, F) -> (R_c, R_mx, R_my)   (jets)

The strong residuals of continuity and of the two momentum equations at the fields `F`
(`vms_fields`). Term by term the strong counterpart of `global_residual`, nonlinear regime
(file header; MIXED_VMS_PLAN.md §3.1).
"""
function vms_strong_residuals(prob::BALFEMProblem, vc, t::Real, trian, F)
    g  = prob.g
    (; η, ηt, Ux, Uy, Uxt, Uyt, ex, ey, dxUx, dyUx, dxUy, dyUy, dxUxt, dyUxt, dxUyt, dyUyt,
       dxDUt, dyDUt, dhx, dhy, hxx, hxy, hyy, H, Hx, Hy, Hxx, Hxy, Hyy, DU, bf, S) = F
    mu_cf  = CellField(prob.mu_sponge, trian)
    src_cf = CellField(x -> prob.wm_src(x, t), trian)

    # ---- continuity:  η_t + ∇·(H ū) − S_src + μ η ------------------------------------------
    Rc = ηt + alg_dot(prob.Φ, S) + mu_cf*η - src_cf

    # ---- the linear pressure operator 𝓛(u̇) and its gradient --------------------------------
    DUt  = dxUxt + dyUyt
    L1   = (-1.0)*(dhx*Uxt + dhy*Uyt)
    L2   = Hx*Uxt + Hy*Uyt
    L3   = (-1.0)*(H*DUt + L2)
    dxL1 = (-1.0)*(hxx*Uxt + dhx*dxUxt + hxy*Uyt + dhy*dxUyt)
    dyL1 = (-1.0)*(hxy*Uxt + dhx*dyUxt + hyy*Uyt + dhy*dyUyt)
    dxL2 = Hxx*Uxt + Hx*dxUxt + Hxy*Uyt + Hy*dxUyt
    dyL2 = Hxy*Uxt + Hx*dyUxt + Hyy*Uyt + Hy*dyUyt
    dxL3 = (-1.0)*(Hx*DUt + H*dxDUt + dxL2)
    dyL3 = (-1.0)*(Hy*DUt + H*dyDUt + dyL2)
    P1, P2, P3 = prob.P
    sP   = alg_mul(P1, L1) + alg_mul(P2, L2) + alg_mul(P3, L3)
    dxsP = alg_mul(P1, dxL1) + alg_mul(P2, dxL2) + alg_mul(P3, dxL3)
    dysP = alg_mul(P1, dyL1) + alg_mul(P2, dyL2) + alg_mul(P3, dyL3)
    LA   = alg_mul(prob.Av[1], L1) + alg_mul(prob.Av[2], L2) + alg_mul(prob.Av[3], L3)
    LK   = alg_mul(prob.Kv[1], L1) + alg_mul(prob.Kv[2], L2) + alg_mul(prob.Kv[3], L3)

    # ---- advection ----------------------------------------------------------------------
    TMx = alg_outer(Ux, dxUx) + alg_outer(Uy, dyUx)
    TMy = alg_outer(Ux, dxUy) + alg_outer(Uy, dyUy)

    #  one momentum direction: (uₐ, u̇ₐ, ∂ₐh, ∂ₐH, ∂ₐη, ∂ₐ(sP), TMₐ)
    function mom(Ua, Uat, dah, Ha, ea, dasP, TMa)
        local R
        R = H*alg_mul(prob.Mv, Uat)                                  # acceleration
        R = R + alg_svec(prob.Φ, g*(H*ea))                           # gravity  g H ∂ₐη Φ
        R = R + (2.0*(H*Ha))*sP + (H*H)*dasP                         # ∂ₐ(H² Σ P𝓛)
        R = R - H*(dah*LA + Ha*LK)                                   # slope package
        R = R + mu_cf*alg_mul(prob.Mv, Ua)                           # sponge
        R = R + H*alg_dc3(prob.M3, TMa) + alg_dc3(prob.G3, alg_outer(S, Ua))   # advection
        return R
    end
    Rx = mom(Ux, Uxt, dhx, Hx, ex, dxsP, TMx)
    Ry = mom(Uy, Uyt, dhy, Hy, ey, dysP, TMy)

    # ---- relaxation zone ------------------------------------------------------------------
    if prob.relax_bc
        tg     = prob.relax_tg
        mug_cf = CellField(prob.relax_mu, trian)
        Rc = Rc + mug_cf*(η - CellField(x -> tg.eta(x, t), trian))
        Rx = Rx + mug_cf*alg_mul(prob.Mv, Ux - CellField(x -> tg.ux(x, t), trian))
        Ry = Ry + mug_cf*alg_mul(prob.Mv, Uy - CellField(x -> tg.uy(x, t), trian))
    end

    # ---- nonlinear pressure 𝓝 (mixed: ∇𝖲 → 𝖦, ∇𝖻 → 𝖥) ---------------------------------------
    if prob.nl_pressure
        Nx, Ny = vms_strong_nlp(prob, vc, trian, F)
        Rx = Rx + Nx
        Ry = Ry + Ny
    end

    # ---- analytic MMS forcing ---------------------------------------------------------------
    if prob.mms_src !== nothing
        Sm = prob.mms_src
        Rc = Rc - CellField(x -> Sm.Seta(x, t), trian)
        Rx = Rx - CellField(x -> Sm.Sx(x, t), trian)
        Ry = Ry - CellField(x -> Sm.Sy(x, t), trian)
    end
    return Rc, Rx, Ry
end

"""
    vms_strong_nlp(prob, vc, trian, F) -> (𝓝x, 𝓝y)   (jets)

Strong form of the nonlinear pressure terms of the momentum equations in the mixed formulation:

    − (∂ₐh NA₆₈ + ∂ₐH NK₆₈) − H (∂ₐh NA₃ + ∂ₐH NK₃) + ∂ₐ(H NP₆₈ + H² NP₃)      components 3,6,7,8
    − H ∂ₐh (W_A ⊙ GU + 𝓐⁽²⁾ ⊙ SD + 𝓐⁽⁴⁾ ⊙ N4)                                 bed-slope half
    − H ∂ₐH NK + ∂ₐ(H² NP)                                                     Class III
"""
function vms_strong_nlp(prob::BALFEMProblem, vc, trian, F)
    (; Ux, Uy, dxUx, dyUx, dxUy, dyUy, dhx, dhy, hxx, hxy, hyy, H, Hx, Hy, DU, bf, S,
       Gx, Gy, Fx, Fy, dxGx, dyGx, dxGy, dyGy, dxFx, dyFx, dxFy, dyFy,
       Uxxx, Uxxy, Uxyy, Uyxx, Uyxy, Uyyy) = F
    A3, K3, P3 = prob.A3, prob.K3, prob.P3
    c3(W, T2, T4, gu, sd, n4) = alg_dc3(W, gu) + alg_dc3(T2, sd) + alg_dc3(T4, n4)

    af   = dhx*Ux + dhy*Uy                                   # 𝖺 = 𝗎·∇h
    dxaf = hxx*Ux + dhx*dxUx + hxy*Uy + dhy*dxUy
    dyaf = hxy*Ux + dhx*dyUx + hyy*Uy + dhy*dyUy
    M6 = (-1.0)*alg_outer(af, S)
    M7 = alg_outer(bf, S)
    M8 = (-1.0)*alg_outer(S, S)
    N3 = (-1.0)*(alg_outer(Ux, dxaf) + alg_outer(Uy, dyaf))
    NA68 = alg_dc3(A3[6], M6) + alg_dc3(A3[7], M7) + alg_dc3(A3[8], M8)
    NK68 = alg_dc3(K3[6], M6) + alg_dc3(K3[7], M7) + alg_dc3(K3[8], M8)
    NP68 = alg_dc3(P3[6], M6) + alg_dc3(P3[7], M7) + alg_dc3(P3[8], M8)
    NA3  = alg_dc3(A3[3], N3);  NK3 = alg_dc3(K3[3], N3);  NP3 = alg_dc3(P3[3], N3)

    dxDU = Uxxx + Uyxy                                       # ∂ₓ(∇·𝗎)
    dyDU = Uxxy + Uyyy
    GU = alg_outer(Gx, Ux) + alg_outer(Gy, Uy)
    SD = alg_outer(S, DU)
    N4 = alg_outer(Ux, Fx) + alg_outer(Uy, Fy)
    NKc = c3(prob.WK3, K3[2], K3[4], GU, SD, N4)
    NAc = c3(vc.WA3,   A3[2], A3[4], GU, SD, N4)
    NPc = c3(prob.WP3, P3[2], P3[4], GU, SD, N4)

    #  third derivatives of the bed enter ∂ₐ(∂_b𝖺) only (variable bed)
    if vc.third === nothing
        t3 = (jzero(), jzero(), jzero(), jzero())
    else
        t3 = map(Jet, vc.third)
    end
    hxxx, hxxy, hxyy, hyyy = t3
    dxdxaf = hxxx*Ux + 2.0*(hxx*dxUx) + dhx*Uxxx + hxxy*Uy + 2.0*(hxy*dxUy) + dhy*Uyxx
    dxdyaf = hxxy*Ux + hxx*dyUx + hxy*dxUx + dhx*Uxxy + hxyy*Uy + hxy*dyUy + hyy*dxUy + dhy*Uyxy
    dydyaf = hxyy*Ux + 2.0*(hxy*dyUx) + dhx*Uxyy + hyyy*Uy + 2.0*(hyy*dyUy) + dhy*Uyyy

    #  one direction: (∂ₐh, ∂ₐH, 𝖦ₐ, 𝖥ₐ, ∂ₐ𝖺, ∂ₐ𝖴x, ∂ₐ𝖴y, ∂ₐ∂ₓ𝖺, ∂ₐ∂ᵧ𝖺, ∂ₐ𝖦x, ∂ₐ𝖦y, ∂ₐ𝖥x, ∂ₐ𝖥y, ∂ₐD)
    function dir(dah, Ha, Ga, Fa, daaf, daUx, daUy, dadxaf, dadyaf, daGx, daGy, daFx, daFy, daDU)
        local daM6, daM7, daM8, daN3, daNP68, daNP3, daGU, daSD, daN4, daNPc, R
        daM6 = (-1.0)*(alg_outer(daaf, S) + alg_outer(af, Ga))
        daM7 = alg_outer(Fa, S) + alg_outer(bf, Ga)
        daM8 = (-1.0)*(alg_outer(Ga, S) + alg_outer(S, Ga))
        daN3 = (-1.0)*(alg_outer(daUx, dxaf) + alg_outer(Ux, dadxaf) +
                       alg_outer(daUy, dyaf) + alg_outer(Uy, dadyaf))
        daNP68 = alg_dc3(P3[6], daM6) + alg_dc3(P3[7], daM7) + alg_dc3(P3[8], daM8)
        daNP3  = alg_dc3(P3[3], daN3)
        daGU = alg_outer(daGx, Ux) + alg_outer(Gx, daUx) + alg_outer(daGy, Uy) + alg_outer(Gy, daUy)
        daSD = alg_outer(Ga, DU) + alg_outer(S, daDU)
        daN4 = alg_outer(daUx, Fx) + alg_outer(Ux, daFx) + alg_outer(daUy, Fy) + alg_outer(Uy, daFy)
        daNPc = c3(prob.WP3, P3[2], P3[4], daGU, daSD, daN4)
        #  components 3,6,7,8
        R = (-1.0)*(dah*NA68 + Ha*NK68) - H*(dah*NA3 + Ha*NK3)
        R = R + Ha*NP68 + H*daNP68 + (2.0*(H*Ha))*NP3 + (H*H)*daNP3
        #  bed-slope half and Class III of components 1,2,4,5
        R = R - (H*dah)*NAc - (H*Ha)*NKc + (2.0*(H*Ha))*NPc + (H*H)*daNPc
        return R
    end
    Nx = dir(dhx, Hx, Gx, Fx, dxaf, dxUx, dxUy, dxdxaf, dxdyaf, dxGx, dxGy, dxFx, dxFy, dxDU)
    Ny = dir(dhy, Hy, Gy, Fy, dyaf, dyUx, dyUy, dxdyaf, dydyaf, dyGx, dyGy, dyFx, dyFy, dyDU)
    return Nx, Ny
end

"""
    vms_strong_contrib(prob, t, u, v, trian, dΩh; vc=nothing)

`∫ q R_c + Wx·R_mx + Wy·R_my`: the strong residuals tested with the physics test functions.
For a smooth state this equals the weak residual of the three physics rows up to the
inter-element jumps (test_mixed_vms.jl G3). Verification only.
"""
function vms_strong_contrib(prob::BALFEMProblem, t::Real, u, v, trian, dΩh; vc = nothing)
    vc = vc === nothing ? (has_vms(prob) ? prob.vms[] : vms_ctx(prob, trian; p_u = 2)) : vc
    F  = vms_fields(prob, trian, u)
    Rc, Rx, Ry = vms_strong_residuals(prob, vc, t, trian, F)
    return jint(jval, v[1]*Rc + (v[2] ⋅ Rx) + (v[3] ⋅ Ry), dΩh)
end

_tau_a(c, k) = (h, u0, s) -> c / sqrt(k + (2.0/h)^2 * (u0*u0 + s))
_dtau_a(c, k) = (h, u0, s) -> -0.5 * c * (2.0/h)^2 * (k + (2.0/h)^2 * (u0*u0 + s))^(-1.5)

"""
    vms_contrib(prob, t, u, v, trian, dΩh; du=nothing, dut=nothing)

The residual-based stabilisation terms (file header). No keyword: the residual. `du`: the exact
Jacobian in the direction of the unknowns. `dut`: the exact Jacobian in the direction of u̇.
Returns `nothing` when the requested part is identically zero.
"""
function vms_contrib(prob::BALFEMProblem, t::Real, u, v, trian, dΩh; du = nothing, dut = nothing)
    vc = prob.vms[]
    vc === nothing && error("vms_contrib: no stabilisation attached (attach_vms!)")
    (du !== nothing && dut !== nothing) &&
        error("vms_contrib: seed either `du` or `dut`, not both")
    pick = (du === nothing && dut === nothing) ? jval : jtan
    F  = vms_fields(prob, trian, u; du = du, dut = dut, static = vc.oss)
    Rc, Rx, Ry = vms_strong_residuals(prob, vc, t, trian, F)
    (; H, Ux, Uy) = F
    q, Wx, Wy = v[1], v[2], v[3]
    g  = prob.g
    r  = nothing

    # ---- OSS: projection rows ∫(Π − Rˢ)·χ, then every term reads P⊥(Rˢ) = Rˢ − Π -------------------
    if vc.oss
        i0 = vc.oss_i0
        Πx = _js(u[i0], _part(du, i0));  Πy = _js(u[i0 + 1], _part(du, i0 + 1))
        δ  = _js(u[i0 + 2], _part(du, i0 + 2))
        r  = jint(pick, ((Πx - Rx) ⋅ v[i0]) + ((Πy - Ry) ⋅ v[i0 + 1]) + (δ - Rc)*v[i0 + 2], dΩh)
        Rx = Rx - Πx;  Ry = Ry - Πy;  Rc = Rc - δ
    end

    #  τ_w, τ_η: still-water scales (no dependence on the state)
    tw = (vc.c_wave * 0.5) * Operation(/)(vc.he, vc.c0)
    te = (vc.c_eta  * 0.5) * Operation(/)(vc.he, vc.c0)

    # ---- momentum SUPG:  ∫ Σₐ (N_b ∂_b Wₐ)·Zₐ ------------------------------------------------
    Zx = jzero(); Zy = jzero()
    if vc.c_wave > 0
        Zx = Zx + tw*alg_mul(vc.Ebar, Rx)
        Zy = Zy + tw*alg_mul(vc.Ebar, Ry)
    end
    if vc.c_adv > 0
        usq = (Ux ⋅ alg_mul(prob.Mv, Ux)) + (Uy ⋅ alg_mul(prob.Mv, Uy))
        u0  = 1.0e-3 * vc.c0
        ta  = jfun(_tau_a(vc.c_adv, vc.kdt), _dtau_a(vc.c_adv, vc.kdt), usq, vc.he, u0)
        Zx  = Zx + ta*alg_mul(vc.Qv, Rx)
        Zy  = Zy + ta*alg_mul(vc.Qv, Ry)
    end
    if vc.c_wave > 0 || vc.c_adv > 0
        Lx = (Ux ⋅ alg_cont1(prob.M3, alg_dx(Wx))) + (Uy ⋅ alg_cont1(prob.M3, alg_dy(Wx)))
        Ly = (Ux ⋅ alg_cont1(prob.M3, alg_dx(Wy))) + (Uy ⋅ alg_cont1(prob.M3, alg_dy(Wy)))
        r  = cadd(r, jint(pick, (Lx ⋅ Zx) + (Ly ⋅ Zy), dΩh))
    end

    # ---- continuity PSPG:  ∫ τ_w Σₐ ∂ₐq (e·R_m,a) ------------------------------------------------
    if vc.c_wave > 0
        r = cadd(r, jint(pick, tw*(alg_dx(q)*alg_dot(vc.evec, Rx) + alg_dy(q)*alg_dot(vc.evec, Ry)), dΩh))
    end

    # ---- continuity SUPG and momentum LSIC, both on R_c ----------------------------------------
    if vc.c_eta > 0
        ubx = alg_dot(prob.Φ, Ux);  uby = alg_dot(prob.Φ, Uy)
        DW  = alg_dx(Wx) + alg_dy(Wy)
        r = cadd(r, jint(pick, te*((ubx*alg_dx(q) + uby*alg_dy(q))*Rc), dΩh))
        r = cadd(r, jint(pick, (g*te)*(H*(alg_dot(prob.Φ, DW)*Rc)), dΩh))
    end
    return r
end

"""
    oss_consistent_ic(u0, prob, U, trian, dΩh; t=0.0) -> (u0c, rel_proj_residual)

`u0` with its OSS projection fields REPLACED by the L² projections of the static residuals at the
state of `u0` (call after `mixed_consistent_ic` on a mixed layout: Rˢ reads 𝖦, 𝖥):

    ∫ Πₐ·χ = ∫ Rˢ_m,a·χ          ∫ δ ϑ = ∫ Rˢ_c ϑ

solved field by field against each projection space's Gram matrix. An inconsistent start is not
absorbed by Crank–Nicolson (the algebraic error alternates sign every step), so a hot start needs
it. `rel_proj_residual` is the projection block of `vms_contrib` at the returned state relative to
its value with the projections zeroed (≈ 1e-12). Sequential only.
"""
function oss_consistent_ic(u0, prob::BALFEMProblem, U, trian, dΩh; t::Real = 0.0)
    has_oss(prob) || error("oss_consistent_ic: no OSS stabilisation attached (attach_vms!(…; oss=true))")
    vc = prob.vms[]
    nf = length(U.spaces)
    i0 = vc.oss_i0
    nf == i0 + OSS_N_FIELDS - 1 ||
        error("oss_consistent_ic: expected $(i0 + OSS_N_FIELDS - 1) fields (projections at $i0…), got $nf")
    F  = vms_fields(prob, trian, u0; static = true)
    Rc, Rx, Ry = vms_strong_residuals(prob, vc, t, trian, F)
    srcs  = (jval(Rx), jval(Ry), jval(Rc))
    nfree = [num_free_dofs(U.spaces[i]) for i in 1:nf]
    offs  = cumsum([0; nfree])
    vals  = copy(get_free_dof_values(u0))
    for (k, src) in enumerate(srcs)
        i  = i0 + k - 1
        Ui = U.spaces[i]
        Vi = Ui isa TrialFESpace ? Ui.space : Ui
        A  = assemble_matrix((a, b) -> ∫(a ⋅ b) * dΩh, Ui, Vi)
        b  = assemble_vector(χ -> ∫(src ⋅ χ) * dΩh, Vi)
        vals[(offs[i] + 1):offs[i + 1]] .= A \ b
    end
    u0c = FEFunction(U, vals)

    #  gate: the projection rows of the assembled stabilisation must now vanish
    V   = MultiFieldFESpace([Ui isa TrialFESpace ? Ui.space : Ui for Ui in U.spaces])
    ip  = (offs[i0] + 1):offs[end]
    tc(uh) = TransientCellField(uh, (uh,))          # u̇ is not read: the OSS terms are static
    rv(uh) = assemble_vector(v -> vms_contrib(prob, t, tc(uh), v, trian, dΩh), V)[ip]
    zv  = copy(vals); zv[ip] .= 0.0
    r0  = norm(rv(FEFunction(U, zv)))
    rel = norm(rv(u0c)) / max(r0, 1e-300)
    return u0c, rel
end
