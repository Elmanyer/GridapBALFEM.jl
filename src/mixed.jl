# ==============================================================
#  mixed.jl — THE MIXED FORMULATION of the Class-III blocks (branch mixed-vms-solver)
#
#  Plan and record: markdown_files/MIXED_VMS_PLAN.md. Ported from v1-solver:src/mixed.jl.
#
#  WHAT IT IS. The 𝓚 (surface-slope) and 𝓟 (leading-pressure) halves of the 𝓝 components
#  {1,2,4,5} need ∇𝖲 and ∇𝖻, second derivatives of the unknowns (𝖲 = ∇·(H𝗎), 𝖻 = 𝗎·∇H). The
#  primal solver (src/broken.jl) takes them cellwise and completes them on the skeleton. Here
#  they are UNKNOWNS:
#
#      𝖦ₐ ≈ ∂ₐ𝖲 ,   𝖥ₐ ≈ ∂ₐ𝖻            MultiField  [η, 𝖴x, 𝖴y, 𝖦x, 𝖦y, 𝖥x, 𝖥y]   (7 fields)
#
#  defined weakly with the integration by parts built in, so no derivative ever lands on 𝖲 or 𝖻:
#
#      ∫ 𝖦ₐ·Ψ + ∫ 𝖲·∂ₐΨ − ∮ 𝖲·Ψ nₐ = 0        ∫ 𝖥ₐ·Ψ + ∫ 𝖻·∂ₐΨ − ∮ 𝖻·Ψ nₐ = 0        ∀Ψ
#
#  Every unknown then carries at most ONE derivative: the weak form is conforming on C⁰ Lagrange
#  elements. There is NO skeleton term and NO penalty on this path.
#
#  ⚠ THE BOUNDARY TERM IS ASSEMBLED, NOT DROPPED. 𝖲 is generally non-zero on a solid wall even
#  where 𝗎·n = 0. v1 measured, against an analytic ∂ₓ𝖲: interior recovery 2.7e-06, boundary cells
#  wrong by a relative 22.4 without it.
#
#  ⚠ CONSISTENT, NOT EXACT. 𝖦, 𝖥 are finite-element approximations of ∇𝖲, ∇𝖻 in their own space.
#
#  ⚠ INDEX-1 DAE. The auxiliary rows carry no time derivative; ∂R/∂u̇ has zero rows there and
#  the stage system stays non-singular through the auxiliary mass block. A non-rest initial state
#  needs consistent auxiliary fields: `mixed_consistent_ic`.
#
#  ⚠ PAIRING. 𝖦, 𝖥 live in a VectorValue{Nσ} Lagrange space of order `p_aux` (default p_u, v1's
#  choice). This pairing has no analysis behind it; quote it with any negative stability result.
#
#  RESIDUAL AND JACOBIAN FROM ONE FORMULA. `mixed_contrib` is written on jets (src/jets.jl): its
#  value is the residual, its tangent the EXACT ∂R/∂(η,𝖴,𝖦,𝖥), including the dependence of the
#  Class-III blocks on η and 𝖴 that v1 left out (quasi-Newton there).
# ==============================================================

"Number of auxiliary VectorValue{Nσ} fields of the mixed layout: 𝖦x, 𝖦y, 𝖥x, 𝖥y."
mixed_n_aux(prob::BALFEMProblem) = prob.formulation === :mixed ? 4 : 0

"""
    attach_mixed!(prob, model; degree) -> prob

Build the boundary objects the auxiliary rows need (`dΓ`, the outward normal split into scalar
CellFields) and store them on `prob.mixed[]`. Called by `build_problem_raw` when
`formulation === :mixed`. A mesh with no boundary (periodic in both directions) carries none.
"""
function attach_mixed!(prob::BALFEMProblem, model; degree::Int)
    prob.formulation === :mixed ||
        error("attach_mixed!: the problem was built with formulation=:$(prob.formulation)")
    Γ  = BoundaryTriangulation(model)
    nb = _num_facets(Γ)
    if nb == 0
        prob.mixed[] = (dΓ = nothing, nx = nothing, ny = nothing, nbfacets = 0)
        return prob
    end
    dΓ = Measure(Γ, degree)
    nΓ = get_normal_vector(Γ)
    nx = Operation(n -> n ⋅ Ex)(nΓ)
    ny = Operation(n -> n ⋅ Ey)(nΓ)
    prob.mixed[] = (dΓ = dΓ, nx = nx, ny = ny, nbfacets = nb)
    return prob
end

# ---- jets of the raw fields (shared with src/vms.jl) ---------------------------------------
"Component `i` of the increment `du`, or `nothing` when no increment is seeded."
_part(du, i::Int) = du === nothing ? nothing : du[i]
"Jet of a raw field `f` with tangent `df` (`nothing` ⇒ zero tangent)."
_js(f, df)  = Jet(f, df === nothing ? JZ : df)
_jdx(f, df) = Jet(alg_dx(f), df === nothing ? JZ : alg_dx(df))
_jdy(f, df) = Jet(alg_dy(f), df === nothing ? JZ : alg_dy(df))
"Jet of the (a,b) element-wise Hessian component of a raw field."
_jh(f, df, a::Int, b::Int) = Jet(alg_hess(f, a, b), df === nothing ? JZ : alg_hess(df, a, b))

"""
    jet_core_fields(prob, trian, u; du=nothing) -> NamedTuple of jets

The fields every mixed and VMS term is built from, exactly as `global_residual` builds them:
`H = d + η`, `∇H`, `D = ∇·𝗎`, `𝖻 = 𝗎·∇H`, `𝖲 = H D + 𝖻`. With `du` (the trial functions of the
unknowns) each jet carries its exact directional derivative. On a flat bed `∂ₐh` is the zero jet,
so every bed-slope term is pruned.
"""
function jet_core_fields(prob::BALFEMProblem, trian, u; du = nothing)
    d_cf = depth_cf(prob, trian)
    d    = Jet(d_cf)
    dhx  = prob.flat_bed ? jzero() : Jet(alg_dx(d_cf))
    dhy  = prob.flat_bed ? jzero() : Jet(alg_dy(d_cf))
    η    = _js(u[1], _part(du, 1))
    Ux   = _js(u[2], _part(du, 2))
    Uy   = _js(u[3], _part(du, 3))
    ex   = _jdx(u[1], _part(du, 1));  ey   = _jdy(u[1], _part(du, 1))
    dxUx = _jdx(u[2], _part(du, 2));  dyUx = _jdy(u[2], _part(du, 2))
    dxUy = _jdx(u[3], _part(du, 3));  dyUy = _jdy(u[3], _part(du, 3))
    H    = d + η
    Hx   = dhx + ex
    Hy   = dhy + ey
    DU   = dxUx + dyUy
    bf   = Hx*Ux + Hy*Uy
    S    = H*DU + bf
    return (; d_cf, d, dhx, dhy, η, Ux, Uy, ex, ey, dxUx, dyUx, dxUy, dyUy, H, Hx, Hy, DU, bf, S)
end

"""
    mixed_contrib(prob, u, v, trian, dΩh; du=nothing)

The terms the mixed formulation adds to the three physics rows of `global_residual`:

  * the four auxiliary rows (volume and boundary parts);
  * the Class-III 𝓚 and 𝓟 blocks fed with the unknowns,
    `GU = Σₐ 𝖦ₐ⊗𝖴ₐ`, `SD = 𝖲⊗D`, `N4 = Σₐ 𝖴ₐ⊗𝖥ₐ`, through the reduced weights `W_K`, `W_P`
    (the same contractions as `_c3_sum`).

`du === nothing`: the residual. With `du` (trial functions of all seven fields): the exact
Jacobian `∂/∂(η,𝖴,𝖦,𝖥)` of these terms in the direction `du`.
"""
function mixed_contrib(prob::BALFEMProblem, u, v, trian, dΩh; du = nothing)
    mx = prob.mixed[]
    mx === nothing &&
        error("mixed_contrib: prob.mixed[] is nothing — build the problem with " *
              "build_problem(…; formulation=:mixed, model, quad_degree)")
    pick = du === nothing ? jval : jtan
    C  = jet_core_fields(prob, trian, u; du = du)
    Gx = _js(u[4], _part(du, 4));  Gy = _js(u[5], _part(du, 5))
    Fx = _js(u[6], _part(du, 6));  Fy = _js(u[7], _part(du, 7))
    Wx, Wy = v[2], v[3]
    Pgx, Pgy, Pfx, Pfy = v[4], v[5], v[6], v[7]
    DW = alg_dx(Wx) + alg_dy(Wy)
    (; H, Hx, Hy, DU, bf, S, Ux, Uy) = C

    # ---- auxiliary rows:  ∫ 𝖦ₐ·Ψ + ∫ 𝖲·∂ₐΨ − ∮ 𝖲·Ψ nₐ   (and the same for 𝖥, 𝖻) -----------------
    r = jint(pick, (Gx ⋅ Pgx) + (S ⋅ alg_dx(Pgx)) + (Gy ⋅ Pgy) + (S ⋅ alg_dy(Pgy)), dΩh)
    r = cadd(r, jint(pick, (Fx ⋅ Pfx) + (bf ⋅ alg_dx(Pfx)) + (Fy ⋅ Pfy) + (bf ⋅ alg_dy(Pfy)), dΩh))
    if mx.dΓ !== nothing
        r = cadd(r, jint(pick, (-1.0)*((S ⋅ Pgx)*mx.nx + (S ⋅ Pgy)*mx.ny), mx.dΓ))
        r = cadd(r, jint(pick, (-1.0)*((bf ⋅ Pfx)*mx.nx + (bf ⋅ Pfy)*mx.ny), mx.dΓ))
    end

    # ---- Class-III blocks on the unknowns --------------------------------------------------
    GU = alg_outer(Gx, Ux) + alg_outer(Gy, Uy)
    SD = alg_outer(S, DU)
    N4 = alg_outer(Ux, Fx) + alg_outer(Uy, Fy)
    NK = alg_dc3(prob.WK3, GU) + alg_dc3(prob.K3[2], SD) + alg_dc3(prob.K3[4], N4)
    NP = alg_dc3(prob.WP3, GU) + alg_dc3(prob.P3[2], SD) + alg_dc3(prob.P3[4], N4)
    r = cadd(r, jint(pick, (-1.0)*(H*(Hx*(NK ⋅ Wx) + Hy*(NK ⋅ Wy))), dΩh))
    r = cadd(r, jint(pick, (-1.0)*((H*H)*(NP ⋅ DW)), dΩh))
    return r
end

"""
    make_initial_conditions_mixed(U, Nσ; eta0_func=nothing, ux0_func=nothing, uy0_func=nothing)

Initial state on the 7-field layout with the auxiliary fields ZERO. That satisfies the auxiliary
rows only from rest; any other start must go through [`mixed_consistent_ic`](@ref).
"""
function make_initial_conditions_mixed(U, Nσ::Int; eta0_func = nothing,
                                       ux0_func = nothing, uy0_func = nothing)
    zvv  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    eta0 = isnothing(eta0_func) ? (x -> 0.0) : eta0_func
    ux0  = isnothing(ux0_func)  ? (x -> zvv) : ux0_func
    uy0  = isnothing(uy0_func)  ? (x -> zvv) : uy0_func
    zf   = x -> zvv
    return interpolate_everywhere([eta0, ux0, uy0, zf, zf, zf, zf], U)
end

"""
    make_initial_conditions_ext(U, prob; eta0_func=nothing, ux0_func=nothing, uy0_func=nothing)

Initial state on any layout of `prob` (3 fields, mixed, OSS): the physics fields from the given
functions, every extra field zero (`extra_field_zeros`). Follow with `mixed_consistent_ic` and/or
`oss_consistent_ic` for a non-rest start.
"""
function make_initial_conditions_ext(U, prob::BALFEMProblem; eta0_func = nothing,
                                     ux0_func = nothing, uy0_func = nothing)
    zvv  = VectorValue(ntuple(_ -> 0.0, prob.Nσ)...)
    eta0 = isnothing(eta0_func) ? (x -> 0.0) : eta0_func
    ux0  = isnothing(ux0_func)  ? (x -> zvv) : ux0_func
    uy0  = isnothing(uy0_func)  ? (x -> zvv) : uy0_func
    return interpolate_everywhere(Any[eta0, ux0, uy0, extra_field_zeros(prob)...], U)
end

"""
    mixed_consistent_ic(u0, prob, U, trian, dΩh) -> (u0c, rel_aux_residual)

`u0` with its auxiliary fields REPLACED by the solution of their own rows at the physical state
`(η, 𝖴x, 𝖴y)` of `u0`:

    ∫ 𝖦ₐ·Ψ = −∫ 𝖲·∂ₐΨ + ∮ 𝖲·Ψ nₐ          ∫ 𝖥ₐ·Ψ = −∫ 𝖻·∂ₐΨ + ∮ 𝖻·Ψ nₐ

solved field by field against the auxiliary Gram matrix. Without it a hot start (the closed box
starts from a wave) violates the constraint by O(1) at t₀ and the first step has to absorb it.

`rel_aux_residual` is the auxiliary block of `mixed_contrib` at the returned state, relative to
its value with the auxiliary fields zeroed (≈ 1e-12): the gate that the field-by-field solve
reproduces the rows the solver integrates. Sequential only.
"""
function mixed_consistent_ic(u0, prob::BALFEMProblem, U, trian, dΩh)
    prob.formulation === :mixed || error("mixed_consistent_ic: the problem is not mixed")
    nf = length(U.spaces)
    nf == 7 + oss_n_fields(prob) ||
        error("mixed_consistent_ic: expected the 7-field layout (+ $(oss_n_fields(prob)) OSS), got $nf fields")
    mx = prob.mixed[]
    C  = jet_core_fields(prob, trian, u0)
    S  = jval(C.S);  bf = jval(C.bf)
    srcs  = [(S, :x), (S, :y), (bf, :x), (bf, :y)]
    nfree = [num_free_dofs(U.spaces[i]) for i in 1:nf]
    offs  = cumsum([0; nfree])
    vals  = copy(get_free_dof_values(u0))
    for (ia, (src, dir)) in enumerate(srcs)
        i  = 3 + ia
        Ua = U.spaces[i]
        Va = Ua isa TrialFESpace ? Ua.space : Ua
        dψ = dir === :x ? alg_dx : alg_dy
        A  = assemble_matrix((gf, ψ) -> ∫(gf ⋅ ψ) * dΩh, Ua, Va)
        b  = if mx.dΓ === nothing
            assemble_vector(ψ -> ∫(-(src ⋅ dψ(ψ))) * dΩh, Va)
        else
            na = dir === :x ? mx.nx : mx.ny
            assemble_vector(ψ -> ∫(-(src ⋅ dψ(ψ))) * dΩh + ∫((src ⋅ ψ) * na) * mx.dΓ, Va)
        end
        vals[(offs[i] + 1):offs[i + 1]] .= A \ b
    end
    u0c = FEFunction(U, vals)

    #  gate: the auxiliary rows of the assembled mixed terms must now vanish
    V   = MultiFieldFESpace([Ua isa TrialFESpace ? Ua.space : Ua for Ua in U.spaces])
    ia  = (offs[4] + 1):offs[8]                     # 𝖦x … 𝖥y (OSS fields, if any, follow)
    rv(uh) = assemble_vector(v -> mixed_contrib(prob, uh, v, trian, dΩh), V)[ia]
    zv  = copy(vals); zv[ia] .= 0.0
    r0  = norm(rv(FEFunction(U, zv)))
    rel = norm(rv(u0c)) / max(r0, 1e-300)
    return u0c, rel
end
