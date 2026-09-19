# ==============================================================
#  test_class3_residual_parity.jl — the reduction must change NOTHING
#
#  The Class-III algebraic reduction (NEW_TREATMENT.md §A.2) is an EXACT identity, so
#  with an identical frozen projection pair the REDUCED assembly must reproduce the
#  DIRECT four-object assembly to round-off, as ASSEMBLED RESIDUAL VECTORS — not merely
#  as a tensor identity (that is test_class3_reduction.jl's job).
#
#  This is a self-consistency test and that is the right kind here (rule 28): the claim
#  under test is an algebraic identity between two implementations of the same operator,
#  not a physical claim. Its value is entirely in being able to FAIL — G4 of
#  test_class3_reduction.jl establishes that the underlying comparison can detect the
#  transpose error, and G3 below adds a mutation control at the residual level.
#
#  Q3/Q2 throughout: at Q2/Q1 the surface Hessian is identically zero in 1-D, so half of
#  what Class III carries would be structurally absent from the comparison.
#
#  RUN:  julia --project=. test/test_class3_residual_parity.jl
# ==============================================================

using GridapBALFEM
using Gridap
using Gridap.TensorValues
using LinearAlgebra, Printf

println("=" ^ 64)
println("  test_class3_residual_parity.jl — reduced ≡ direct (assembled vectors)")
println("=" ^ 64)

n_pass = 0; n_fail = 0
function check(name, cond)
    global n_pass, n_fail
    if cond; println("  PASS  $name"); n_pass += 1
    else;    println("  FAIL  $name"); n_fail += 1; end
end

"""
Assemble the 𝓚 and 𝓟 Class-III blocks BOTH ways on one state and return
(‖reduced − direct‖∞ / ‖direct‖∞, ‖direct‖∞).
"""
function parity(M, p, d_fun, label)
    vert = assemble_vertical_tensors(M, p, Vector{Float64}(resolve_cbdy(M, nothing, p)))
    Nσ   = vert.N_dof
    Lx, Ly = 4.0, 2.0
    model, trian = build_horizontal_model(((0.0, Lx), (0.0, Ly)), (8, 4))
    dΩh  = Measure(trian, 14)
    U, V = build_fe_spaces(model, 3, Nσ; y_wall_bc=:open, p_eta=2)

    prob = build_problem(vert; h_bathy=d_fun, regime=:nonlinear,
                               nl_pressure=:full, flat_bed=false)

    # a smooth, deliberately asymmetric state (asymmetry catches k/j slot swaps)
    ax = [0.030, -0.011, 0.021, 0.014, -0.019][1:Nσ]
    ay = [0.017,  0.026, -0.013, -0.022, 0.011][1:Nσ]
    ujx(j) = x -> ax[j]*x[1]*(Lx - x[1])
    ujy(j) = x -> ay[j]*x[2]*(Ly - x[2])
    eta_f  = x -> 0.004*x[1]^2 - 0.003*x[1]*x[2] + 0.002*x[2]^2 + 0.006*x[1]
    stackf(fs) = x -> VectorValue(ntuple(j -> fs[j](x), Nσ)...)
    uh = interpolate_everywhere([eta_f, stackf([ujx(j) for j in 1:Nσ]),
                                        stackf([ujy(j) for j in 1:Nσ])], U)

    # ONE projection pair, shared by both assemblies — the comparison is of the
    # ASSEMBLY, so the projections must be bit-identical on both sides.
    ctx = build_nlp_ctx(model, 3, Nσ, trian, dΩh)
    update_nlp_state!(prob, ctx, uh)
    st = prob.nlp_state[]

    η = uh[1]; Ux = uh[2]; Uy = uh[3]
    d_cf = CellField(d_fun, trian); H = d_cf + η
    dHx = alg_dx(d_cf) + alg_dx(η); dHy = alg_dy(d_cf) + alg_dy(η)
    DU  = alg_dx(Ux) + alg_dy(Uy)
    S   = H*DU + (dHx*Ux + dHy*Uy)

    rdir = assemble_vector(V) do v
        N1, N2, N4, N5 = nlp_frozen_N(Ux, Uy, S, DU, st.piS, st.pib)
        nlp_gradH_frozen_contrib(prob, H, dHx, dHy, v[2], v[3], N1, N2, N4, N5, dΩh) +
        nlp_P_frozen_contrib(prob, H, alg_dx(v[2]) + alg_dy(v[3]), N1, N2, N4, N5, dΩh)
    end
    rred = assemble_vector(V) do v
        GU, SD, N4 = nlp_class3_reduced_fields(Ux, Uy, S, DU, st.piS, st.pib)
        nlp_gradH_reduced_contrib(prob, H, dHx, dHy, v[2], v[3], GU, SD, N4, dΩh) +
        nlp_P_reduced_contrib(prob, H, alg_dx(v[2]) + alg_dy(v[3]), GU, SD, N4, dΩh)
    end

    scale = maximum(abs, rdir)
    rel   = maximum(abs, rred .- rdir) / scale
    @printf("  %-10s Nσ=%d  ‖direct‖∞=%.3e  rel diff=%.3e\n", label, Nσ, scale, rel)
    return rel, scale, prob, (Ux, Uy, S, DU, H, dHx, dHy, st, dΩh, V)
end

# ---- G1: flat-topped (constant) bed --------------------------------------------
rel1, sc1, _, _ = parity(2, 1, x -> 3.5, "P1LFE-2 flat")
check("G1 reduced ≡ direct, constant bed, Q3/Q2  (rel $(round(rel1, sigdigits=3)))",
      rel1 < 1e-12 && sc1 > 1e-12)

# ---- G2: sloping bed — the 𝓐 IBP path is live alongside and must be undisturbed --
rel2, sc2, prob2, ctxs = parity(2, 1, x -> 3.5 - 0.012*x[1]^2 + 0.008*x[1]*x[2],
                                "P1LFE-2 slope")
check("G2 reduced ≡ direct, sloping bed  (rel $(round(rel2, sigdigits=3)))",
      rel2 < 1e-12 && sc2 > 1e-12)

# ---- G3: a different layer count — the identity is index algebra ----------------
rel3, sc3, _, _ = parity(3, 1, x -> 3.5 - 0.012*x[1]^2, "P1LFE-3 slope")
check("G3 reduced ≡ direct at Nσ=4  (rel $(round(rel3, sigdigits=3)))",
      rel3 < 1e-12 && sc3 > 1e-12)

# ---- G4 MUTATION CONTROL: perturb W and the parity MUST break -------------------
#  Without this, G1–G3 would also pass if both sides were silently zero, or if the
#  comparison had no resolution (rule 30: a bounds check on the right configuration
#  is not a value check; rule 38g: a null result needs an effect size).
let (Ux, Uy, S, DU, H, dHx, dHy, st, dΩh, V) = ctxs
    Nσ = num_free_dofs(V[2]) == 0 ? 3 : 3      # P1LFE-2
    Wbad = TensorValues.ThirdOrderTensorValue{Nσ,Nσ,Nσ}(
               (1.01*prob2.WK3[i,k,j] for i in 1:Nσ, k in 1:Nσ, j in 1:Nσ)...)
    GU, SD, N4 = nlp_class3_reduced_fields(Ux, Uy, S, DU, st.piS, st.pib)
    good = assemble_vector(v -> nlp_gradH_reduced_contrib(prob2, H, dHx, dHy,
                                    v[2], v[3], GU, SD, N4, dΩh), V)
    bad  = assemble_vector(V) do v
        NK = alg_dc3(Wbad, GU) + alg_dc3(prob2.K3[2], SD) + alg_dc3(prob2.K3[4], N4)
        ∫( (-1.0)*H*( dHx*(v[2] ⋅ NK) + dHy*(v[3] ⋅ NK) ) ) * dΩh
    end
    d = maximum(abs, bad .- good) / maximum(abs, good)
    @printf("  mutation (W → 1.01·W) moves the block by rel %.3e\n", d)
    check("G4 CONTROL: a 1%% mutation of W is DETECTED (rel $(round(d, sigdigits=3)))",
          d > 1e-4)
end

println("\n" * "=" ^ 64)
@printf("  %d passed, %d failed\n", n_pass, n_fail)
println("=" ^ 64)
n_fail == 0 || error("test_class3_residual_parity: $n_fail gate(s) failed")
