# ==============================================================
#  test_mixed_jacobian.jl — the coupled mixed Jacobian against AD
#
#  Gates NEW_TREATMENT.md Part G. The mixed Jacobian is
#
#      J = [ A  B ;  C  M ]     A = physics (hand),  M = auxiliary mass,
#                               C = dR_G/d(eta,u),   B = dR_phys/dG
#
#  and C, B were previously OMITTED (block-diagonal), costing 26-30 Newton iterations at
#  Q2/Q1 and 64-88 at Q3/Q2 against the projected path's 6-8.
#
#  ORACLE: AD of `global_residual_mixed` itself (rule 29 — AD is an oracle for the
#  JACOBIAN, never for the residual; here that is exactly the right use).
#
#  G1  the COUPLED Jacobian matches AD far better than the block-diagonal one does.
#  G2  ⚠ the C and B blocks are compared INDIVIDUALLY, not only in the total. A total that
#      agrees while two blocks are wrong with cancelling errors is precisely what a summed
#      comparison cannot see (rule 36).
#  G3  the remaining gap is the DELIBERATE omission of G.3 (Class-III d/d(eta,u)) and must
#      VANISH WITH AMPLITUDE; a gap that does not scale is a bug, not the design (rule 5).
#  G4  coupling=false still reproduces the old block-diagonal matrix exactly (no regression).
#
#  RUN:  julia --project=. test/test_mixed_jacobian.jl
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf

println("="^70); println("  test_mixed_jacobian.jl — coupled mixed Jacobian vs AD"); println("="^70)
np=0; nf=0
chk(n,c)=(global np,nf; c ? (println("  PASS  $n"); np+=1) : (println("  FAIL  $n"); nf+=1))

const M,PV,PU,PE = 2,1,3,2
const D0 = 3.5
const LX,NX,LY = 6.0,8,0.75
vert = assemble_vertical_tensors(M,PV,Vector{Float64}(resolve_cbdy(M,nothing,PV)))
const Nσ = vert.N_dof
model,trian = build_horizontal_model(((0.0,LX),(0.0,LY)),(NX,1))
dΩh = Measure(trian, 2*PU+4)
Γ   = BoundaryTriangulation(model); dΓ = Measure(Γ,10); nΓ = get_normal_vector(Γ)

prob = build_problem(vert; h_bathy=(x->D0 - 0.02*x[1]), regime=:nonlinear,
                           nl_pressure=:full, flat_bed=false, c3_mask=(true,false))
naux = mixed_n_aux(prob)
U,V  = build_fe_spaces(model,PU,Nσ; y_wall_bc=:wall, p_eta=PE, n_aux=naux)

"State at amplitude a, as a free-dof VECTOR; auxiliary fields seeded non-zero so B is exercised."
function statevec(a)
    zvv = VectorValue(ntuple(_->0.0,Nσ)...)
    uh = interpolate_everywhere(
        [x -> a*0.05*cos(1.1x[1]),
         x -> VectorValue(ntuple(j -> a*0.03*j*sin(0.7x[1]+0.2j), Nσ)...),
         x -> zvv,
         x -> VectorValue(ntuple(j -> a*0.02*j*cos(0.9x[1]), Nσ)...),
         x -> zvv], U)
    return copy(get_free_dof_values(uh))
end

#  ⚠ u̇ MUST BE HELD FIXED, INDEPENDENT OF x. A first version built
#  `TransientCellField(uh,(uh,))`, i.e. u̇ = u, so perturbing x perturbed u̇ as well and the
#  finite difference measured (∂R/∂u + ∂R/∂u̇)·d while the hand matrix is ∂R/∂u alone. The
#  C and B gates still passed — the auxiliary rows carry no time derivative — but the FULL
#  directional check then showed a 93 % "error" that was entirely the mass/dispersion terms
#  of ∂R/∂u̇. Rule 32's lesson in miniature: verify what the comparison actually isolates.
const UDOT = FEFunction(U, 0.37 .* statevec(1.0))     # fixed, non-zero, never perturbed

"Residual VECTOR at free-dof vector x, with u̇ frozen."
function resid(x)
    uh = FEFunction(U, x)
    tu = TransientCellField(uh, (UDOT,))
    return assemble_vector(v -> global_residual_mixed(0.0, tu, v, prob, trian, dΩh, dΓ, nΓ), V)
end

"Hand ∂R/∂u MATRIX at x, with the same frozen u̇."
function jac_hand(x; coupling)
    uh = FEFunction(U, x); tu = TransientCellField(uh, (UDOT,))
    assemble_matrix((du,v) -> jacobian_u(0.0,tu,du,v,prob,trian,dΩh) +
                              mixed_coupling_jacobian(prob,tu,du,v,trian,dΩh,dΓ,nΓ;
                                                      coupling=coupling), U, V)
end

#  ⚠ THE ORACLE IS A CENTRAL FINITE DIFFERENCE OF THE RESIDUAL, not AD.
#  A first attempt used Gridap's AD and failed on plumbing (TransientCellField cannot be
#  rebuilt from a tuple of single fields). The FD directional derivative tests exactly the
#  same thing — J·δ against dR/dε — with no dependence on the transient-AD API, and it is
#  an INDEPENDENT oracle rather than a differentiation of the same assembled code.
function fd_dir(x, d; h=1e-7)
    (resid(x .+ h.*d) .- resid(x .- h.*d)) ./ (2h)
end

nu = num_free_dofs(V[1]) + num_free_dofs(V[2]) + num_free_dofs(V[3])   # physics dofs come first
ntot = num_free_dofs(U)
x0 = statevec(1.0)

using Random; Random.seed!(20260922)
dphys = zeros(ntot); dphys[1:nu]        .= randn(nu)            # perturb PHYSICS only -> tests C
daux  = zeros(ntot); daux[nu+1:ntot]    .= randn(ntot-nu)       # perturb AUX only     -> tests B

Jd = jac_hand(x0; coupling=false)
Jc = jac_hand(x0; coupling=true)

println("\n  directional check:  J*d  vs  dR/de  (central FD, h=1e-7)")
for (lbl, d, rows, blk) in (("C = dR_G/d(eta,u)", dphys, (nu+1):ntot, "G rows, physics direction"),
                            ("B = dR_phys/dG",   daux,  1:nu,        "physics rows, aux direction"))
    fd = fd_dir(x0, d)
    ed = maximum(abs, (Jd*d)[rows] .- fd[rows])
    ec = maximum(abs, (Jc*d)[rows] .- fd[rows])
    sc = maximum(abs, fd[rows])
    @printf("    %-20s |dR/de|=%.4e  omitted err=%.4e (rel %.2e)  COUPLED err=%.4e (rel %.2e)\n",
            lbl, sc, ed, ed/max(sc,1e-300), ec, ec/max(sc,1e-300))
    chk("G2 $lbl matches the FD oracle (rel $(round(ec/max(sc,1e-300), sigdigits=3)))",
        sc > 1e-8 && ec/sc < 1e-5)
    chk("G2b $lbl was genuinely MISSING before (omitted err $(round(ed/max(sc,1e-300), sigdigits=3)))",
        ed/max(sc,1e-300) > 0.1)
end

#  ---- G1: the whole Jacobian, both directions at once -----------------------
dall = dphys .+ daux
fd   = fd_dir(x0, dall)
ed   = maximum(abs, Jd*dall .- fd); ec = maximum(abs, Jc*dall .- fd); sc = maximum(abs, fd)
@printf("\n  FULL J*d: |dR/de|=%.4e  block-diagonal err=%.4e  COUPLED err=%.4e  improvement %.1fx\n",
        sc, ed, ec, ed/max(ec,1e-300))
chk("G1 the coupled Jacobian is much closer to the true Jacobian", ec < 0.25*ed)

#  ---- G3: the remaining gap is the DELIBERATE Class-III omission -> scales ---
println("\n  remaining gap vs amplitude (must shrink — it is the Class-III d/d(eta,u) omission):")
ratios = let prev = nothing, rs = Float64[]
    for a in (1.0, 0.5, 0.25)
        x = statevec(a); J = jac_hand(x; coupling=true)
        f = fd_dir(x, dall)
        g = maximum(abs, J*dall .- f)/max(maximum(abs,f),1e-300)
        @printf("    a=%.2f  rel gap = %.4e%s\n", a, g,
                prev === nothing ? "" : @sprintf("   ratio %.2f", prev/g))
        prev === nothing || push!(rs, prev/g)
        prev = g
    end
    rs
end
chk("G3 the residual gap shrinks with amplitude (the deliberate omission, not a bug)",
    all(r -> r > 1.2, ratios))

#  ---- G4: coupling=false reproduces the previous block-diagonal matrix -------
chk("G4 coupling=false is still block-diagonal (no coupling entries)",
    maximum(abs, Jd[(nu+1):ntot, 1:nu]) == 0.0 && maximum(abs, Jd[1:nu, (nu+1):ntot]) == 0.0)

println("\n" * "="^70); @printf("  %d passed, %d failed\n", np, nf); println("="^70)
nf == 0 || error("test_mixed_jacobian: $nf gate(s) failed")
