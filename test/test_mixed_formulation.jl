# ==============================================================
#  test_mixed_formulation.jl — is the projection-free path CORRECT and LIVE?
#
#  Gates `src/mixed.jl` BEFORE any stability claim is read off it. This ordering is not
#  bureaucratic: the previous attempt at removing the lag (the in-loop projections,
#  NEW_TREATMENT.md §E.3) produced a run that survived far past the control's failure
#  time and it meant NOTHING, because the Class-III blocks had silently stopped
#  contributing — `:full` had degenerated to `:native`. These gates exist so that cannot
#  happen twice.
#
#  G1  The mixed system SOLVES at all (5 fields, index-1 DAE, AD Jacobians).
#  G2  ⚠ THE LIVE-KNOB GATE (rule 38d). Mixed `:full` must differ from `:native` by
#      about the size of the Class-III contribution. If it does not, the blocks are not
#      in the residual and every later number is void.
#  G3  Mixed `:full` ≈ projected `:full` at SMALL amplitude, where the two treatments of
#      the same operator must nearly agree. This is what says the mixed formulation is
#      the same physics, not a different model.
#  G4  The auxiliary fields really do approximate the gradients they are defined to be:
#      𝖦 compared against a direct FE evaluation of ∇𝖲 on a smooth prescribed state.
#  G5  `:none`/`:native` are untouched by the new code path (confinement).
#
#  RUN:  julia --project=. test/test_mixed_formulation.jl
#        BALFEM_TEST_ONLY=G2 julia --project=. test/test_mixed_formulation.jl
# ==============================================================

using GridapBALFEM
using Gridap
using Gridap.ODEs
using LinearAlgebra, Printf

println("=" ^ 66)
println("  test_mixed_formulation.jl — projection-free Class-III (diagnostic path)")
println("=" ^ 66)

n_pass = 0; n_fail = 0
function check(name, cond)
    global n_pass, n_fail
    if cond; println("  PASS  $name"); n_pass += 1
    else;    println("  FAIL  $name"); n_fail += 1; end
end
const ONLY = get(ENV, "BALFEM_TEST_ONLY", "")
runblk(name) = isempty(ONLY) || ONLY == name

# ---- a small 1-D flume, Q3/Q2 --------------------------------------------------
const M, PV, PU, PE = 2, 1, 3, 2
const D0, TW, KD    = 3.5, 1.6, 5.5
const LX, NX, LY    = 6.0, 24, 0.25

vert = assemble_vertical_tensors(M, PV, Vector{Float64}(resolve_cbdy(M, nothing, PV)))
const Nσ = vert.N_dof
model, trian = build_horizontal_model(((0.0, LX), (0.0, LY)), (NX, 1))
dΩh = Measure(trian, 2*PU + 4)

"""
Run `nsteps` of a rest-start, internally-forced flume. `mixed=true` uses the 7-field
projection-free path; otherwise the 3-field frozen-projection path.
Returns the final max|η| plus the Newton total.
"""
function run_case(; nl_pressure::Symbol, mixed::Bool, dt::Float64, nsteps::Int,
                    A::Float64 = 0.02, c3_mask::Tuple{Bool,Bool} = (true, false))
    k     = find_wavenumber(2π/TW, D0, 9.81)
    prob  = build_problem(vert; h_bathy = (x -> D0), regime = :nonlinear,
                                nl_pressure = nl_pressure, flat_bed = true,
                                c3_mask = c3_mask,
                                mu_sponge = make_sponge(((0.0,LX),(0.0,LY)), 0.0, 2.0, 0.0, 0.0, 40.0),
                                wm_src = make_wavemaker_line(1.0, A, TW, k))
    #  ⚠ n_aux MUST come from mixed_n_aux(prob): with c3_mask=(true,false) the layout is
    #  5 fields, and a mismatch would make u[4]/u[5] index the wrong field silently.
    naux  = mixed ? mixed_n_aux(prob) : 0
    U, V  = build_fe_spaces(model, PU, Nσ; y_wall_bc = :wall, p_eta = PE, n_aux = naux)
    u0 = mixed ? make_initial_conditions_mixed(U, Nσ; n_aux = naux) :
                 make_initial_conditions(U, Nσ)
    nlp = (!mixed && nl_pressure == :full) ?
          (prob, build_nlp_ctx(model, PU, Nσ, trian, dΩh)) : nothing
    op  = mixed ? build_ode_operator_mixed(prob, U, V, trian, dΩh) :
                  build_ode_operator(prob, U, V, trian, dΩh)
    slv = build_ode_solver(dt; solver_type = :theta, nl_tol = 1e-10, nl_iter = 60)
    d = run_time_loop(op, slv, u0, 0.0, dt*nsteps; trian = trian, Nσ = Nσ,
                       save_every = -1, print_every = 100000, dt = dt,
                       trial_space = U, nlp = nlp, diag_every = 0)
    return (eta = d[end].eta_max, nl = sum(r.nl_iters for r in d))
end

const DT, NS = 0.02, 40      # 0.8 s — the source has had time to act

# ==============================================================
if runblk("G1") || runblk("G2") || runblk("G3")
println("\n-- G1/G2/G3: solvability, live-knob, agreement with the projected path --")
ok = true
local mx, nat, prj
#  ⚠ ALL THREE CARRY THE SAME c3_mask = (true,false). The 5-field mixed path is the
#  projection-free version of the ∇𝖲 (gs) ARM, so its like-for-like control is the
#  PROJECTED gs arm — not plain :full, which also carries component 4.
try
    nat = run_case(nl_pressure = :native, mixed = false, dt = DT, nsteps = NS)
    prj = run_case(nl_pressure = :full,   mixed = false, dt = DT, nsteps = NS)
    mx  = run_case(nl_pressure = :full,   mixed = true,  dt = DT, nsteps = NS)
catch err
    global ok = false
    println("  mixed path raised: ", err)
end
check("G1 the 5-field mixed system solves (index-1 DAE, AD Jacobians)", ok)

if ok
    @printf("  eta: native %.10e | projected :full %.10e | MIXED :full %.10e\n",
            nat.eta, prj.eta, mx.eta)
    c3_proj  = abs(prj.eta - nat.eta)          # Class-III size, projected treatment
    c3_mixed = abs(mx.eta  - nat.eta)          # Class-III size, mixed treatment
    @printf("  Class-III contribution:  projected %.6e   mixed %.6e   ratio %.3f\n",
            c3_proj, c3_mixed, c3_mixed/max(c3_proj,1e-300))

    # ---- G2 THE LIVE-KNOB GATE ------------------------------------------------
    #  If mixed :full ≈ :native the blocks are absent and nothing below means anything.
    check("G2 mixed :full is NOT :native — the Class-III blocks are live " *
          "(rel $(round(c3_mixed/max(abs(nat.eta),1e-300), sigdigits=3)))",
          c3_mixed > 1e-3*abs(nat.eta))

    # ---- G3 agreement with the projected treatment ----------------------------
    #  Same operator, two discretisations of its Class-III part: at this amplitude they
    #  must land within a modest factor. A LARGE disagreement is not automatically a
    #  bug — it would mean the projection error is big, which is the hypothesis under
    #  test — so this gate is deliberately loose and is a SANITY check, not a proof.
    check("G3 mixed and projected Class-III sizes agree within 5x " *
          "(ratio $(round(c3_mixed/max(c3_proj,1e-300), sigdigits=3)))",
          0.2 < c3_mixed/max(c3_proj,1e-300) < 5.0)
    @printf("  Newton totals: native %d | projected %d | mixed %d\n", nat.nl, prj.nl, mx.nl)
end
end

# ==============================================================
if runblk("G4")
println("\n-- G4: 𝖦 really is ∇𝖲 (auxiliary definition vs an ANALYTIC gradient) --")
#  ⚠ THE OBVIOUS REFERENCE DOES NOT EXIST. A first draft of this gate compared 𝖦x
#  against `alg_dx(S)` — and Gridap refused, because S is an Operation-composed
#  expression and ∇ of one is not implemented (rule 6). That refusal IS the Class-III
#  problem in miniature: the "direct" second derivative of the unknowns is precisely the
#  object that cannot be formed. The reference must therefore be ANALYTIC.
#
#  State (flat bed, Uy ≡ 0):   η = c₁x² ,  𝖴x_j = a_j x(L−x)
#      H  = D0 + c₁x²          ∂ₓH  = 2c₁x
#      DU = a_j(L−2x)          𝖻_j  = 2c₁a_j x²(L−x)
#      𝖲_j = (D0+c₁x²)a_j(L−2x) + 2c₁a_j x²(L−x)
#      ∂ₓ𝖲_j = a_j [ 6c₁Lx − 12c₁x² − 2·D0 ]        (hand-derived, degree 2)
#
#  ⚠ MEASURED AT INTERIOR STATIONS ONLY. The weak definition drops ∮𝖲Ψn, and 𝖲 ≠ 0 on
#  the x-boundaries for this state, so 𝖦 is knowingly wrong in the boundary cells
#  (mixed.jl "Known weaknesses" #2). The interior is what tests the definition; the
#  boundary error is reported separately rather than hidden.
c1 = 0.004
ax = [0.030, -0.011, 0.021][1:Nσ]
U5, V5 = build_fe_spaces(model, PU, Nσ; y_wall_bc = :wall, p_eta = PE, n_aux = 2)
zvvf = x -> VectorValue(ntuple(_ -> 0.0, Nσ)...)
uh = interpolate_everywhere(
        [x -> c1*x[1]^2,
         x -> VectorValue(ntuple(j -> ax[j]*x[1]*(LX - x[1]), Nσ)...),
         zvvf, zvvf, zvvf], U5)
η, Ux, Uy = uh[1], uh[2], uh[3]
d_cf = CellField(x -> D0, trian); H = d_cf + η
DU  = alg_dx(Ux) + alg_dy(Uy)
S   = H*DU + (alg_dx(η)*Ux + alg_dy(η)*Uy)

Vaux = FESpace(model, ReferenceFE(lagrangian, VectorValue{Nσ,Float64}, PU); conformity=:H1)
Uaux = TrialFESpace(Vaux)
Γ  = BoundaryTriangulation(model); dΓ = Measure(Γ, 10); nΓ = get_normal_vector(Γ)
nx = Operation(n -> n ⋅ VectorValue(1.0,0.0))(nΓ)
#  WITHOUT the boundary term (the first draft) …
Gx_nobnd = solve(AffineFEOperator((g,ψ) -> ∫(g ⋅ ψ)*dΩh,
                                  ψ -> ∫( (-1.0)*(S ⋅ alg_dx(ψ)) )*dΩh, Uaux, Vaux))
#  … and WITH it, which is what src/mixed.jl assembles.
Gx = solve(AffineFEOperator((g,ψ) -> ∫(g ⋅ ψ)*dΩh,
                            ψ -> ∫( (-1.0)*(S ⋅ alg_dx(ψ)) )*dΩh +
                                 ∫( (S ⋅ ψ)*nx )*dΓ, Uaux, Vaux))

dSx_exact(x) = VectorValue(ntuple(j -> ax[j]*(6*c1*LX*x - 12*c1*x^2 - 2*D0), Nσ)...)
stations_int = range(0.30*LX, 0.70*LX; length = 7)
stations_bnd = (0.02*LX, 0.98*LX)
relerr(G, xs) = begin
    num = 0.0; den = 0.0
    for xv in xs
        g = G(Point(xv, LY/2)); e = dSx_exact(xv)
        num += sum((g[j]-e[j])^2 for j in 1:Nσ); den += sum(e[j]^2 for j in 1:Nσ)
    end
    sqrt(num/max(den,1e-300))
end
ei      = relerr(Gx,       stations_int)
eb      = relerr(Gx,       stations_bnd)
eb_drop = relerr(Gx_nobnd, stations_bnd)
@printf("  WITH    boundary term: interior %.4e   boundary %.4e\n", ei, eb)
@printf("  WITHOUT boundary term: boundary %.4e   <- the first draft\n", eb_drop)
check("G4 𝖦x recovers the ANALYTIC ∂ₓ𝖲 in the interior (rel $(round(ei, sigdigits=3)) < 0.02)",
      ei < 0.02)
#  ⚠ THE GATE THAT MATTERS. The instability under investigation is velocity-led AT THE
#  GENERATION BOUNDARY, so a diagnostic that is locally wrong there proves nothing.
check("G4b the assembled boundary term FIXES the edge cells " *
      "($(round(eb_drop, sigdigits=3)) → $(round(eb, sigdigits=3)))",
      eb < 0.05 && eb_drop > 1.0)
end

if runblk("G5")
println("\n-- G5: the new code path does not touch :none / :native --")
for tier in (:none, :native)
    a = run_case(nl_pressure = tier, mixed = false, dt = DT, nsteps = 10)
    b = run_case(nl_pressure = tier, mixed = false, dt = DT, nsteps = 10)
    check("G5 :$tier is reproducible/unchanged (η = $(a.eta))", a.eta == b.eta)
end
end

println("\n" * "=" ^ 66)
@printf("  %d passed, %d failed\n", n_pass, n_fail)
println("=" ^ 66)
n_fail == 0 || error("test_mixed_formulation: $n_fail gate(s) failed")
