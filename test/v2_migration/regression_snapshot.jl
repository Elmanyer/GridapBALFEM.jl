# ==============================================================
#  regression_snapshot.jl — the v1 → v2 migration gate (V2_SOLVER_PLAN.md step 0)
#
#  Assembles, at one fixed non-trivial state (u, u̇), the residual VECTOR and the two Jacobian
#  MATRICES (∂R/∂u, ∂R/∂u̇) for every configuration that survives into v2, and either RECORDS them
#  (on the v1 code) or COMPARES against the recording (on the v2 code). Entry by entry, so it is
#  far stronger than a convergence rate: any changed term in any surviving path shows up.
#
#  Configurations (v1 API → v2 API):
#    lin_flat / lin_bed          regime=:linear                         (models 1, 2)
#    nl0_flat / nl0_bed          :nonlinear, nl_pressure=:none → false  (models 3, 4)
#    nl1_flat / nl1_bed          :nonlinear, nl_pressure=:full + broken, BOTH Class-III arms → true
#    nl1_bed_ghost               nl1_bed + :ghostvolume γ = 0.01 (both fields)
#    nl1_bed_jg2                 nl1_bed + :jumpgrad order 2, γ = 2e-3
#    nl0_bed_jg1                 nl0_bed + :jumpgrad order 1, γ = 0.3 (a penalty on a model without 𝓝)
#    nl1_per_ghost               x-periodic box, flat bed, nl1 + ghost 0.01 (the wrap facet)
#    lin_bed_q21 / nl1_bed_q21   the Q2/Q1 pairing
#
#  The API is DETECTED (v1 has the field `nl_pressure_full`), so this one file serves both sides.
#
#  RUN:  julia --project=. test/v2_migration/regression_snapshot.jl record    # on v1 code (step 0)
#        julia --project=. test/v2_migration/regression_snapshot.jl compare   # after each v2 step
#  DATA: output/v2_baseline/*.jls  (gitignored; the recording is only needed during the migration)
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf, Random, Serialization

const MODE = length(ARGS) ≥ 1 ? ARGS[1] : "compare"
MODE in ("record", "compare") || error("usage: regression_snapshot.jl record|compare")
const V1  = :nl_pressure_full in fieldnames(BALFEMProblem)
const DIR = joinpath(@__DIR__, "..", "..", "output", "v2_baseline")
mkpath(DIR)
println("="^74)
@printf("  regression_snapshot.jl — mode %s, API detected: %s\n", MODE, V1 ? "v1" : "v2")
println("="^74)

const vert = assemble_vertical_tensors(2, 1, Vector{Float64}(resolve_cbdy(2, nothing, 1)))
const Nσ   = vert.N_dof
bed(x)  = 3.5 - 0.08 * x[1] + 0.05 * x[2] + 0.03 * sin(1.3 * x[1]) * cos(0.9 * x[2])   # bed Hessian ≠ 0
flat(x) = 3.5

function setup(pu; periodic = false)
    pe = pu - 1
    model, trian = build_horizontal_model(((0.0, 3.0), (0.0, 1.5)), (6, 3); x_periodic = periodic)
    deg = 2 * pu + 4
    dΩ  = Measure(trian, deg)
    U, V = build_fe_spaces(model, pu, Nσ; y_wall_bc = :open, p_eta = pe)
    return (model = model, trian = trian, dΩ = dΩ, U = U, V = V, pu = pu, pe = pe, deg = deg)
end

function state(S, a; seed, kx)
    uh = interpolate_everywhere(
        [x -> a * (0.05 * cos(1.1kx * x[1]) + 0.03 * sin(0.8x[2] + 0.3kx * x[1])),
         x -> VectorValue(ntuple(j -> a * 0.03j * sin(0.7kx * x[1] + 0.2j + 0.4x[2]), Nσ)...),
         x -> VectorValue(ntuple(j -> a * 0.02 * (Nσ + 1 - j) * cos(0.6x[2] + 0.5kx * x[1] - 0.1j), Nσ)...)],
        S.U)
    x = copy(get_free_dof_values(uh))
    x .+= 0.01 * a .* randn(MersenneTwister(seed), length(x))       # rough part: exercises the jumps
    return FEFunction(S.U, x)
end

"Build the problem for a configuration, in whichever API is loaded."
function mkprob(S, c)
    bedf = c.flat ? flat : bed
    if V1
        nlp = c.regime === :linear ? :none : (c.nl ? :full : :none)
        p = build_problem(vert; h_bathy = bedf, regime = c.regime, nl_pressure = nlp, flat_bed = c.flat)
        if c.nl || c.gamma > 0
            attach_skeleton!(p, S.model; broken = c.nl, cip_gamma_u = c.gamma, cip_gamma_eta = c.gamma,
                             cip_order = c.order, stabilization = c.stab, p_u = S.pu, p_eta = S.pe,
                             degree = S.deg)
        end
    else
        p = build_problem(vert; h_bathy = bedf, regime = c.regime, nl_pressure = c.nl, flat_bed = c.flat,
                          model = S.model, quad_degree = S.deg)
        c.gamma > 0 && attach_skeleton!(p, S.model; cip_gamma_u = c.gamma, cip_gamma_eta = c.gamma,
                                        cip_order = c.order, stabilization = c.stab, p_u = S.pu,
                                        p_eta = S.pe, degree = S.deg)
    end
    return p
end

C(name; regime = :nonlinear, nl = false, flat = false, gamma = 0.0, order = 1, stab = :jumpgrad,
  pu = 3, periodic = false) =
    (name = name, regime = regime, nl = nl, flat = flat, gamma = gamma, order = order, stab = stab,
     pu = pu, periodic = periodic)

const CASES = [
    C("lin_flat"; regime = :linear, flat = true),
    C("lin_bed";  regime = :linear),
    C("nl0_flat"; flat = true),
    C("nl0_bed"),
    C("nl1_flat"; nl = true, flat = true),
    C("nl1_bed";  nl = true),
    C("nl1_bed_ghost"; nl = true, gamma = 0.01, stab = :ghostvolume),
    C("nl1_bed_jg2";   nl = true, gamma = 2e-3, order = 2),
    C("nl0_bed_jg1";   gamma = 0.3, order = 1),
    C("nl1_per_ghost"; nl = true, flat = true, gamma = 0.01, stab = :ghostvolume, periodic = true),
    C("lin_bed_q21";   regime = :linear, pu = 2),
    C("nl1_bed_q21";   nl = true, pu = 2),
]

#  SNAP_CASES="lin_flat,nl0_bed" restricts the run to a subset (comma list of names).
const ONLY = haskey(ENV, "SNAP_CASES") ? Set(strip.(split(ENV["SNAP_CASES"], ","))) : nothing
ONLY === nothing || println("  subset: ", join(sort!(collect(ONLY)), ", "))

relerr(a, b) = norm(a - b, Inf) / max(norm(b, Inf), 1e-300)

nfail = 0
for c in CASES
    ONLY === nothing || c.name in ONLY || continue
    t0 = time()
    S  = setup(c.pu; periodic = c.periodic)
    p  = mkprob(S, c)
    uh = state(S, 1.0; seed = 1, kx = 1.0)
    ud = state(S, 0.5; seed = 2, kx = 1.3)
    tu = TransientCellField(uh, (ud,))
    R  = assemble_vector(v -> global_residual(0.0, tu, v, p, S.trian, S.dΩ), S.V)
    Ju = assemble_matrix((du, v) -> jacobian_u(0.0, tu, du, v, p, S.trian, S.dΩ), S.U, S.V)
    Jt = assemble_matrix((du, v) -> jacobian_u_t(0.0, tu, du, v, p, S.trian, S.dΩ), S.U, S.V)
    f  = joinpath(DIR, c.name * ".jls")
    if MODE == "record"
        serialize(f, (R = R, Ju = Ju, Jt = Jt))
        @printf("  recorded %-15s  |R|∞=%.4e  |Ju|∞=%.4e  |Jt|∞=%.4e   [%.0f s]\n",
                c.name, norm(R, Inf), norm(Ju, Inf), norm(Jt, Inf), time() - t0)
    else
        ref = deserialize(f)
        eR = relerr(R, ref.R); eJ = relerr(Ju, ref.Ju); eT = relerr(Jt, ref.Jt)
        ok = max(eR, eJ, eT) < 1e-12
        global nfail += !ok
        @printf("  %s  %-15s  rel diff  R %.2e   ∂R/∂u %.2e   ∂R/∂u̇ %.2e   [%.0f s]\n",
                ok ? "PASS" : "FAIL", c.name, eR, eJ, eT, time() - t0)
    end
    flush(stdout)
end
MODE == "compare" && (println("="^74); @printf("  %d configuration(s) differ\n", nfail))
nfail == 0 || error("regression_snapshot: $nfail configuration(s) changed")
