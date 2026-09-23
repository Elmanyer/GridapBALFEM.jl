# ==============================================================
#  test_diagnostics_mixed.jl — the run diagnostics must work on the MIXED layout
#
#  `build_run_diagnostics` used to interpolate onto a hardcoded 3-field list, so on the
#  5-/7-field mixed layout (src/mixed.jl) it threw — taking `x_at_max`, diagnostics.csv
#  and the divergence guard with it. A mixed failure then had a TIME but no PLACE, which
#  is the one quantity that would let it be compared with the projected failure (that one
#  pins at the inflow, CLAUDE.md §5.2b).
#
#  G1  3-field layout still works and is UNCHANGED (no regression).
#  G2  5-field (𝖦-only) layout builds diagnostics and returns finite values.
#  G3  7-field (𝖦+𝖥) layout likewise.
#  G4  ⚠ x_at_max is MEANINGFUL, not just finite: on a state whose |η| peaks at a known
#      station, the reported location must be that station. A diagnostic that returns a
#      number is not the same as one that returns the RIGHT number (rule 30).
#  G5  the 3-field and 5-field paths agree on the SAME physical state — the auxiliary
#      fields must not perturb the diagnostics of η and 𝗎.
#
#  RUN:  julia --project=. test/test_diagnostics_mixed.jl
# ==============================================================
using GridapBALFEM, Gridap, LinearAlgebra, Printf

println("="^66); println("  test_diagnostics_mixed.jl — diagnostics on the mixed layout"); println("="^66)
np = 0; nf = 0
chk(n,c) = (global np, nf; c ? (println("  PASS  $n"); np+=1) : (println("  FAIL  $n"); nf+=1))

M, p, PU, PE = 2, 1, 3, 2
vert = assemble_vertical_tensors(M, p, Vector{Float64}(resolve_cbdy(M, nothing, p)))
Nσ = vert.N_dof
LX, LY = 6.0, 0.5
model, trian = build_horizontal_model(((0.0,LX),(0.0,LY)), (12,1))
dΩh = Measure(trian, 8)
prob = build_problem(vert; h_bathy = (x -> 3.5), regime = :nonlinear,
                           nl_pressure = :full, flat_bed = true)

#  η peaks at a KNOWN station: a narrow bump centred on x = 4.0
const XPEAK = 4.0
etaf = x -> 0.05*exp(-((x[1]-XPEAK)/0.4)^2)
zvv  = VectorValue(ntuple(_->0.0, Nσ)...)
uxf  = x -> VectorValue(ntuple(j -> 0.01*j*sin(pi*x[1]/LX), Nσ)...)

function probe(naux)
    U, V = build_fe_spaces(model, PU, Nσ; y_wall_bc=:wall, p_eta=PE, n_aux=naux)
    zf = x -> zvv
    uh = interpolate_everywhere([etaf, uxf, zf, fill(zf, naux)...], U)
    rd = build_run_diagnostics(prob, U, trian, dΩh; eta_ref=0.05, div_factor=20.0)
    return field_diagnostics(rd, uh)
end

r3 = try probe(0) catch e; println("  3-field raised: ", e); nothing end
chk("G1 3-field layout still builds diagnostics (no regression)", r3 !== nothing)

r5 = try probe(2) catch e; println("  5-field raised: ", e); nothing end
chk("G2 5-field (𝖦-only) layout builds diagnostics", r5 !== nothing)

r7 = try probe(4) catch e; println("  7-field raised: ", e); nothing end
chk("G3 7-field (𝖦+𝖥) layout builds diagnostics", r7 !== nothing)

if r5 !== nothing
    @printf("  5-field: eta_max=%.6e  x_at_max=%.4f  u_max=%.6e  mass=%.6e\n",
            r5.eta_max, r5.x_at_max, r5.u_max, r5.mass)
    #  ⚠ THE GATE THAT MATTERS: the location must be RIGHT, not merely finite.
    chk("G4 x_at_max finds the known peak at x=$XPEAK (got $(round(r5.x_at_max,digits=3)))",
        isfinite(r5.x_at_max) && abs(r5.x_at_max - XPEAK) < 0.5)
end

if r3 !== nothing && r5 !== nothing
    d_eta = abs(r5.eta_max - r3.eta_max); d_u = abs(r5.u_max - r3.u_max)
    d_x   = abs(r5.x_at_max - r3.x_at_max)
    @printf("  3-field vs 5-field on the same state: Δeta=%.3e Δu=%.3e Δx=%.3e\n", d_eta, d_u, d_x)
    chk("G5 the auxiliary fields do not perturb the η/𝗎 diagnostics",
        d_eta < 1e-14 && d_u < 1e-14 && d_x < 1e-12)
end

println("\n" * "="^66); @printf("  %d passed, %d failed\n", np, nf); println("="^66)
nf == 0 || error("test_diagnostics_mixed: $nf gate(s) failed")
