# ==============================================================
#  test_taylor_hood.jl — the horizontal element-pairing gate
#
#  WHY THIS TEST EXISTS. For three days a "nonlinear instability of the model" was
#  chased through ten refuted hypotheses — Jacobian, sponge, boundaries, domain
#  length, CFL, Benjamin–Feir, quadrature aliasing, advection energy conservation.
#  It was none of them. Every failing run used EQUAL-ORDER Q2/Q2 horizontal spaces,
#  while the entire MMS verification campaign ran Q3/Q2: the configuration that blew
#  up was never the configuration that was verified.
#
#  η enters the momentum equation undifferentiated, via ∇·v after the integration by
#  parts, so it plays the pressure role of a Stokes system and the pairing is subject
#  to the inf-sup (LBB) condition. Equal-order continuous spaces admit a spurious
#  checkerboard at λ ≈ 2·dx that refinement does not remove — and the measured growth
#  spectrum peaked at λ ≈ 2–3·dx. On Taylor-Hood the mode does not exist and the
#  refinement signature inverts. CLAUDE.md rules 2b and 12b.
#
#  ⚠ THIS GATE MUST RAISE, NOT WARN. A warning is a line of scrollback in a run that
#  takes hours; the whole point is that a non-Taylor-Hood configuration cannot be
#  launched at all. The test therefore checks the REJECTIONS as carefully as the
#  acceptances, and checks that the message explains the CAUSE — a gate that merely
#  states its constraint teaches the next person nothing.
#
#  RUN:  julia --project=. test/test_taylor_hood.jl
# ==============================================================


using GridapBALFEM
using Gridap
np=0; nf=0
chk(n,c)=(global np,nf; c ? (println("  PASS  $n");np+=1) : (println("  FAIL  $n");nf+=1))
raises(f) = try; f(); false; catch e; e isa ErrorException; end

println("="^62); println("  check_taylor_hood — the horizontal pairing gate"); println("="^62)

# valid Taylor-Hood pairings
for (pu,pe) in ((2,1),(3,2),(4,3))
    chk("Q$pu/Q$pe accepted (Taylor-Hood)", check_taylor_hood(pu,pe) === nothing)
end
# invalid: equal order — the defect that caused the 2026-09 instability
for pu in (2,3,4)
    chk("Q$pu/Q$pu REJECTED (equal order, inf-sup deficient)", raises(()->check_taylor_hood(pu,pu)))
end
# invalid: gap != 1, and eta order 0
chk("Q4/Q2 REJECTED (gap of two)",        raises(()->check_taylor_hood(4,2)))
chk("Q2/Q3 REJECTED (eta above velocity)", raises(()->check_taylor_hood(2,3)))
chk("p_eta=0 REJECTED",                    raises(()->check_taylor_hood(1,0)))

# the gate is reachable through build_fe_spaces, on the real path
vert = assemble_vertical_tensors(2,1,[0.0,0.728,1.0]); Nσ = vert.N_dof
model, trian = build_horizontal_model(((0.0,1.0),(0.0,1.0)),(2,2))
UV = build_fe_spaces(model,2,Nσ)          # default p_eta = p_u-1 = 1 -> must not raise
chk("build_fe_spaces default is Taylor-Hood (p_eta = p_u-1)", UV isa Tuple)
chk("build_fe_spaces REJECTS equal order",
    raises(()->build_fe_spaces(model,2,Nσ; p_eta=2)))

# the error message must name the cause, not just the constraint
msg = try; check_taylor_hood(2,2); "" ; catch e; sprint(showerror,e); end
chk("error message explains WHY (inf-sup / checkerboard)",
    occursin("inf-sup", msg) && occursin("checkerboard", msg))
chk("error message suggests a fix", occursin("Use p_u", msg))

println(); println("  $np PASS, $nf FAIL")
nf>0 && error("check_th_gate: $nf failed")
