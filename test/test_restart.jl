# ==============================================================
#  test_restart.jl — gate for checkpoint + resume (src/checkpoint.jl), sequential
#
#  A run of 2N steps straight through (checkpointing every N) against the same run resumed from
#  its step-N checkpoint. The checkpoint stores u_n only (author's choice, 2026-10-31), so the
#  resumed run is NOT bit-exact: Gridap starts the first resumed Newton solve from a different
#  guess (ThetaMethod: u_{n−1}; DIRK: the previous stage slopes), and Newton stops at nl_tol.
#  This test MEASURES that difference and gates it at the Newton-tolerance level.
#
#  Cases (closed x-periodic box, one wavelength, 8 cells, Q3/Q2, A = 0.10, dt = 0.04):
#    model 3 under CN, model 3 under SDIRK_2_2, model 5 (𝓝 + facet terms + CIP) under CN.
#  Gates per case:
#    R1  the straight run wrote checkpoints at N and 2N, with the right step and t;
#    R2  the resumed run starts at step N, t_N, and ends at step 2N, t_2N;
#    R3  ‖u_resumed − u_straight‖∞ / ‖u_straight‖∞ ≤ 1e-4 (Newton-tolerance level, reported);
#    R4  loading the step-N checkpoint reproduces the stored values exactly (Σv, Σv² checked);
#  Refusals (model 3): R5 a different mesh, R6 a different dt without restart_allow_change,
#    R7 the same dt change WITH restart_allow_change resumes; R8 T_final ≤ t_checkpoint.
#
#  RUN:  julia --project=. test/test_restart.jl     (~30 min cold: model 5 compiles)
#  Verdict from the PASS/FAIL lines (rule 35), never from the exit code.
# ==============================================================
using GridapBALFEM, Gridap, Printf, LinearAlgebra

const OUT  = mktempdir(; prefix = "balfem_restart_")
const N    = 5                       # steps per half
const DT   = 0.04
const L    = 3.9720                  # ≈ one model wavelength at T = 1.6 s, d = 3.5 m (value irrelevant)
np, nf = 0, 0
chk(n, c) = (global np, nf; c ? (println("  PASS  $n"); np += 1) : (println("  FAIL  $n"); nf += 1); flush(stdout))

k0 = 2π / L
base_kw(; nl_pressure, solver_type, ncell = 8, dt = DT) = (
    M = 2, p_vertical = 1, p_u = 3, p_eta = 2,
    domain = ((0.0, L), (0.0, L / ncell)), partition = (ncell, 1),
    x_periodic = true, y_wall_bc = :wall, x_wall_bc = false,
    h_val = 3.5, T_wave = 1.6, A_wave = 0.10,
    sponge_wL = 0.0, sponge_wR = 0.0, sponge_wB = 0.0, sponge_wT = 0.0,
    eta0_func = x -> 0.10 * cos(k0 * x[1]),
    ux0_func  = x -> VectorValue(0.05 * cos(k0 * x[1]), 0.10 * cos(k0 * x[1]), 0.20 * cos(k0 * x[1])),
    dt = dt, regime = :nonlinear, nl_pressure = nl_pressure, flat_bed = true,
    solver_type = solver_type, nl_tol = 1e-5,
    save_every = 0, print_every = 1, check_every = 0, diag_every = 1, eta_ref = 0.10)

ckvals(dir, step) = begin
    h = read_checkpoint_header(joinpath(dir, "checkpoints", @sprintf("ck_%09d.toml", step)))
    v = open(h["base"] * ".f64") do io; ltoh.(reinterpret(Float64, read(io))) end
    h, v
end

for (tag, nlp, st) in (("m3_theta", false, :theta), ("m3_sdirk", false, :sdirk), ("m5_theta", true, :theta))
    println("\n=== case $tag ==="); flush(stdout)
    kw = base_kw(; nl_pressure = nlp, solver_type = st)
    dA = joinpath(OUT, tag * "_straight"); dB = joinpath(OUT, tag * "_resumed")
    setup_and_run(; kw..., T_final = 2N * DT, output_dir = dA, checkpoint_every = N, checkpoint_keep = 0)
    hN, vN   = ckvals(dA, N)
    h2N, v2N = ckvals(dA, 2N)
    chk("R1 $tag straight run checkpointed at steps $N and $(2N) (t = $(hN["t"]), $(h2N["t"]))",
        parse(Int, hN["step"]) == N && parse(Int, h2N["step"]) == 2N &&
        isapprox(parse(Float64, hN["t"]), N * DT; atol = 1e-12) &&
        isapprox(parse(Float64, h2N["t"]), 2N * DT; atol = 1e-12))

    setup_and_run(; kw..., T_final = 2N * DT, output_dir = dB, checkpoint_every = N, checkpoint_keep = 0,
                  restart_from = joinpath(dA, "checkpoints", @sprintf("ck_%09d.toml", N)))
    hB, vB = ckvals(dB, 2N)
    chk("R2 $tag resumed run ends at step $(hB["step"]), t = $(hB["t"]) (expected $(2N), $(2N * DT))",
        parse(Int, hB["step"]) == 2N && isapprox(parse(Float64, hB["t"]), 2N * DT; atol = 1e-12) &&
        !isfile(joinpath(dB, "checkpoints", @sprintf("ck_%09d.toml", N))))
    rel = norm(vB - v2N, Inf) / norm(v2N, Inf)
    @printf("  measured: ‖u_resumed − u_straight‖∞/‖u_straight‖∞ = %.2e at step %d (u_n-only resume)\n", rel, 2N)
    chk(@sprintf("R3 %s resumed ≈ straight at the Newton-tolerance level (rel %.1e ≤ 1e-4)", tag, rel), rel <= 1e-4)

    U, _ = build_fe_spaces(build_horizontal_model(kw.domain, kw.partition; x_periodic = true)[1],
                           3, 3; y_wall_bc = :wall, p_eta = 2)
    uN = load_checkpoint_state(hN, U)
    chk("R4 $tag reloaded step-$N state is bit-identical to the stored values",
        get_free_dof_values(uN) == vN)
end

println("\n=== refusals (model 3) ==="); flush(stdout)
kw = base_kw(; nl_pressure = false, solver_type = :theta)
ckN = joinpath(OUT, "m3_theta_straight", "checkpoints", @sprintf("ck_%09d.toml", N))
refused(f) = try f(); false catch e; println("    refused: ", first(sprint(showerror, e), 160)); true end
chk("R5 a different mesh is refused",
    refused(() -> setup_and_run(; base_kw(; nl_pressure = false, solver_type = :theta, ncell = 16)...,
                                T_final = 2N * DT, output_dir = joinpath(OUT, "r5"), restart_from = ckN)))
chk("R6 a different dt is refused without restart_allow_change",
    refused(() -> setup_and_run(; base_kw(; nl_pressure = false, solver_type = :theta, dt = 0.02)...,
                                T_final = 2N * DT, output_dir = joinpath(OUT, "r6"), restart_from = ckN)))
chk("R7 the same dt change resumes WITH restart_allow_change",
    !refused(() -> setup_and_run(; base_kw(; nl_pressure = false, solver_type = :theta, dt = 0.02)...,
                                 T_final = N * DT + 2 * 0.02, output_dir = joinpath(OUT, "r7"),
                                 restart_from = ckN, restart_allow_change = true)))
chk("R8 T_final ≤ t_checkpoint is refused",
    refused(() -> setup_and_run(; kw..., T_final = N * DT, output_dir = joinpath(OUT, "r8"), restart_from = ckN)))

println("\n  $np PASS, $nf FAIL   (scratch: $OUT)")
