# ==============================================================
#  test_restart_distributed.jl — gate for checkpoint + resume on the DISTRIBUTED driver
#
#  The distributed twin of test_restart.jl: a forced flume (interior source, sponges, model 3,
#  Q3/Q2, CN, GMRES + :schwarz) run 2N steps straight through, checkpointing every N, against
#  the same run resumed from its step-N checkpoint. The forcing depends on t, so this also checks
#  that the resumed run evaluates the source at ABSOLUTE time. u_n-only resume: not bit-exact
#  (the Newton starting guess differs; src/checkpoint.jl), gated at the Newton-tolerance level.
#    D1  the straight run wrote one part file per rank and a header at N and 2N;
#    D2  the resumed run ends at step 2N, t_2N;
#    D3  ‖u_resumed − u_straight‖∞ / ‖u_straight‖∞ ≤ 1e-4 (reported);
#    D4  a checkpoint from a different rank count is refused (sequential resume of a 2-part one).
#
#  RUN:  ~/.julia/bin/mpiexecjl --project=. -n 2 julia --project=. test/test_restart_distributed.jl
#  Verdict from the PASS/FAIL lines on rank 0 (rule 35).
# ==============================================================
using GridapBALFEM, Printf, LinearAlgebra, MPI

MPI.Initialized() || MPI.Init()
const RANK0 = MPI.Comm_rank(MPI.COMM_WORLD) == 0
const NR    = MPI.Comm_size(MPI.COMM_WORLD)
NR == 2 || error("run with mpiexecjl -n 2")
const OUT = joinpath(@__DIR__, "..", "output", "restart_gate", "dist")
const N, DT = 5, 0.04
np, nf = 0, 0
chk(n, c) = RANK0 && (global np, nf; c ? (println("  PASS  $n"); np += 1) : (println("  FAIL  $n"); nf += 1); flush(stdout))

kw = (cpu_grid = (2, 1), M = 2, p_vertical = 1, p_u = 3, p_eta = 2,
      domain = (0.0, 20.0, 0.0, 0.25), partition = (80, 1),
      h_val = 3.5, T_wave = 1.6, A_wave = 0.01, x_wm = 8.0, sponge_wL = 4.0, sponge_wR = 4.0, mu_max = 40.0,
      dt = DT, solver_type = :theta, regime = :nonlinear, nl_pressure = false, flat_bed = true,
      y_wall_bc = :wall, precond = :schwarz, nl_tol = 1e-8,
      save_every = 0, print_every = 1, check_every = 0, diag_every = 1, checkpoint_keep = 0)

dA, dB = joinpath(OUT, "straight"), joinpath(OUT, "resumed")
RANK0 && (rm(OUT; force = true, recursive = true); mkpath(OUT))
MPI.Barrier(MPI.COMM_WORLD)

ckpath(d, s) = joinpath(d, "checkpoints", @sprintf("ck_%09d", s))
readall(d, s) = vcat([open(io -> ltoh.(reinterpret(Float64, read(io))), @sprintf("%s.part%d.f64", ckpath(d, s), r))
                      for r in 1:NR]...)

setup_and_run_distributed(; kw..., T_final = 2N * DT, output_dir = dA, checkpoint_every = N)
MPI.Barrier(MPI.COMM_WORLD)
chk("D1 straight run: one part file per rank + header at steps $N and $(2N)",
    all(isfile(ckpath(dA, s) * ".toml") && all(isfile(@sprintf("%s.part%d.f64", ckpath(dA, s), r)) for r in 1:NR)
        for s in (N, 2N)))

setup_and_run_distributed(; kw..., T_final = 2N * DT, output_dir = dB, checkpoint_every = N,
                          restart_from = ckpath(dA, N) * ".toml")
MPI.Barrier(MPI.COMM_WORLD)
hB = read_checkpoint_header(ckpath(dB, 2N) * ".toml")
chk("D2 resumed run ends at step $(hB["step"]), t = $(hB["t"])",
    parse(Int, hB["step"]) == 2N && isapprox(parse(Float64, hB["t"]), 2N * DT; atol = 1e-12))
if RANK0
    vA, vB = readall(dA, 2N), readall(dB, 2N)
    rel = norm(vB - vA, Inf) / norm(vA, Inf)
    @printf("  measured: ‖u_resumed − u_straight‖∞/‖u_straight‖∞ = %.2e at step %d (u_n-only resume)\n", rel, 2N)
    chk(@sprintf("D3 resumed ≈ straight at the Newton-tolerance level (rel %.1e ≤ 1e-4)", rel), rel <= 1e-4)
    refused = try
        setup_and_run(; M = 2, p_vertical = 1, p_u = 3, p_eta = 2, domain = ((0.0, 20.0), (0.0, 0.25)),
                      partition = (80, 1), h_val = 3.5, T_wave = 1.6, A_wave = 0.01, x_wm = 8.0,
                      sponge_wL = 4.0, sponge_wR = 4.0, mu_max = 40.0, dt = DT, T_final = 3N * DT,
                      solver_type = :theta, regime = :nonlinear, nl_pressure = false, flat_bed = true,
                      output_dir = joinpath(OUT, "d4"), restart_from = ckpath(dA, N) * ".toml")
        false
    catch e
        msg = sprint(showerror, e)
        println("    refused: ", first(msg, 160)); occursin("part", msg)
    end
    chk("D4 a 2-part checkpoint is refused by a sequential resume", refused)
    println("\n  $np PASS, $nf FAIL")
end
