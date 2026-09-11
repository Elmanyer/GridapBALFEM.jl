# ==============================================================
#  mpi_smoke.jl — minimal DISTRIBUTED smoke test (2 ranks, seconds of physics)
#
#  Covers `setup_and_run_distributed` → GMRES + Jacobi + Newton and the MPI stack.
#  No sequential test reaches that path, and it is the one the cluster uses, so a
#  green sequential suite says nothing about it.
#
#  It also exercises, on the distributed side:
#    * the Taylor-Hood gate (check_taylor_hood runs on EVERY rank before the
#      σ-tensors, so a bad pairing fails identically everywhere rather than
#      deadlocking one rank against the others),
#    * both forks — Gridap @ fix-transient-multifield-ad and WaveSpec's
#      change_seed! — since `_dist_common.jl` pulls in the sea-state helpers.
#
#  RUN:  ~/.julia/bin/mpiexecjl --project=. -n 2 julia --project=. examples/mpi_smoke.jl
#
#  ⚠ Use mpiexecjl, never the system mpiexec (PMIx version mismatch on this
#  machine), and expect a benign OFI error from MPI_Finalize with exit 143.
# ==============================================================

using GridapBALFEM, Gridap, Printf
include(joinpath(@__DIR__, "distributed", "_dist_common.jl"))

is_rank0() && println("=== MPI smoke: 2 ranks, Q2/Q1 Taylor-Hood ===")

for (reg, nlp) in ((:linear, :none), (:nonlinear, :native))
    diags, _, _ = setup_and_run_distributed(
        cpu_grid = (2, 1), M = 2, p_vertical = 1,
        domain = (0.0, 4.0, 0.0, 4.0), partition = (4, 4),
        p_u = 2, p_eta = 1,                      # Taylor-Hood (rule 2b)
        h_val = 3.5, T_wave = 1.6, A_wave = 0.001,
        x_wm = 1.0, y_wm = nothing,
        sponge_wL = 0.5, sponge_wR = 0.5, mu_max = 5.0,
        T_final = 0.1, dt = 0.05,
        regime = reg, nl_pressure = nlp, flat_bed = true,
        save_every = 0, print_every = 10^6,
        diag_every = 1, diag_csv = false,
        output_dir = mktempdir())
    if is_rank0()
        em = maximum(x.eta_max for x in diags)
        ok = isfinite(em) && em < 1.0
        @printf("  %s  %-9s / %-6s   steps=%d   max|eta|=%.3e\n",
                ok ? "PASS" : "FAIL", reg, nlp, length(diags), em)
    end
end

is_rank0() && println("=== MPI smoke complete ===")
