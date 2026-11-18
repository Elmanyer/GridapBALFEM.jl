# ==============================================================
#  run_ring_wave_dist.jl — DISTRIBUTED point-source ring waves (algebraic)
#
#  Gaussian POINT-source wavemaker at the basin centre → circular wavefronts,
#  4-side sponge absorption. Full field set to VTK (η, u/v per σ-node, w_s*,
#  p_s*). Validation targets: radial symmetry; A ∝ 1/√r cylindrical spreading.
#
#  LAUNCH (px·py MUST equal -n):
#    BABALFEM_M=2 BABALFEM_PX=4 BABALFEM_PY=4 \
#    ~/.julia/bin/mpiexecjl --project=. -n 16 julia --project=. \
#        GridapBABALFEM.jl/examples/distributed/run_ring_wave_dist.jl
#
#  Case-specific env vars:
#    BABALFEM_L        basin side length [m]      200
#    BABALFEM_D        still-water depth [m]      3.5
#    BABALFEM_TWAVE    wave period [s]            1.6
#    BABALFEM_AWAVE    amplitude [m]              0.001
#    BABALFEM_SPONGE   sponge width, all sides    20
#    BABALFEM_MUMAX    sponge strength            10
#    BABALFEM_PERIODS  run length in periods      20
# ==============================================================

include(joinpath(@__DIR__, "_dist_common.jl"))

M        = genv_i("BABALFEM_M", 2)
px, py   = genv_i("BABALFEM_PX", 2), genv_i("BABALFEM_PY", 2)
L        = genv_f("BABALFEM_L", 200.0)
nx, ny   = genv_i("BABALFEM_NX", 200), genv_i("BABALFEM_NY", 200)
feord    = genv_i("BABALFEM_FE_ORDER", 2)
#  p_eta = 0 keeps the historical EQUAL-ORDER spaces (unchanged default).
#  Set BABALFEM_P_ETA = BABALFEM_FE_ORDER-1 for the Taylor-Hood-like pairing, which is
#  the only one measured to reach the theoretical order in BOTH fields. It is
#  NOT automatically the better production choice: at a GIVEN mesh the
#  equal-order spaces were 40x more accurate, because eta sits in a richer
#  space. Compare error-vs-DOF before switching.
p_eta    = genv_i("BABALFEM_P_ETA", 0)
d        = genv_f("BABALFEM_D", 3.5)
Twave    = genv_f("BABALFEM_TWAVE", 1.6)
Awave    = genv_f("BABALFEM_AWAVE", 0.001)
sponge   = genv_f("BABALFEM_SPONGE", 20.0)
mumax    = genv_f("BABALFEM_MUMAX", 10.0)
dt       = genv_f("BABALFEM_DT", 0.02)
periods  = genv_f("BABALFEM_PERIODS", 20.0)
Tfinal   = haskey(ENV, "BABALFEM_TFINAL") ? genv_f("BABALFEM_TFINAL", 0.0) : periods * Twave
save_ev  = genv_i("BABALFEM_SAVE_EVERY", 20)
outdir   = genv("BABALFEM_OUTDIR", joinpath(ROOT, "output", "ring_wave_dist_M$(M)"))

banner("RING WAVE (point source) — distributed", M, (px,py), (nx,ny), nx*ny, outdir)

diags, vert, prob = setup_and_run_distributed(
    cpu_grid=(px,py), M=M, c_bdy=cbdy_override(), p_horizontal=feord, p_eta=p_eta,
    domain=(0.0,L,0.0,L), partition=(nx,ny),
    h_val=d, T_wave=Twave, A_wave=Awave,
    x_wm=L/2, y_wm=L/2,                                   # point source → ring waves
    sponge_wL=sponge, sponge_wR=sponge, sponge_wB=sponge, sponge_wT=sponge, mu_max=mumax,
    T_final=Tfinal, dt=dt,
    regime=regime_sym(), nl_pressure=nl_pressure_sym(), flat_bed=flat_bed_flag(),          # flat sea bed (∇h≡0); set BABALFEM_FLAT_BED=0 for variable bathymetry
    y_wall_bc=:wall, x_wall_bc=false,
    output_dir=outdir, save_every=save_ev,
    write_w=write_w_flag(), write_pressure=write_p_flag(), rho=rho_val(),
    solver_type=solver_sym(), tableau=tableau_sym(),
    nl_iter=nl_iter_val(), nl_tol=nl_tol_val(),
    ls_rtol=ls_rtol_val(), ls_maxiter=ls_maxiter_val(), krylov_m=krylov_m_val(), precond=precond_sym(),
    diag_every=diag_every_val(), diag_csv=diag_csv_flag(),
    div_factor=div_factor_val(), eta_ref=eta_ref_val(),
    print_every=genv_i("BABALFEM_PRINT_EVERY", 10))

is_rank0() && @printf("ring_wave_dist done: %d steps, %d snapshots to %s\n",
                      length(diags), length(diags) ÷ max(save_ev,1), outdir)
