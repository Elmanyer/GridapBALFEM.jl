# ==============================================================
#  run_plane_wave_dist.jl — DISTRIBUTED long-crested plane wave (algebraic)
#
#  Gaussian LINE-source wavemaker → plane waves down a long flume, absorbed by
#  sponge layers at both x-ends. Writes η, per-node u/v components AND the
#  reconstructed vertical velocity (w_s*) and total pressure (p_s*) fields.
#  Runs the FULL nonlinear physics distributed (one Gridap path).
#
#  LAUNCH (px·py MUST equal -n):
#    BABALFEM_M=2 BABALFEM_PX=8 BABALFEM_PY=1 \
#    ~/.julia/bin/mpiexecjl --project=. -n 8 julia --project=. \
#        GridapBABALFEM.jl/examples/distributed/run_plane_wave_dist.jl
#
#  Case-specific env vars (shared ones in _dist_common.jl):
#    BABALFEM_LX,BABALFEM_LY  domain size [m]         default 400 × 20
#    BABALFEM_D           still-water depth [m]   3.5
#    BABALFEM_TWAVE       wave period [s]         1.6   (→ kd≈5.5 at d=3.5)
#    BABALFEM_AWAVE       wave amplitude [m]      0.001
#    BABALFEM_XWM         wavemaker x-position    40
#    BABALFEM_SPONGE      sponge width [m]        40
#    BABALFEM_MUMAX       sponge strength         5
#    BABALFEM_PERIODS     run length in wave periods (if BABALFEM_TFINAL unset)  50
# ==============================================================

include(joinpath(@__DIR__, "_dist_common.jl"))

M        = genv_i("BABALFEM_M", 2)
px, py   = genv_i("BABALFEM_PX", 32), genv_i("BABALFEM_PY", 4)
nx, ny   = genv_i("BABALFEM_NX", 2000), genv_i("BABALFEM_NY", 100)
feord    = genv_i("BABALFEM_FE_ORDER", 2)
#  p_eta = 0 keeps the historical EQUAL-ORDER spaces (unchanged default).
#  Set BABALFEM_P_ETA = BABALFEM_FE_ORDER-1 for the Taylor-Hood-like pairing, which is
#  the only one measured to reach the theoretical order in BOTH fields. It is
#  NOT automatically the better production choice: at a GIVEN mesh the
#  equal-order spaces were 40x more accurate, because eta sits in a richer
#  space. Compare error-vs-DOF before switching.
p_eta    = genv_i("BABALFEM_P_ETA", 0)
Lx, Ly   = genv_f("BABALFEM_LX", 400.0), genv_f("BABALFEM_LY", 20.0)
d        = genv_f("BABALFEM_D", 3.5)
Twave    = genv_f("BABALFEM_TWAVE", 1.6)
Awave    = genv_f("BABALFEM_AWAVE", 0.001)
x_wm     = genv_f("BABALFEM_XWM", 40.0)
sponge   = genv_f("BABALFEM_SPONGE", 40.0)
mumax    = genv_f("BABALFEM_MUMAX", 5.0)
dt       = genv_f("BABALFEM_DT", 0.02)
periods  = genv_f("BABALFEM_PERIODS", 50.0)
Tfinal   = haskey(ENV, "BABALFEM_TFINAL") ? genv_f("BABALFEM_TFINAL", 0.0) : periods * Twave
save_ev  = genv_i("BABALFEM_SAVE_EVERY", 10)
outdir   = genv("BABALFEM_OUTDIR", joinpath(ROOT, "output", "plane_wave_dist_M$(M)"))

banner("PLANE WAVE (line source) — distributed", M, (px,py), (nx,ny), nx*ny, outdir)

diags, vert, prob = setup_and_run_distributed(
    cpu_grid=(px,py), M=M, c_bdy=cbdy_override(), p_horizontal=feord, p_eta=p_eta,
    domain=(0.0,Lx,0.0,Ly), partition=(nx,ny),
    h_val=d, T_wave=Twave, A_wave=Awave, x_wm=x_wm, y_wm=nothing,   # line source
    sponge_wL=sponge, sponge_wR=sponge, sponge_wB=0.0, sponge_wT=0.0, mu_max=mumax,
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

is_rank0() && @printf("plane_wave_dist done: %d steps, %d snapshots to %s\n",
                      length(diags), length(diags) ÷ max(save_ev,1), outdir)
