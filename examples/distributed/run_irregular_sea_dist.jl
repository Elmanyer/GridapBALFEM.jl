# ==============================================================
#  run_irregular_sea_dist.jl — DISTRIBUTED long-crested JONSWAP sea
#                              via Dirichlet boundary generation (PRODUCTION)
#
#  A stochastic sea state (WaveSpec.jl: JONSWAP spectrum, seeded phases —
#  identical component table on every rank) enters through time-varying
#  Dirichlet data (η, 𝖴x) on the LEFT boundary — NO interior wavemaker — and
#  is absorbed by the right sponge. The :model discrete-eigenmode polarization
#  makes the boundary data an exact linearised discrete solution.
#
#  LAUNCH (px·py MUST equal -n):
#    BABALFEM_M=2 BABALFEM_PX=32 BABALFEM_PY=4 \
#    ~/.julia/bin/mpiexecjl --project=. -n 128 julia --project=. \
#        GridapBABALFEM.jl/examples/distributed/run_irregular_sea_dist.jl
#  (SLURM: see run/run_snellius.sh — point it at this script.)
#
#  Case-specific env vars (sea-state + shared ones in _dist_common.jl):
#    BABALFEM_LX,BABALFEM_LY  domain size [m]           default 400 × 20
#    BABALFEM_D           still-water depth [m]     3.5   (kd_p≈3.6 at Tp=1.6)
#    BABALFEM_SPONGE      right-sponge width [m]    40    (≥ 4 peak wavelengths)
#    BABALFEM_MUMAX       sponge strength           5
#    BABALFEM_PERIODS     run length in Tp units (if BABALFEM_TFINAL unset)  50
#
#  Guidance (see README): keep the spectral band inside the model — fmax such
#  that kd(fmax) ≲ kd_app(M) — and resolve the shortest component with ≥6
#  cells/λ; the defaults (LFE-2, band [1/(2.5Tp), 1/(0.75Tp)], 2000×100 mesh)
#  satisfy both. Field output (w_s/p_s) defaults OFF here (large runs).
# ==============================================================

include(joinpath(@__DIR__, "_dist_common.jl"))

M        = genv_i("BABALFEM_M", 3)
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
sponge   = genv_f("BABALFEM_SPONGE", 40.0)
mumax    = genv_f("BABALFEM_MUMAX", 5.0)
dt       = genv_f("BABALFEM_DT", 0.02)
periods  = genv_f("BABALFEM_PERIODS", 20.0)
Tp       = tp_val()
Tfinal   = haskey(ENV, "BABALFEM_TFINAL") ? genv_f("BABALFEM_TFINAL", 0.0) : periods * Tp
save_ev  = genv_i("BABALFEM_SAVE_EVERY", 1)     # SAVE A LOT OF SCREENSHOTS FOR VISUALISATION
outdir   = genv("BABALFEM_OUTDIR", joinpath(ROOT, "output", "irregular_sea_dist_M$(M)"))

state = build_airy_state(d; directional=false)

banner("IRREGULAR SEA (JONSWAP, Dirichlet BC generation) — distributed",
       M, (px,py), (nx,ny), nx*ny, outdir)
is_rank0() && @printf("#   Hs=%.4g m Tp=%.3g s | bc_side=%s profile=%s seed=%d\n",
                      hs_val(), Tp, string(bc_side_sym()),
                      string(bc_profile_sym()), genv_i("BABALFEM_SEED", 20260723))

diags, vert, prob = setup_and_run_distributed(
    cpu_grid=(px,py), M=M, c_bdy=cbdy_override(), p_horizontal=feord, p_eta=p_eta,
    domain=(0.0,Lx,0.0,Ly), partition=(nx,ny),
    h_val=d, T_wave=Tp, A_wave=hs_val()/2,
    wave_bc=state, bc_side=bc_side_sym(), bc_profile=bc_profile_sym(),
    T_ramp=tramp_val(), relax_bc=relax_flag(), relax_width=relax_w_val(),
    sponge_wL=0.0, sponge_wR=sponge, sponge_wB=0.0, sponge_wT=0.0, mu_max=mumax,
    T_final=Tfinal, dt=dt,
    regime=regime_sym(), nl_pressure=nl_pressure_sym(), flat_bed=flat_bed_flag(),          # flat sea bed (∇h≡0); set BABALFEM_FLAT_BED=0 for variable bathymetry
    y_wall_bc=:wall, x_wall_bc=false,
    output_dir=outdir, save_every=save_ev,
    write_w=genv_b("BABALFEM_WRITE_W", 0), write_pressure=genv_b("BABALFEM_WRITE_PRESSURE", 0),
    rho=rho_val(),
    solver_type=solver_sym(), tableau=tableau_sym(),
    nl_iter=nl_iter_val(), nl_tol=nl_tol_val(),
    ls_rtol=ls_rtol_val(), ls_maxiter=ls_maxiter_val(), krylov_m=krylov_m_val(), precond=precond_sym(),
    diag_every=diag_every_val(), diag_csv=diag_csv_flag(),
    div_factor=div_factor_val(), eta_ref=eta_ref_val(),
    print_every=genv_i("BABALFEM_PRINT_EVERY", 10))

is_rank0() && @printf("irregular_sea_dist done: %d steps, %d snapshots to %s\n",
                      length(diags), save_ev > 0 ? length(diags) ÷ save_ev : 0, outdir)
