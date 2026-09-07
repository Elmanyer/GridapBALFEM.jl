# ==============================================================
#  run_periodic_plane_wave_dist.jl — DISTRIBUTED plane wave, PERIODIC y edges
#
#  Long-crested plane wave (Gaussian LINE source, y-invariant) in a flume whose
#  lateral (y) boundaries are PERIODIC (y_wall_bc=:periodic) — the numerical
#  infinite-width flume. The top/bottom edges are identified in the mesh, so no
#  lateral walls and no lateral sponges are used; the x-ends keep sponges.
#
#  LAUNCH (px·py MUST equal -n):
#    BALFEM_M=2 BALFEM_PX=8 BALFEM_PY=1 \
#    ~/.julia/bin/mpiexecjl --project=. -n 8 julia --project=. \
#        GridapBALFEM.jl/examples/distributed/run_periodic_plane_wave_dist.jl
#
#  Case-specific env vars (shared ones in _dist_common.jl):
#    BALFEM_LX,BALFEM_LY  domain size [m]         default 400 × 20
#    BALFEM_D           still-water depth [m]   3.5
#    BALFEM_TWAVE       wave period [s]         1.6   (→ kd≈5.5 at d=3.5)
#    BALFEM_AWAVE       wave amplitude [m]      0.001
#    BALFEM_XWM         wavemaker x-position    40
#    BALFEM_SPONGE      x-end sponge width [m]  40
#    BALFEM_MUMAX       sponge strength         5
#    BALFEM_PERIODS     run length in periods   50
# ==============================================================

include(joinpath(@__DIR__, "_dist_common.jl"))

M        = genv_i("BALFEM_M", 2)

#  Vertical BASIS ORDER. The model is named P{p_vert}LFE-{M}: `Pp` is the

#  vertical Lagrange order and `M` the number of vertical elements, so the run

#  says which member of the BALFE-M family it actually exercises. Default p=1

#  reproduces the piecewise-linear models of Yang & Liu.

p_vert  = genv_i("BALFEM_P_VERT", 1)

model_name = "P$(p_vert)LFE-$(M)"
px, py   = genv_i("BALFEM_PX", 32), genv_i("BALFEM_PY", 4)
nx, ny   = genv_i("BALFEM_NX", 2000), genv_i("BALFEM_NY", 100)
feord    = genv_i("BALFEM_FE_ORDER", 2)
#  ⚠ TAYLOR-HOOD IS MANDATORY: p_eta = BALFEM_FE_ORDER-1 (the default here).
#  check_taylor_hood ERRORS on any other pairing. eta plays the pressure role of a
#  Stokes system, so equal order is inf-sup deficient. CLAUDE.md rule 2b.
p_eta    = genv_i("BALFEM_P_ETA", feord - 1)
Lx, Ly   = genv_f("BALFEM_LX", 400.0), genv_f("BALFEM_LY", 20.0)
d        = genv_f("BALFEM_D", 3.5)
Twave    = genv_f("BALFEM_TWAVE", 1.6)
Awave    = genv_f("BALFEM_AWAVE", 0.001)
x_wm     = genv_f("BALFEM_XWM", 40.0)
sponge   = genv_f("BALFEM_SPONGE", 40.0)
mumax    = genv_f("BALFEM_MUMAX", 5.0)
dt       = genv_f("BALFEM_DT", 0.02)
periods  = genv_f("BALFEM_PERIODS", 50.0)
Tfinal   = haskey(ENV, "BALFEM_TFINAL") ? genv_f("BALFEM_TFINAL", 0.0) : periods * Twave
save_ev  = genv_i("BALFEM_SAVE_EVERY", 10)
#  Standardised name — building_files/OUTPUT_NAMING_PROPOSAL.md; generator in
#  src/utilities.jl. ⚠ y_wall_bc must match the setup_and_run_distributed call
#  below: it is what makes the <domain> token truthful.
_name  = output_dir_name(; M=M, p_vert=p_vert, ny=ny, y_wall_bc=:periodic,
                           wave_kind="plane", wave_gen=(:inner),
                           regime=regime_sym(), nl_pressure=nl_pressure_sym(),
                           bed="flat", p_u=feord, p_eta=p_eta,
                           amplitude=Awave, period=Twave, irregular=false)
outdir = haskey(ENV, "BALFEM_OUTDIR") ? genv("BALFEM_OUTDIR", "") :
         unique_output_dir(joinpath(ROOT, "output"), _name)

banner("PERIODIC PLANE WAVE (line source, periodic y) — distributed",
       M, (px,py), (nx,ny), nx*ny, outdir)

diags, vert, prob = setup_and_run_distributed(
    cpu_grid=(px,py), M=M, p_vertical=p_vert, c_bdy=cbdy_override(), p_u=feord, p_eta=p_eta,
    domain=(0.0,Lx,0.0,Ly), partition=(nx,ny),
    h_val=d, T_wave=Twave, A_wave=Awave, x_wm=x_wm, y_wm=nothing,   # line source
    sponge_wL=sponge, sponge_wR=sponge, sponge_wB=0.0, sponge_wT=0.0, mu_max=mumax,
    T_final=Tfinal, dt=dt,
    regime=regime_sym(), nl_pressure=nl_pressure_sym(), flat_bed=flat_bed_flag(),          # flat sea bed (∇h≡0); set BALFEM_FLAT_BED=0 for variable bathymetry
    y_wall_bc=:periodic, x_wall_bc=false,   # ★ periodic lateral (y) boundaries
    output_dir=outdir, save_every=save_ev,
    write_w=write_w_flag(), write_pressure=write_p_flag(), rho=rho_val(),
    solver_type=solver_sym(), tableau=tableau_sym(),
    nl_iter=nl_iter_val(), nl_tol=nl_tol_val(),
    ls_rtol=ls_rtol_val(), ls_maxiter=ls_maxiter_val(), krylov_m=krylov_m_val(), precond=precond_sym(),
    diag_every=diag_every_val(), diag_csv=diag_csv_flag(),
    div_factor=div_factor_val(), eta_ref=eta_ref_val(),
    print_every=genv_i("BALFEM_PRINT_EVERY", 10))

is_rank0() && @printf("periodic_plane_wave_dist done: %d steps, %d snapshots to %s\n",
                      length(diags), length(diags) ÷ max(save_ev,1), outdir)
