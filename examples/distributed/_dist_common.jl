# ==============================================================
#  _dist_common.jl — shared config for the distributed example scripts
#                        (algebraic solver, GridapBALFEM)
#
#  All example scripts read their parameters from environment variables so the
#  SAME script serves any M, any core count, any mesh size — ideal for a
#  cluster job array. Every knob has a sensible default; override via `export`.
#
#  Common environment variables (case-specific ones in each script):
#    BALFEM_M           vertical layers (elements)         default per script
#    BALFEM_PX,BALFEM_PY  MPI process grid (px×py)           MUST satisfy px·py == mpiexec -n
#    BALFEM_NX,BALFEM_NY  horizontal elements                (nx divisible by px, ny by py)
#    BALFEM_FE_ORDER    VELOCITY Q-order (≥2 required)     2
#    BALFEM_P_ETA       SURFACE Q-order; MUST be FE_ORDER-1   FE_ORDER-1
#                     Set to BALFEM_FE_ORDER−1 for the Taylor-Hood-like pairing.
#                     η enters momentum undifferentiated (via ∇·v after IBP), so
#                     equal-order spaces are inf-sup deficient and the analytic
#                     MMS measures order p there, not p+1. But a better RATE is
#                     not a better ANSWER at a given mesh — equal order was 40x
#                     more accurate at nx=24. Compare error-vs-DOF first.
#    BALFEM_DT          time step [s]                      0.02
#    BALFEM_TFINAL      final time [s] (overrides periods) —
#    BALFEM_SAVE_EVERY  VTK snapshot every N steps         (script default)
#    BALFEM_OUTDIR      output directory                   (script default)
#    BALFEM_WRITE_W         write w_s<σ> fields (1/0)       1
#    BALFEM_WRITE_PRESSURE  write p_s<σ> fields (1/0)       1
#    BALFEM_CBDY        comma-sep σ node boundaries        (else optimised M≤4 / uniform)
#    BALFEM_RHO         water density [kg/m³]              1025
#    BALFEM_PRINT_EVERY solver progress report every N steps    (default 10)
#    BALFEM_REGIME      physics regime: linear|nonlinear   nonlinear
#    BALFEM_NL_PRESSURE nonlinear pressure: none|native|full  none
#    BALFEM_FLAT_BED    flat sea bed ∇h≡0 (1/0)  1 flat scripts / 0 bathymetry script
#
#  The algebraic solver runs the FULL physics distributed through the one
#  Gridap path (GMRES+Jacobi+Newton) — no owned V⊗H loop, no linear-core
#  restriction (unlike the old solver's distributed drivers).
#
#  c_bdy: paper-optimised vertical nodes exist for M≤4 (Yang & Liu 2024
#  Table 1); for larger M the driver falls back to a uniform σ-grid.
# ==============================================================

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
using GridapBALFEM
using Printf

genv(k, d)   = get(ENV, k, string(d))
genv_i(k, d) = parse(Int,     genv(k, d))
genv_f(k, d) = parse(Float64, genv(k, d))
genv_b(k, d) = lowercase(genv(k, d)) in ("1", "true", "yes", "on")

# c_bdy override (BALFEM_CBDY="0,0.7,1"), else nothing → driver picks optimised/uniform.
cbdy_override() = haskey(ENV, "BALFEM_CBDY") ?
    parse.(Float64, split(ENV["BALFEM_CBDY"], ",")) : nothing

# shared knobs
write_w_flag()  = genv_b("BALFEM_WRITE_W", 1)
write_p_flag()  = genv_b("BALFEM_WRITE_PRESSURE", 1)
rho_val()       = genv_f("BALFEM_RHO", 1025.0)
regime_sym()         = Symbol(genv("BALFEM_REGIME", "nonlinear"))
nl_pressure_sym()    = Symbol(genv("BALFEM_NL_PRESSURE", "none"))
# Sea-bed geometry: false = variable bathymetry (∇h≠0), true = flat bed (∇h≡0).
# Default per script (BALFEM_FLAT_BED): flat-bed cases pass 1, the bathymetry case passes 0.
flat_bed_flag(default::Int=1) = genv_b("BALFEM_FLAT_BED", default)

# ---- Time integrator + nonlinear/linear solver controls ---------------------
#  Default integrator is the fully-implicit RungeKutta :SDIRK_2_2 (L-stable,
#  more robust than Crank–Nicolson). Override caps/tolerances per cluster job:
#    BALFEM_SOLVER      integrator: sdirk|theta        sdirk
#    BALFEM_TABLEAU     RK tableau (solver=sdirk)       SDIRK_2_2
#    BALFEM_NL_ITER     max Newton iterations / stage    50
#    BALFEM_NL_TOL      Newton tolerance (production)    1e-5
#    BALFEM_LS_MAXITER  max GMRES iterations / Newton    1000   (TIME bound)
#    BALFEM_LS_RTOL     GMRES relative tolerance         1e-5  (measured: 1e-9→1e-6→1e-5
#                       each cut GMRES iterations with max|η| unmoved. NOTE this now
#                       EQUALS BALFEM_NL_TOL rather than sitting one order below it —
#                       see the ladder discussion in markdown_files/CONFIGURATION.md §4)
#    BALFEM_KRYLOV_M    GMRES basis size, restart=true    100   (MEMORY bound)
#      NOTE BALFEM_KRYLOV_M and BALFEM_LS_MAXITER are different things: the basis
#      costs m+1 distributed vectors plus a dense (m+1)×m Hessenberg PER RANK at
#      every numerical setup, so it drives memory; BALFEM_LS_MAXITER only decides
#      how long GMRES may iterate. Conflating them (passing the iteration cap as
#      the basis size) is what OOM-killed the 32-rank runs.
# ---- Field diagnostics (monitor.jl) -----------------------------------------
#    BALFEM_DIAG_EVERY  sample every N steps            0 → follow BALFEM_PRINT_EVERY
#                                                     −1 → disable the whole block
#    BALFEM_DIAG_CSV    write <outdir>/diagnostics.csv  1
#    BALFEM_DIV_FACTOR  abort at div_factor·eta_ref     20
#    BALFEM_ETA_REF     override the amplitude scale    (unset → inferred from the
#                                                      forcing: A_wave / Hs / peak η₀)
#  The diagnostics are ON by default and cost within noise (measured −0.8 %/+2.2 %).
#  They are what makes a failed cluster run diagnosable from its log alone: WHERE
#  max|η| sits, its interior/damped split, |u|/|η|, mass & energy, GMRES
#  saturation, and per-rank RSS. Read one with examples/inspect_run.jl.
diag_every_val() = genv_i("BALFEM_DIAG_EVERY", 0)
diag_csv_flag()  = genv_b("BALFEM_DIAG_CSV", 1)
div_factor_val() = genv_f("BALFEM_DIV_FACTOR", 20.0)
eta_ref_val()    = haskey(ENV, "BALFEM_ETA_REF") ? genv_f("BALFEM_ETA_REF", 0.0) : nothing

solver_sym()     = Symbol(genv("BALFEM_SOLVER", "sdirk"))
tableau_sym()    = Symbol(genv("BALFEM_TABLEAU", "SDIRK_2_2"))
nl_iter_val()    = genv_i("BALFEM_NL_ITER", 50)
nl_tol_val()     = genv_f("BALFEM_NL_TOL", 1e-5)
ls_maxiter_val() = genv_i("BALFEM_LS_MAXITER", 1000)
ls_rtol_val()    = genv_f("BALFEM_LS_RTOL", 1e-5)
krylov_m_val()   = genv_i("BALFEM_KRYLOV_M", 100)
#    BALFEM_PRECOND   GMRES preconditioner: jacobi|schwarz|gs   jacobi
#      Jacobi is weak for this operator (450-760 iters/solve measured); :schwarz
#      does an exact LU on each rank's own block. See build_preconditioner.
precond_sym()    = Symbol(genv("BALFEM_PRECOND", "jacobi"))

#  Extra output-name tokens (the trailing `[_<extra>...]` of the grammar in
#  markdown_files/OUTPUT_NAMING_PROPOSAL.md). Use it whenever a launcher varies
#  something the grammar has no field for -- dt, tolerances, rank count -- so the
#  run does NOT land on its sibling's name and get silently suffixed _v2.
#  Comma-separated, e.g.  BALFEM_NAME_EXTRA=dt0.01
name_extra() = String[strip(t) for t in split(genv("BALFEM_NAME_EXTRA", ""), ",") if !isempty(strip(t))]

# rank-0 detection BEFORE MPI.Init (OpenMPI / MPICH set these per rank).
is_rank0() = get(ENV, "OMPI_COMM_WORLD_RANK",
                 get(ENV, "PMI_RANK", "0")) == "0"
my_rank()  = parse(Int, get(ENV, "OMPI_COMM_WORLD_RANK",
                            get(ENV, "PMI_RANK", "0")))

# ---- Dirichlet boundary wave generation (WaveSpec.jl sea states) -------------
#  Environment variables (used by run_irregular_sea_dist / run_directional_sea_dist):
#    BALFEM_HS          significant wave height [m]        0.002 (linear regime)
#    BALFEM_TP          peak period [s]                    1.6
#    BALFEM_GAMMA       JONSWAP peakedness (0 = estimate)  3.3
#    BALFEM_NFREQ       frequency samples (bins = nf-1)    21
#    BALFEM_FMIN_FAC    fmin = 1/(FMIN_FAC·Tp)             2.5
#    BALFEM_FMAX_FAC    fmax = 1/(FMAX_FAC·Tp)             0.75  (keep kd ≲ kd_app AND
#                                                        ≥6 cells per shortest λ — see README)
#    BALFEM_SAMPLING    bin spacing: uniform|energy        uniform (leakage-free analysis)
#    BALFEM_NTHETA      angle samples (≤1 → long-crested)  0
#    BALFEM_SPREAD_STD  cosine-power spreading σθ [deg]    20
#    BALFEM_THETA_MAX   angular truncation ±θmax [deg]     60
#    BALFEM_SEED        phase seed (reproducible)          20260723
#    BALFEM_BC_SIDE     generation boundary left|right     left
#    BALFEM_BC_PROFILE  vertical polarization model|airy   model
#    BALFEM_TRAMP       Hann ramp [s] (unset → 2·Tp)       —
#    BALFEM_RELAX       relaxation zone at the inflow      0
#    BALFEM_RELAX_W     zone width [m] (0 → one peak λ)    0
#  The WaveInput conversion snapshots the SEEDED phases into plain arrays, so
#  every MPI rank builds an identical component table — no communication needed.
hs_val()         = genv_f("BALFEM_HS", 0.002)
tp_val()         = genv_f("BALFEM_TP", 1.6)
bc_side_sym()    = Symbol(genv("BALFEM_BC_SIDE", "left"))
bc_profile_sym() = Symbol(genv("BALFEM_BC_PROFILE", "model"))
tramp_val()      = haskey(ENV, "BALFEM_TRAMP") ? genv_f("BALFEM_TRAMP", 0.0) : nothing
relax_flag()     = genv_b("BALFEM_RELAX", 0)
relax_w_val()    = genv_f("BALFEM_RELAX_W", 0.0)

"""
    build_airy_state(h_val; directional=false) → WaveSpec.AiryWaves.AiryState

Env-configured stochastic sea state: JONSWAP spectrum, uniform-frequency (or
equal-energy) bins, optional cosine-power angular spreading, seeded phases.
`directional=true` switches the default spreading on (`BALFEM_NTHETA` ≥ 2).

⚠ **TWO SEEDS HAVE TO BE PINNED, NOT ONE, OR EVERY MPI RANK BUILDS A DIFFERENT
SEA.** `AiryState` seeds the *phases*; `DiscreteAngularSpreading` separately
seeds the *propagation directions* — `get_angles` is `sort(rand(rng(sm.seed),
sm.distribution, sm.nθ))` (`WaveSpec/src/AngularSpreading/AngularSpreading.jl`
:157), and both of its constructors default that seed to `abs(rand(Int64))`
drawn from the PROCESS-GLOBAL RNG, which Julia seeds from system entropy per
process. Under `mpiexec -n N` every rank runs this function independently, so
each one drew its own θⱼ *and* its own Δθⱼ (hence its own component amplitudes,
`A ∝ √(D·Δθ)`).

Measured consequence, before this fix (2026-09-11, `output/snellius/
dist_small_new`): on the 42-rank directional run, η on the Dirichlet inflow edge
was smooth inside each rank block and jumped at exactly the six `cpu_grid=(7,6)`
y-cuts — y = 3.5, 6.75, 10.0, 13.25, 16.5 — flipping sign at 16.5
(+2.36e-02 → −2.61e-02). The domain carried six uncorrelated seas side by side.
The long-crested path is hit too, weakly but not harmlessly: its
`DiscreteAngularSpreading(θ::Real)` draws two angles from `Uniform(θ±0.001)`, so
the per-rank central angle moved ±0.04° and Δθ by 3×, which is the 1e-3
transverse asymmetry those runs carried from step 1 — the seed the transverse
velocity instability then amplified to 0.89 before the run died.

The spread must be re-seeded BEFORE the `AiryState` is built: its constructor
calls `get_central_angles(spread)` to fill `state.θ`, and `get_amplitudes` later
re-reads `get_weights`/`get_bandwidths` from the same seed.
"""
function build_airy_state(h_val; directional::Bool=false)
    Hs, Tp = hs_val(), tp_val()
    γ      = genv_f("BALFEM_GAMMA", 3.3)
    nf     = genv_i("BALFEM_NFREQ", 21)
    fmin   = 1.0/(genv_f("BALFEM_FMIN_FAC", 2.5)*Tp)
    fmax   = 1.0/(genv_f("BALFEM_FMAX_FAC", 0.75)*Tp)
    dom    = lowercase(genv("BALFEM_SAMPLING", "uniform")) == "energy" ?
             WaveSpec.SpectralSampling.Energy : WaveSpec.SpectralSampling.Frequency
    spec   = γ > 0 ? WaveSpec.ContinuousSpectrums.JONSWAP(Hs, Tp, γ) :
                     WaveSpec.ContinuousSpectrums.JONSWAP(Hs, Tp)
    dspec  = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
                 spec, WaveSpec.SpectralSampling.UniformSampling(),
                 fmin, fmax, nf; domain=dom, mess=is_rank0())
    seed = genv_i("BALFEM_SEED", 20260723)
    nθ = genv_i("BALFEM_NTHETA", directional ? 7 : 0)
    spread = nθ <= 1 ? WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0; seed=seed) :
             WaveSpec.AngularSpreading.DiscreteAngularSpreading(
                 :cosinepow, 0.0, genv_f("BALFEM_SPREAD_STD", 20.0)*pi/180,
                 -genv_f("BALFEM_THETA_MAX", 60.0)*pi/180,
                  genv_f("BALFEM_THETA_MAX", 60.0)*pi/180, nθ)
    #  ⚠ LOAD-BEARING, NOT TIDYING (see the docstring): the :cosinepow branch has
    #  no `seed` keyword, so its directions come from the process-global RNG until
    #  this line replaces them. Re-seed BEFORE AiryState is constructed — that
    #  constructor reads get_central_angles(spread) to fill state.θ.
    spread = WaveSpec.AngularSpreading.change_seed!(spread, seed)
    state = WaveSpec.AiryWaves.AiryState(dspec, spread, h_val)
    state = WaveSpec.AiryWaves.change_seed!(state, seed)
    #  Rank-invariance is an INVARIANT of this function, so assert it rather than
    #  trusting it: a silent regression here is invisible in a serial test and
    #  cost a 42-rank cluster run last time. θ is what the partition corrupted.
    @assert issorted(state.θ) "spread angles unsorted — WaveSpec contract changed"
    #  ONE LINE PER RANK, deliberately. The whole failure mode was that a
    #  per-rank difference is invisible from rank 0, so a rank-0-only print would
    #  reproduce exactly the blind spot. Check a finished job with:
    #      grep seastate job.out | awk '{print $NF}' | sort -u | wc -l    # must be 1
    @printf("[seastate] rank %-4d nω=%-3d nθ=%-3d θ(deg)=[%s] fp=%s\n",
            my_rank(), state.nω, state.nθ,
            join((@sprintf("%.3f", rad2deg(t)) for t in state.θ), ","),
            airy_state_fingerprint(state))
    flush(stdout)
    return state
end

"""
    airy_state_fingerprint(state) → String

Rank-comparable digest of the sea state actually built. Print it from every rank
(not just rank 0) when a `:bc_gen` run starts: identical strings across ranks is
the cheap proof that `build_airy_state` is partition-independent. Divergent
strings mean the two-seed bug is back, and the run is generating one sea per
rank — see `build_airy_state`.
"""
function airy_state_fingerprint(state)
    #  Hash EXACTLY what the boundary evaluation consumes, not a proxy for it.
    #  `WaveInput(vert, ::AiryState)` (src/waveinput.jl:192) reads precisely
    #  get_amplitudes, get_random_phases, state.ω and state.θ, and nothing else
    #  in src/ draws a random number. So these four arrays being equal across
    #  ranks IS the statement that every rank evaluates the same η(y,t), u(y,t).
    #  (state.k is deliberately absent: WaveInput DISCARDS it and re-solves the
    #  wavenumbers with the solver's own g and the model dispersion.)
    h = hash(round.(WaveSpec.AiryWaves.get_amplitudes(state),     digits=12))
    h = hash(round.(WaveSpec.AiryWaves.get_random_phases(state),  digits=12), h)
    h = hash(round.(state.ω, digits=12), h)
    h = hash(round.(state.θ, digits=12), h)
    return string(h, base=16, pad=16)
end

function banner(title, M, cpu_grid, partition, ncells, outdir)
    is_rank0() || return
    @printf("############################################################\n")
    @printf("# %s  [ALGEBRAIC solver, stacked layout]\n", title)
    @printf("#   M=%d layers | cpu_grid=%s (%d ranks) | mesh=%s = %d cells\n",
            M, string(cpu_grid), prod(cpu_grid), string(partition), ncells)
    @printf("#   regime=%s nl_pressure=%s flat_bed=%s\n",
            string(regime_sym()), string(nl_pressure_sym()), string(flat_bed_flag()))
    @printf("#   write_w=%s write_pressure=%s | out=%s\n",
            string(write_w_flag()), string(write_p_flag()), outdir)
    @printf("############################################################\n")
    flush(stdout)
end
