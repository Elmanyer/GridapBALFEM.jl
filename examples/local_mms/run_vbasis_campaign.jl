# ==============================================================
#  run_vbasis_campaign.jl — THE VERTICAL-BASIS CONVERGENCE CAMPAIGN
#
#  Plan (read it first): building_files/MMS_VBASIS_CAMPAIGN.md
#
#  Phase 1  optimise the σ-mesh for each (M,p) by maximising the applicable kd,
#           exactly as Yang & Liu (2024) did — and CALIBRATE the optimiser against
#           their published p=1 table before trusting anything it says about p≥2.
#  Phase 2  the MMS convergence matrix: all 8 model configurations × every basis,
#           SPATIAL (fix dt, refine h) and TEMPORAL (fix h, refine dt).
#
#  ⚠ THE TWO PHASES ANSWER DIFFERENT QUESTIONS AND MUST NOT BE CONFLATED. The
#  node positions change the error CONSTANT, never the ORDER: Phase 1 is a
#  DISPERSION-accuracy result, Phase 2 a DISCRETISATION-order result. Neither is
#  evidence for the other.
#
#  Output: ONE file, output/local/mms_campaign/campaign_results.csv, one row per
#  mesh (phase=vmesh) or per refinement level (phase=space|time), plus a rendered
#  campaign_report.md. Workers checkpoint into their own files as they go, so an
#  interrupted campaign keeps everything already measured.
#
#  ENV
#    VBC_WORKERS   worker processes                     8
#    VBC_PHASE     1 | 2 | both                          both
#    VBC_LEVELS_S  spatial refinement levels             4
#    VBC_LEVELS_T  temporal refinement levels            4
#    VBC_NX0       coarsest nx (spatial)                 8
#    VBC_NXT       fixed nx (temporal)                   36
#    VBC_NSTEPS_S  steps per STATIC spatial level        20   (was 100; 15x, measured)
#    VBC_MODELS    comma list of model numbers 1-8       1,2,3,4,5,6,7,8
#    VBC_OUT       output directory                      output/local/mms_campaign
#
#  RUN:  julia --project=. examples/local_mms/run_vbasis_campaign.jl
# ==============================================================

using Distributed, Printf, Dates, LinearAlgebra

gs(k,d)=get(ENV,k,d); gi(k,d)=parse(Int,get(ENV,k,string(d)))

const NWORK   = gi("VBC_WORKERS", 8)
const PHASE   = gs("VBC_PHASE", "both")
const LEVELS_S= gi("VBC_LEVELS_S", 4)
const LEVELS_T= gi("VBC_LEVELS_T", 4)
const NX0     = gi("VBC_NX0", 8)
const NXT     = gi("VBC_NXT", 36)
const NSTEPS_S= gi("VBC_NSTEPS_S", 20)   # static-study steps; see the note in run_one
const OUT     = gs("VBC_OUT", "output/local/mms_campaign")
const MODELSEL= parse.(Int, split(gs("VBC_MODELS","1,2,3,4,5,6,7,8"), ","))
mkpath(OUT)
const CKPT = joinpath(OUT, "checkpoints"); mkpath(CKPT)
#  Checkpoint files are APPENDED to, so a stale one from an earlier campaign would
#  be merged into this campaign's results file and silently misattributed. Clear
#  them at the start, not at the end: a crashed run must leave its evidence behind.
for f in readdir(CKPT); endswith(f,".csv") && rm(joinpath(CKPT,f)); end

#  Physical case. d ≠ 1 IS DELIBERATE: d = 1 makes multiplication by h the
#  identity and the whole h-weighting of the momentum equation unobservable.
const LX, LY, DEPTH, A_BED, KBX = 1.7, 1.1, 2.5, 0.2, 1.3
const G_ACC = 9.81
const P_U, P_ETA = 3, 2          # Q3/Q2 — equal order is inf-sup deficient
const NL_TOL = 1e-12             # 8 orders tighter than production, NOT negotiable

#  The five bases of the plan: matched-Nσ pairs, so each comparison isolates basis
#  SHAPE from basis SIZE.  (M=1 has no free interior boundary ⇒ nothing to optimise.)
const BASES = [(M=2,p=1), (M=1,p=2), (M=3,p=1), (M=4,p=1), (M=2,p=2)]

#  The eight model configurations. rate_u=false exactly where the MMS cannot reach
#  the operator: :full's velocity error is capped by frozen L² projections, not by
#  the mesh, so a rate gate there would assert something false.
const MODELS = Dict(
    1=>(regime=:linear,   flat_bed=true,  nlp=:none,   rate_u=true),
    2=>(regime=:linear,   flat_bed=false, nlp=:none,   rate_u=true),
    3=>(regime=:nonlinear,flat_bed=true,  nlp=:none,   rate_u=true),
    4=>(regime=:nonlinear,flat_bed=false, nlp=:none,   rate_u=true),
    5=>(regime=:nonlinear,flat_bed=true,  nlp=:native, rate_u=true),
    6=>(regime=:nonlinear,flat_bed=false, nlp=:native, rate_u=true),
    7=>(regime=:nonlinear,flat_bed=true,  nlp=:full,   rate_u=false),
    8=>(regime=:nonlinear,flat_bed=false, nlp=:full,   rate_u=false),
)

#  Yang & Liu (2024) Table 1 — the CALIBRATION STANDARD for the Phase-1 optimiser.
#  `props` = (C, C_g, γ) applicable ranges from StokesWaveFourierAnalysis.tex
#  Table 4.1 (the BALFE-M column, i.e. OUR functionals at THEIR nodes — the right
#  standard to calibrate against, since the LFE-M column carries their tolerance
#  conventions, which differ on C_g and γ).
const YL_TABLE = Dict(
    2 => (c=[0.0,0.728,1.0],             kd=10.9,  props=(10.84,  6.42,  5.48)),
    3 => (c=[0.0,0.726,0.925,1.0],       kd=39.2,  props=(39.23, 24.49, 21.87)),
    4 => (c=[0.0,0.745,0.923,0.977,1.0], kd=127.9, props=(127.92,76.20, 68.76)),
)

#  ---- workers up FIRST, so master and workers share ONE definition of every
#  ---- shared name. Defining the row builder twice is how a column set drifts
#  ---- between the two writers of the same file.
addprocs(NWORK)
@everywhere using GridapBALFEM, Printf, LinearAlgebra

@everywhere begin
    const COLS = [:phase,:basis,:M,:p_vert,:Nsigma,:c_bdy,:kd_app,:kd_per_prop,:model,:regime,
                  :flat_bed,:nl_pressure,:integrator,:p_u,:p_eta,:level,:h,:dt,:ndofs,:e_eta,
                  :e_u,:pw_eta,:pw_u,:fit_eta,:fit_u,:opt_eta,:opt_u,:rate_gated_u,:verdict,:note]

    """
        csv_row(; kwargs...) → String

    One results row, built BY NAME from `COLS`, with anything unsupplied left empty.

    Deliberately not a `@sprintf` format string: the rows here have 30 fields and
    three different shapes (mesh / study / error), and a positional format string for
    that is a silent-corruption machine — miscount the commas by one and every column
    after it shifts, producing a file that parses cleanly and means something else.
    Naming the fields makes the miscount unrepresentable, and an unknown name is an
    error rather than a dropped value. Commas are stripped; nothing here is quoted.
    """
    csv_row(; kwargs...) = begin
        d = Dict(kwargs)
        for k in keys(d); k in COLS || error("csv_row: unknown column :$k"); end
        join((replace(string(get(d, c, "")), "," => " ") for c in COLS), ",")
    end

    cbstr(c) = "[" * join(map(x -> @sprintf("%.4f", x), c), " ") * "]"

    const _LX, _LY, _DEPTH, _AB, _KBX = 1.7, 1.1, 2.5, 0.2, 1.3
    const _PU, _PETA = 3, 2
    #  ⚠ THE NEWTON TOLERANCE IS PER-REGIME, AND 1e-12 IS UNREACHABLE FOR THE
    #  NONLINEAR MODELS. The nonlinear Jacobians are quasi-Newton BY DESIGN (the
    #  pressure blocks' η-dependence is frozen), so Newton converges LINEARLY and
    #  stalls around 1e-10 — `test_mms_convergence_nonlinear.jl` pins 1e-9 for
    #  exactly this reason and says so in its header. The campaign's first attempt
    #  specified a flat 1e-12 on the "8 orders tighter than production" rule and
    #  every nonlinear study would have burned its full 400-iteration budget at
    #  every step, forever. The linear models DO reach 1e-12 in one step.
    #
    #  1e-9 cannot contaminate the rates here: the FINEST spatial error in this
    #  campaign is ~4e-6 (Q3/Q2 at nx=64), three orders above it.
    #  ⚠ NEVER loosen this further to make a study pass — the algebraic error would
    #  then sit inside the discretisation error being measured. Raise nl_iter, or
    #  fix the window, instead.
    nltol_for(regime) = regime === :linear ? 1e-12 : 1e-9
end
const HDR = join(COLS, ",")

println("#"^92)
println("#  VERTICAL-BASIS CONVERGENCE CAMPAIGN   started $(Dates.now())")
println("#    plan: building_files/MMS_VBASIS_CAMPAIGN.md")
println("#    bases: ", join(["P$(b.p)LFE-$(b.M)(Nσ=$(b.M*b.p+1))" for b in BASES], " "))
println("#    models: $MODELSEL   workers: $NWORK   out: $OUT")
println("#"^92); flush(stdout)

# ===========================================================================
#  PHASE 1 — optimised σ-meshes
# ===========================================================================
#  Run on the MASTER only: it is pure vertical algebra (no PDE), it is cheap, and
#  Phase 2 needs its output before it can start. (Workers are already up and warm;
#  they idle for the few minutes this takes, which is cheaper than starting them
#  after Phase 1 and paying `using GridapBALFEM` x8 on the critical path.)

"""
    cbdy_from_z(z, M) → Vector{Float64}

Unconstrained ⇒ constrained. The M element widths are `softmax([0; z])`, so they
are positive and sum to 1 BY CONSTRUCTION and the ordering `0 < c_1 < … < 1` can
never be violated — no penalty terms, no repaired iterates, no optimiser wandering
into a mesh that `assemble_vertical_tensors` would reject.
"""
function cbdy_from_z(z::Vector{Float64}, M::Int)
    w = exp.([0.0; z]); w ./= sum(w)
    c = [0.0; cumsum(w)]
    c[end] = 1.0                       # kill accumulated round-off at the endpoint
    return c
end

"""
    nelder_mead(f, z0; ...) → (z_best, f_best)

Plain Nelder–Mead. Local, in-file, and deliberately so: the alternative is a new
dependency for forty lines of standard simplex bookkeeping.
"""
function nelder_mead(f, z0::Vector{Float64}; step=0.35, maxiter=400, ftol=1e-10)
    n = length(z0)
    n == 0 && return (z0, f(z0))
    simplex = [copy(z0)]
    for i in 1:n
        z = copy(z0); z[i] += step; push!(simplex, z)
    end
    fv = [f(z) for z in simplex]
    for _ in 1:maxiter
        o = sortperm(fv); simplex = simplex[o]; fv = fv[o]
        (abs(fv[end] - fv[1]) < ftol * (abs(fv[1]) + ftol)) && break
        centroid = sum(simplex[1:end-1]) / n
        zr = centroid + (centroid - simplex[end]); fr = f(zr)
        if fr < fv[1]
            ze = centroid + 2.0*(centroid - simplex[end]); fe = f(ze)
            fe < fr ? (simplex[end]=ze; fv[end]=fe) : (simplex[end]=zr; fv[end]=fr)
        elseif fr < fv[end-1]
            simplex[end]=zr; fv[end]=fr
        else
            zc = centroid + 0.5*(simplex[end] - centroid); fc = f(zc)
            if fc < fv[end]
                simplex[end]=zc; fv[end]=fc
            else
                for i in 2:n+1
                    simplex[i] = simplex[1] + 0.5*(simplex[i]-simplex[1]); fv[i]=f(simplex[i])
                end
            end
        end
    end
    o = argmin(fv)
    return (simplex[o], fv[o])
end

#  Tolerances of eq: applicable range definition — RELATIVE 2 % on C and C_g,
#  ABSOLUTE 0.02 on γ (which changes sign near kd ≈ 1.2, where a relative measure
#  is singular). Dividing each error by its own tolerance puts all three on one
#  scale, so "≤ 1" is the feasibility statement for the whole property set.
const TOLS = (C=0.02, Cg=0.02, gamma=0.02)
enorm(v, kd) = maximum(property_errors(v, kd) ./ (TOLS.C, TOLS.Cg, TOLS.gamma))

"""
    minimax_error(v, K; n=120) → max_{0<kd≤K} normalised error

The inner objective's integrand: the WORST normalised error anywhere in `[0,K]`,
over all three properties.

⚠ **NOTE THE QUANTIFIER.** The applicable range requires the tolerance to hold
THROUGHOUT `[0,K]`, not merely AT `K` — the error curves are not monotone, and a
mesh can leave the band and re-enter. Sampling geometrically puts points where the
structure is: γ's difficulty is at small `kd`, C's at large.
"""
minimax_error(v, K::Float64; n::Int=400) =
    maximum(enorm(v, kd) for kd in exp.(range(log(1e-3), log(K); length=n)))

"""
    applicable_range_multi(v; kd_max=400.0, n=4000, refine=1e-4) → Float64

The MULTI-PROPERTY applicable range: `max{K : enorm(kd) <= 1 for all kd <= K}`,
i.e. the smallest of the three per-property ranges, from ONE fine geometric scan
with the crossing bisected.

⚠ **This, not the inner objective's sample grid, is what may be reported.** The
inner minimax samples `[0,K]` at a fixed resolution and can step over a narrow
breach — the coarse-scan trap again, one level up. The first run of this campaign
duly reported `kd_app = 8.61` for a mesh whose gamma range is 6.11, which is
self-contradictory: if the tolerance held throughout `[0, 8.61]` the gamma range
could not be smaller than 8.61. The outer bisection now accepts a band only if THIS
function confirms it, and the reported number is THIS function's value.
"""
function applicable_range_multi(v; kd_max::Float64=400.0, n::Int=4000,
                                refine::Float64=1e-4)
    kd0 = 1e-3
    enorm(v, kd0) > 1.0 && return 0.0
    lo = kd0
    for kd in exp.(range(log(kd0), log(kd_max); length=n))
        if enorm(v, kd) > 1.0
            hi = kd
            while hi - lo > refine
                mid = 0.5*(lo+hi)
                enorm(v, mid) > 1.0 ? (hi = mid) : (lo = mid)
            end
            return lo
        end
        lo = kd
    end
    return kd_max
end

"""
    optimise_mesh(M, p; K) → NamedTuple

The σ-mesh design problem of `StokesWaveFourierAnalysis.tex` §sec: vertical grid
optimisation, in the form that section says is the correct one:

    (inner)  E(K) = min_{c₁<…<c_{M−1}}  max_{0≤kd≤K} |e(kd)|      Chebyshev, continuous
    (outer)  kd_app = max{K : E(K) ≤ tol}                          monotone bisection

**MAXIMISING `kd_app` DIRECTLY IS ILL-POSED, AND THIS SCRIPT ORIGINALLY DID IT.**
`kd_app` is a DISCONTINUOUS functional of the mesh: for `M=2`, `c₁=0.860` has an
interior dip that stays just inside 2 % and gives `kd_app=21.0`, while `c₁=0.8696`
has the same dip just breaching it and collapses to `7.7` — a factor of three from a
1 % change in the design variable. An optimiser pointed at that objective chases the
cliffs and returns FRAGILE designs, and this one duly returned `c₁=0.8696`, the
exact value the derivation names as the cliff. The minimax inner problem is
continuous in the design variables and equioscillates at the optimum; the outer
problem is a scalar bisection.

**AND THE OBJECTIVE MUST BE MULTI-PROPERTY.** Optimising phase celerity alone drives
the interface towards the free surface and DESTROYS the group velocity: for `M=2`,
`kd_app^(C)` rises monotonically to 15.5 at `c₁=0.81` while `kd_app^(Cg)` collapses
from 10.5 to 3.3 between `c₁=0.81` and `0.82`. A design excellent in celerity can be
useless for energy propagation. `enorm` therefore takes the worst of `C`, `C_g`, `γ`.

⚠ **`a₂` (the second-order bound harmonic) is NOT in the constraint set** — the
transfer function is not implemented in the solver. This matters for interpretation:
the published `c₁ = 0.728` is a band-`K≈6` design, chosen because `kd_app^(a₂)=6.0`
and extending the linear band past the nonlinear one buys accuracy the model cannot
use. A band-maximal design therefore does NOT reproduce 0.728, and should not be
expected to. Both are reported.
"""
function optimise_mesh(M::Int, p::Int; K_cap::Float64=400.0)
    if M == 1
        #  One element ⇒ no interfaces ⇒ NOTHING TO OPTIMISE. And by (P1) the
        #  interior nodes of an element do not affect any linear wave property at
        #  all — the design variables are the element interfaces and nothing else.
        c = [0.0, 1.0]
        v = assemble_dispersion_tensors(M, p, c)
        return (c_bdy=c, kd_app=applicable_range_multi(v),
                ranges=(applicable_range(v; prop=:C,     tol=TOLS.C),
                        applicable_range(v; prop=:Cg,    tol=TOLS.Cg),
                        applicable_range(v; prop=:gamma, tol=TOLS.gamma)),
                Emin=NaN, note="M=1: no interface, mesh fixed by construction")
    end
    z_from_c(c) = (w = diff(c); log.(w[2:end] ./ w[1]))
    #  inner problem at a fixed design band K
    function inner(K)
        obj(z) = (cc = cbdy_from_z(z, M);
                  vv = try assemble_dispersion_tensors(M, p, cc) catch; return Inf end;
                  minimax_error(vv, K))
        starts = Vector{Vector{Float64}}()
        haskey(YL_TABLE, M) && push!(starts, z_from_c(YL_TABLE[M].c))
        push!(starts, zeros(M-1))
        best_z, best_f = zeros(M-1), Inf
        for z0 in starts
            z, fv = nelder_mead(obj, collect(Float64, z0); step=0.30, maxiter=250)
            fv < best_f && (best_f = fv; best_z = z)
        end
        return (best_z, best_f)
    end
    #  Outer problem: bisect on the design band K. A band is ACCEPTED only if the
    #  mesh the inner problem returns for it actually achieves it under the FINE
    #  multi-property scan — not merely under the inner objective's sample grid,
    #  which can step over a narrow breach and certify a band the mesh lacks.
    lo, hi = 1.0, K_cap
    zbest  = inner(lo)[1]
    kdbest = applicable_range_multi(assemble_dispersion_tensors(M, p, cbdy_from_z(zbest, M)))
    for _ in 1:14                       # 14 geometric halvings of [1,400]
        mid = sqrt(lo*hi)
        zm  = inner(mid)[1]
        vm  = assemble_dispersion_tensors(M, p, cbdy_from_z(zm, M))
        kdm = applicable_range_multi(vm)
        kdm >= mid ? (lo = mid) : (hi = mid)
        #  Even a REJECTED band can hand back the best mesh seen: the inner optimum
        #  for an over-ambitious K is often still a fine design, and discarding it
        #  would make the result depend on the bisection path rather than the physics.
        kdm > kdbest && (zbest = zm; kdbest = kdm)
    end
    c = cbdy_from_z(zbest, M)
    v = assemble_dispersion_tensors(M, p, c)
    return (c_bdy=c, kd_app=kdbest,
            ranges=(applicable_range(v; prop=:C,     tol=TOLS.C),
                    applicable_range(v; prop=:Cg,    tol=TOLS.Cg),
                    applicable_range(v; prop=:gamma, tol=TOLS.gamma)),
            Emin=minimax_error(v, kdbest), note="")
end

const cfmt = cbstr        # one formatter, shared with the workers

meshes = Dict{Tuple{Int,Int},Any}()
vmesh_rows = String[]
if PHASE in ("1","both")
    println("\n" * "="^92)
    println("  PHASE 1 — optimised σ-meshes  (objective: maximise first-crossing applicable kd)")
    println("="^92); flush(stdout)

    # ---- THE CALIBRATION. Two standards, and BOTH must hold, because they fail
    #      independently and mean different things.
    #
    #      (a) THE PROPERTIES. Our C / C_g / γ applicable ranges, evaluated AT the
    #          published nodes, must reproduce Table 4.1. If this fails, our
    #          functionals are not theirs and nothing downstream has provenance.
    #      (b) THE OPTIMISER. The INNER minimax problem at a fixed design band must
    #          reproduce the band-dependent optima the derivation states:
    #          c₁ = 0.702 at K=5, 0.802 at K=10, 0.86 at K=20 for M=2.
    #
    #      ⚠ WHAT IS *NOT* A CALIBRATION: expecting a BAND-MAXIMAL design to return
    #      the published c₁=0.728. It will not, and should not. 0.728 is a design
    #      for band K≈6, chosen because the model's SECOND-ORDER range is
    #      kd_app^(a₂)=6.0 and extending the linear band beyond the nonlinear one
    #      buys accuracy a nonlinear computation cannot use. That is an engineering
    #      compromise, not a band-maximal optimum. An earlier version of this script
    #      treated the mismatch as a failure; it is a difference of objective.
    println("\n  --- calibration (a): our wave properties at the PUBLISHED nodes ---")
    calib_ok = true
    for M in (2,3,4)
        global calib_ok
        yl = YL_TABLE[M]
        vy = assemble_dispersion_tensors(M, 1, yl.c)
        got = (applicable_range(vy; prop=:C,     tol=TOLS.C),
               applicable_range(vy; prop=:Cg,    tol=TOLS.Cg),
               applicable_range(vy; prop=:gamma, tol=TOLS.gamma))
        okp = all(abs.(got .- yl.props) ./ yl.props .< 0.01)
        calib_ok &= okp
        @printf("    M=%d  C=%8.2f (%6.2f)  Cg=%7.2f (%6.2f)  gamma=%7.2f (%6.2f)  %s\n",
                M, got[1], yl.props[1], got[2], yl.props[2], got[3], yl.props[3],
                okp ? "OK" : "**MISMATCH**")
    end

    println("\n  --- calibration (b): inner minimax vs the published band-dependent optima ---")
    #  StokesWaveFourierAnalysis.tex §subsec: opt wellposedness: solving the INNER
    #  problem alone gives c₁ = 0.702 (K=5), 0.802 (K=10), 0.86 (K=20). Those are
    #  celerity-only figures, so this check uses the celerity-only objective — the
    #  point is that the OPTIMISER works, not that the constraint set is complete.
    for (K, cpub) in ((5.0, 0.702), (10.0, 0.802), (20.0, 0.86))
        global calib_ok
        objC(z) = (cc = cbdy_from_z(z, 2);
                   vv = try assemble_dispersion_tensors(2, 1, cc) catch; return Inf end;
                   maximum(property_errors(vv, kd)[1]
                           for kd in exp.(range(log(1e-3), log(K); length=120))))
        z, _ = nelder_mead(objC, [0.0]; step=0.30, maxiter=250)
        c1   = cbdy_from_z(z, 2)[2]
        okb  = abs(c1 - cpub) < 0.02
        calib_ok &= okb
        @printf("    K=%4.1f  ours c1=%.4f   published %.3f   |Δ|=%.4f  %s\n",
                K, c1, cpub, abs(c1-cpub), okb ? "OK" : "**MISMATCH**")
    end

    if !calib_ok
        println("\n  ⚠⚠ CALIBRATION FAILED — the optimiser or the wave properties do not")
        println("     reproduce the published standard, so the p ≥ 2 meshes below have no")
        println("     provenance. They are recorded and LABELLED UNCALIBRATED. Phase 2 is")
        println("     unaffected: node positions change the error CONSTANT, not the ORDER.")
    else
        println("\n  Calibration passed on both standards. The optimiser is a measuring")
        println("  instrument with a known reading against a known input.")
    end
    flush(stdout)

    println("\n  --- band-maximal multi-property optimum for every basis ---")
    println("     (kd_app is the MULTI-PROPERTY range: the worst of C, Cg, gamma)")
    @printf("    %-9s %3s  %-40s %9s | %8s %8s %8s\n",
            "basis","Nσ","c_bdy","kd_app","C","Cg","gamma")
    for b in BASES
        haskey(meshes, (b.M,b.p)) || (meshes[(b.M,b.p)] = optimise_mesh(b.M, b.p))
        r = meshes[(b.M,b.p)]
        @printf("    P%dLFE-%-3d %3d  %-40s %9.2f | %8.2f %8.2f %8.2f  %s\n",
                b.p, b.M, b.M*b.p+1, cfmt(r.c_bdy), r.kd_app,
                r.ranges[1], r.ranges[2], r.ranges[3], r.note)
        push!(vmesh_rows, csv_row(; phase="vmesh", basis="P$(b.p)LFE-$(b.M)",
                M=b.M, p_vert=b.p, Nsigma=b.M*b.p+1, c_bdy=cfmt(r.c_bdy),
                kd_app=@sprintf("%.4f", r.kd_app),
                kd_per_prop=@sprintf("C=%.2f Cg=%.2f g=%.2f", r.ranges...),
                verdict = calib_ok ? "OK" : "UNCALIBRATED",
                note = isempty(r.note) ? @sprintf("Emin=%.4f", r.Emin) : r.note))
        flush(stdout)
    end
    open(joinpath(OUT,"phase1_meshes.csv"),"w") do io
        println(io, "M,p_vert,Nsigma,c_bdy,kd_app_multi,kd_app_C,kd_app_Cg,kd_app_gamma,note")
        for b in BASES
            r = meshes[(b.M,b.p)]
            @printf(io, "%d,%d,%d,%s,%.4f,%.4f,%.4f,%.4f,%s\n", b.M, b.p, b.M*b.p+1,
                    cfmt(r.c_bdy), r.kd_app, r.ranges[1], r.ranges[2], r.ranges[3], r.note)
        end
    end
    println("\n  wrote $(joinpath(OUT,"phase1_meshes.csv"))")
end
PHASE == "1" && (println("\n  Phase 1 only — stopping here."); exit(0))

# ===========================================================================
#  PHASE 2 — the MMS convergence matrix
# ===========================================================================
println("\n" * "="^92)
println("  PHASE 2 — MMS convergence matrix  (spatial: fix dt, refine h | temporal: fix h, refine dt)")
println("="^92); flush(stdout)

#  The mesh each basis runs on. Phase 1's optimum when we have it; resolve_cbdy
#  otherwise — and it makes NO DIFFERENCE TO THE RATES either way, since node
#  positions change the error constant, not the order. Recorded per row regardless,
#  so every rate is traceable to the mesh it used without a join.
cbdy_for(M,p) = haskey(meshes,(M,p)) ? meshes[(M,p)].c_bdy : resolve_cbdy(M)
kdf_for(M,p)  = haskey(meshes,(M,p)) ? @sprintf("%.4f", meshes[(M,p)].kd_app) : ""
kdl_for(M,p)  = haskey(meshes,(M,p)) ?
                @sprintf("C=%.2f Cg=%.2f g=%.2f", meshes[(M,p)].ranges...) : ""

#  ---- the case list ------------------------------------------------------
#  SORTED BY (Nsigma, nl_pressure, model): Nσ lives in the FE value type, so each
#  distinct Nσ recompiles the whole residual/Jacobian stack, as does each pressure
#  tier. Compilation, not arithmetic, is this campaign's dominant cost — grouping
#  like with like is worth more than any scheduling cleverness.
cases = NamedTuple[]
for b in BASES, mno in MODELSEL
    m = MODELS[mno]
    push!(cases, (kind=:space, M=b.M, p=b.p, model=mno, integrator=:sdirk))
    push!(cases, (kind=:time,  M=b.M, p=b.p, model=mno, integrator=:theta))
    #  The production integrator is L-stable ⇒ DISSIPATIVE BY CONSTRUCTION, and the
    #  documented G7 measurements show it contaminating exactly this temporal
    #  window. :theta is therefore the primary temporal measurement — but leaving
    #  :sdirk entirely unmeasured would be a different blind spot, so it is run on
    #  the REFERENCE basis across all models, and reported separately, never merged.
    (b.M == 2 && b.p == 1) && push!(cases,
        (kind=:time, M=b.M, p=b.p, model=mno, integrator=:sdirk))
end
#  Sort by (Nσ, COST TIER, model, kind). ⚠ The tier key must be an explicit COST
#  RANK, not the symbol's name: sorting on `String(nlp)` orders them
#  "full" < "native" < "none" ALPHABETICALLY, i.e. it hands the workers the most
#  expensive tier FIRST and the cheapest last — the exact reverse of what you want
#  when the run may need to be cut short, and it is what the first attempt did.
costrank(m) = (m.regime === :linear ? 0 : 1) +
              (m.nlp === :none ? 0 : m.nlp === :native ? 2 : 4)
sort!(cases, by = c -> (c.M*c.p+1, costrank(MODELS[c.model]), c.model, String(c.kind)))
println("  $(length(cases)) studies queued  ($(count(c->c.kind===:space,cases)) spatial, " *
        "$(count(c->c.kind===:time,cases)) temporal)"); flush(stdout)

println("  $(nworkers()) worker processes ready"); flush(stdout)


#  Push the campaign's constants and the resolved meshes to every worker.
@everywhere const CFG = $(Dict(
    :levels_s=>LEVELS_S, :levels_t=>LEVELS_T, :nx0=>NX0, :nxt=>NXT,
    :nsteps_s=>NSTEPS_S, :ckpt=>CKPT))
@everywhere const CBDY = $(Dict((b.M,b.p)=>cbdy_for(b.M,b.p) for b in BASES))
@everywhere const MODELDEF = $MODELS
@everywhere const KDF = $(Dict((b.M,b.p)=>kdf_for(b.M,b.p) for b in BASES))
@everywhere const KDL = $(Dict((b.M,b.p)=>kdl_for(b.M,b.p) for b in BASES))

@everywhere function run_one(c)
    m   = MODELDEF[c.model]
    Nsg = c.M*c.p + 1
    cb  = CBDY[(c.M,c.p)]
    tag = "P$(c.p)LFE-$(c.M)"
    #  Nonlinear ⇒ quasi-Newton ⇒ LINEAR convergence ⇒ needs BUDGET, not a looser
    #  tolerance. Loosening nl_tol would bury the algebraic error inside the
    #  discretisation error the rate is measuring.
    nliter = m.regime === :linear ? 50 : 400
    hf = m.flat_bed ? nothing : bathymetry_field(; d0=_DEPTH, a_b=_AB, kbx=_KBX, kby=0.0)
    t0 = time()
    #  ⚠ LOG EVERY LEVEL, NOT JUST EVERY STUDY. A worker is otherwise a black box:
    #  the first two campaign attempts each ran >4 h on 8 workers and completed ONE
    #  study, and with `verbose=false` everywhere there was no way to tell a slow
    #  study from a stuck one without killing the run. Worker stdout is forwarded to
    #  the master, so these land in the campaign log with a "From worker N:" prefix.
    @printf("  [w%d] START %s M%d %s/%s\n", myid(), tag, c.model, c.kind, c.integrator)
    flush(stdout)
    local hs, ee, eu, nd, dts
    try
        hs = Float64[]; ee = Float64[]; eu = Float64[]; nd = Int[]; dts = Float64[]
        for l in 0:(c.kind === :space ? CFG[:levels_s] : CFG[:levels_t]) - 1
            #  SPATIAL: mode=:static ⇒ ω=0 ⇒ ∂ₜu* ≡ 0 and the temporal error is
            #  IDENTICALLY zero. No dt-independence guard needed — there is no dt
            #  error to guard against. dt/nsteps carried over unchanged from the
            #  verified nonlinear studies.
            #  TEMPORAL: ω=1.3 ⇒ period 4.83 s, and T_final=2.4 s is HALF A PERIOD.
            #  A short window is how G7 came to measure the spatial floor sitting
            #  still and call it a temporal rate.
            if c.kind === :space
                nx, ny = CFG[:nx0]*2^l, 3
                #  nsteps = 20, NOT the 100 inherited from `run_conv_study`'s default.
                #  MEASURED 2026-08-21 (M1, nx=16, ladder 100/20/10/5/2/1):
                #    e_eta  3.365584740738e-05 -> 3.365573857446e-05   (2.5e-6 rel)
                #    e_u    8.186772790758e-07 -> 7.592653028321e-07   (7.3e-2 rel)
                #    cost   302 s -> 1.2 s
                #  The problem is STATIC (ω=0), so this is a relaxation from the exact
                #  CONTINUOUS initial state to the DISCRETE steady state — and `e_eta`,
                #  the field with the tighter gate, is converged to SIX DIGITS after a
                #  single step. `e_u` is not: it drifts ~7 % and is still moving at
                #  nsteps=20, so the short run is slightly UNDER-RELAXED.
                #  ⚠ That 7 % is a real value change, accepted deliberately for a 15x
                #  saving, on the grounds that it is a smooth monotone offset applied
                #  IDENTICALLY at every refinement level and so cancels from the rate.
                #  It is NOT verified that it cancels — spot-check one basis at
                #  nsteps=100 before quoting any e_u CONSTANT from this campaign.
                dt, Tf, ω = 1e-5, 1e-5*CFG[:nsteps_s], 0.0
            else
                nx, ny = CFG[:nxt], 3
                dt = 0.15/2^l; Tf = 2.4; ω = 1.3
            end
            f = MMSField(Nsg; Lx=_LX, Ly=_LY, omega=ω, ky=0.0)   # ky=0 ⇒ quasi-1D
            r = run_mms_case(; nx=nx, ny=ny, dt=dt, T_final=Tf, Lx=_LX, Ly=_LY,
                               d=_DEPTH, M=c.M, p_vert=c.p, c_bdy=cb,
                               p_horizontal=_PU, p_eta=_PETA, field=f,
                               regime=m.regime, nl_pressure=m.nlp, flat_bed=m.flat_bed,
                               hfun=hf, solver_type=c.integrator, theta=0.5,
                               nl_tol=nltol_for(m.regime), nl_iter=nliter, verbose=false)
            #  ⚠ The spatial abscissa is Lx/nx, NOT `r.h`. `run_mms_case` returns
            #  `h = max(Lx/nx, Ly/ny)`, and this is a QUASI-1D study: ny is pinned at
            #  3 while nx refines, so Ly/ny = 0.367 dominates every level and `r.h`
            #  would be CONSTANT across the whole sequence — `convergence_rate` on a
            #  constant abscissa returns a meaningless slope from a division by zero
            #  in log-space. (`run_conv_study` pushes `Lx/nx` for exactly this reason.)
            push!(hs, c.kind === :space ? _LX/nx : dt); push!(dts, dt)
            push!(ee, r.e_eta); push!(eu, r.e_u); push!(nd, r.ndofs)
            @printf("  [w%d]   %s M%d %s L%d nx=%d dt=%.5g -> e_eta=%.3e e_u=%.3e (%.0f s)\n",
                    myid(), tag, c.model, c.kind, l, nx, dt, r.e_eta, r.e_u, time()-t0)
            flush(stdout)
        end
    catch e
        msg = first(split(sprint(showerror, e), '\n'))
        open(joinpath(CFG[:ckpt], "worker_$(myid()).csv"), "a") do io
            println(io, csv_row(; phase=String(c.kind), basis=tag, M=c.M, p_vert=c.p,
                Nsigma=Nsg, c_bdy=cbstr(cb), kd_app=KDF[(c.M,c.p)], kd_per_prop=KDL[(c.M,c.p)],
                model=c.model, regime=m.regime, flat_bed=m.flat_bed, nl_pressure=m.nlp,
                integrator=c.integrator, p_u=_PU, p_eta=_PETA,
                rate_gated_u=m.rate_u, verdict="ERROR", note=first(msg,160)))
        end
        @printf("  [w%d] ERROR %s M%d %s/%s after %.0f s: %s\n",
                myid(), tag, c.model, c.kind, c.integrator, time()-t0, first(msg,90))
        flush(stdout)
        return (ok=false, tag="$tag M$(c.model) $(c.kind)", msg=msg, secs=time()-t0)
    end
    fe, pwe = convergence_rate(hs, ee)
    fu, pwu = convergence_rate(hs, eu)
    #  Optimal rates follow the PAIRING, and the two fields differ: u ∈ Q_p ⇒ p+1,
    #  η ∈ Q_{p-1} ⇒ p. Temporal: 2 for both (SDIRK_2_2 and CN are both 2nd order).
    opte, optu = c.kind === :space ? (Float64(_PETA+1), Float64(_PU+1)) : (2.0, 2.0)
    #  READ THE LAST PAIRWISE RATE, not the fitted slope: the fit averages over
    #  pre-asymptotic levels, and a saturated study and a wrong coefficient give
    #  the same fitted number.
    okη = abs(pwe[end]-opte) < 0.3
    oku = m.rate_u ? abs(pwu[end]-optu) < 0.3 : true
    verdict = (okη && oku) ? "PASS" : "CHECK"
    open(joinpath(CFG[:ckpt], "worker_$(myid()).csv"), "a") do io
        for i in eachindex(hs)
            println(io, csv_row(; phase=String(c.kind), basis=tag, M=c.M, p_vert=c.p,
                Nsigma=Nsg, c_bdy=cbstr(cb), kd_app=KDF[(c.M,c.p)], kd_per_prop=KDL[(c.M,c.p)],
                model=c.model, regime=m.regime, flat_bed=m.flat_bed, nl_pressure=m.nlp,
                integrator=c.integrator, p_u=_PU, p_eta=_PETA, level=i-1,
                h = c.kind === :space ? hs[i] : _LX/CFG[:nxt], dt=dts[i], ndofs=nd[i],
                e_eta=@sprintf("%.10e",ee[i]), e_u=@sprintf("%.10e",eu[i]),
                #  pairwise rate i is the rate BETWEEN level i-1 and i, so level 0
                #  has none. Writing 0.0 there would be a number that means "no rate"
                #  and reads as "rate zero" — the one value this campaign must never
                #  fabricate, since slope 0 is the signature of a wrong coefficient.
                pw_eta = i == 1 ? "" : @sprintf("%.4f", pwe[i-1]),
                pw_u   = i == 1 ? "" : @sprintf("%.4f", pwu[i-1]),
                fit_eta=@sprintf("%.4f",fe), fit_u=@sprintf("%.4f",fu),
                opt_eta=opte, opt_u=optu, rate_gated_u=m.rate_u, verdict=verdict))
        end
    end
    @printf("  [w%d] DONE  %s M%d %s/%s  fit_eta=%.3f fit_u=%.3f  %s  (%.0f s)\n",
            myid(), tag, c.model, c.kind, c.integrator, fe, fu, verdict, time()-t0)
    flush(stdout)
    return (ok=true, tag="$tag M$(c.model) $(c.kind)/$(c.integrator)",
            fe=fe, fu=fu, pwe=pwe[end], pwu=pwu[end], opte=opte, optu=optu,
            e_u_last=eu[end], rate_u=m.rate_u, verdict=verdict, secs=time()-t0, Nsg=Nsg)
end

t_start = time()
#  ⚠ `batch_size` IS THE SINGLE BIGGEST COST LEVER HERE, AND THE DEFAULT IS WRONG
#  FOR THIS WORKLOAD. `Nσ` lives in the FE value type and each `nl_pressure` tier is
#  a different code path, so every distinct (Nσ, tier) combination triggers a FULL
#  JIT of the residual and both Jacobians — measured at ~10 min, against ~2 s/step
#  once warm. Compilation, not arithmetic, dominates.
#
#  With the default (dynamic, one case at a time) every worker walks the whole
#  sorted list and therefore compiles EVERY combination: ~12 combos × 8 workers ≈ 96
#  compilations. Handing each worker a CONTIGUOUS BLOCK of the (Nσ, tier)-sorted list
#  instead means a worker meets one or two combinations, not twelve — roughly 15
#  compilations total.
#
#  Two batches per worker rather than one: a single block per worker maximises reuse
#  but leaves no slack, so the worker holding the `:full` cases would run alone for
#  hours after the others finish. Halving the blocks buys back most of that.
batch = max(1, cld(length(cases), 2*nworkers()))
println("  batch_size=$batch (contiguous blocks of like (Nσ,tier) ⇒ ~1-2 JIT events/worker)")
flush(stdout)
results = pmap(cases; batch_size = batch,
               on_error = e -> (ok=false, tag="?", msg=sprint(showerror,e), secs=0.0)) do c
    r = run_one(c)
    return r
end

# ---- merge every checkpoint into THE single results file --------------------
merged = joinpath(OUT, "campaign_results.csv")
open(merged, "w") do io
    println(io, HDR)
    for r in vmesh_rows; println(io, r); end
    for f in sort(readdir(CKPT))
        endswith(f, ".csv") || continue
        for line in eachline(joinpath(CKPT, f)); println(io, line); end
    end
end

# ---- report -----------------------------------------------------------------
nok = count(r -> r.ok, results)
println("\n" * "="^116)
@printf("  CAMPAIGN COMPLETE — %d/%d studies returned, %.1f min wall\n",
        nok, length(results), (time()-t_start)/60)
println("="^116)
@printf("  %-28s %4s %-9s %-13s %-13s %-12s %s\n",
        "study","Nσ","verdict","eta fit/opt","u fit/opt","e_u(finest)","last pairwise η|u")
for r in results
    if !r.ok
        @printf("  %-28s %4s %-9s %s\n", r.tag, "-", "ERROR", get(r,:msg,""))
    else
        @printf("  %-28s %4d %-9s %5.2f/%-7.0f %5.2f/%-7s %.4e  %5.2f | %5.2f\n",
                r.tag, r.Nsg, r.verdict, r.fe, r.opte, r.fu,
                r.rate_u ? string(Int(r.optu)) : "floor", r.e_u_last, r.pwe, r.pwu)
    end
end
println("\n  single results file: $merged")
println("  checkpoints kept in : $CKPT")
println("\n  ⚠ Reading this table:")
println("    * the LAST PAIRWISE rate is the asymptotic one — the fit averages over")
println("      pre-asymptotic levels, and a saturated study and a wrong coefficient")
println("      produce the SAME fitted number;")
println("    * check the error MAGNITUDE before trusting a fine-level high-order rate;")
println("    * :full rows (models 7/8) are NOT rate-gated on u — compare e_u(finest)")
println("      ACROSS Nσ instead. That, not a rate, is what tier 3 asked;")
println("    * Phase 1 (dispersion) and Phase 2 (order) answer DIFFERENT questions.")
rmprocs(workers())
