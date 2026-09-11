# ==============================================================
#  run_vbasis_shard.jl — ONE SHARD of the Phase-2 convergence matrix
#
#  Plan: markdown_files/MMS_VBASIS_CAMPAIGN.md
#
#  WHY THIS EXISTS INSTEAD OF `pmap`. The `Distributed` version of this campaign
#  (`run_vbasis_campaign.jl`) was launched three times and completed 1 study in
#  ~4 h each time, while the IDENTICAL `run_mms_case` calls run at their expected
#  speed in a plain single process (measured: linear spatial level, nx=16,
#  20 steps → 19 s). The unit of work is not the problem; the harness around it is.
#  So: no `Distributed`, no worker stdout relay, no pmap batching. Each shard is an
#  ORDINARY Julia process that runs its slice of the case list and writes its own
#  CSV — the exact configuration whose throughput is known.
#
#  Two further things this buys, both of which bit the pmap version:
#    * stdout is THIS process's own, so `flush` actually reaches the log (under
#      `pmap` the workers' prints are relayed through the master, whose buffer is
#      never flushed — 35 min of a run with no way to tell progress from a hang);
#    * a shard that wedges can be killed and relaunched alone, without discarding
#      the other seven.
#
#  Shards take CONTIGUOUS blocks of the (Nσ, cost-tier)-sorted list, so each meets
#  one or two type combinations and pays JIT once or twice rather than twelve times.
#
#  ENV
#    VBC_SHARD    0-based shard index                     (required)
#    VBC_NSHARD   number of shards                        (required)
#    VBC_MESHES   phase1_meshes.csv from Phase 1          output/local/mms/vbasis_campaign_2026-08-30/phase1_meshes.csv
#    VBC_OUT      output directory                        output/local/mms/vbasis_campaign_2026-08-30
#    VBC_LEVELS_S / VBC_LEVELS_T / VBC_NX0 / VBC_NXT / VBC_NSTEPS_S / VBC_MODELS
#
#  RUN (all shards, one per core):
#    for i in $(seq 0 7); do
#      VBC_SHARD=$i VBC_NSHARD=8 julia --project=. \
#        examples/local_mms/run_vbasis_shard.jl > shard_$i.log 2>&1 &
#    done
# ==============================================================
using GridapBALFEM, Printf, Dates
const T_START = time()

gs(k,d)=get(ENV,k,d); gi(k,d)=parse(Int,get(ENV,k,string(d)))
const SHARD  = parse(Int, ENV["VBC_SHARD"])
const NSHARD = parse(Int, ENV["VBC_NSHARD"])
const OUT    = gs("VBC_OUT","output/local/mms/vbasis_campaign_2026-08-30")
const MESHCSV= gs("VBC_MESHES", joinpath(OUT,"phase1_meshes.csv"))
const LEVELS_S, LEVELS_T = gi("VBC_LEVELS_S",4), gi("VBC_LEVELS_T",4)
const NX0, NXT, NSTEPS_S = gi("VBC_NX0",8), gi("VBC_NXT",36), gi("VBC_NSTEPS_S",20)
const MODELSEL = parse.(Int, split(gs("VBC_MODELS","1,2,3,4,5,6"), ","))
#  ⚠ BOUND THE PROCESS LIFETIME. Measured 2026-08-29: shard processes left running
#  for days grew to 3.6-4.5 GB and fell to 7-14 % CPU — GC-bound, not compute-bound.
#  One spent 6.2 DAYS on an nx=8 level that takes seconds when fresh. A long-lived
#  Julia+Gridap process degrades catastrophically, so after each completed study a
#  worker that has exceeded this age exits cleanly and the supervisor restarts it.
#  Cost is one JIT (~3-5 min) per restart, against studies of 30 min to many hours.
const MAX_LIFE = parse(Float64, gs("VBC_MAX_LIFE_S","10800"))
#  ...and an RSS cap, which is the DIRECT trigger: the pathology is memory-driven,
#  and a time cap only bounds it indirectly (a 2.5 h cap was still letting processes
#  reach 3.4 GB). Whichever fires first recycles the worker.
const MAX_RSS_GB = parse(Float64, gs("VBC_MAX_RSS_GB","2.6"))
self_rss_gb() = try
    for ln in eachline("/proc/self/status")
        startswith(ln, "VmRSS:") && return parse(Float64, split(ln)[2]) / 1048576
    end
    0.0
catch; 0.0 end
#  :sjf maximises the NUMBER of studies finished per unit time — the right choice
#  when the deliverable is coverage. :lpt minimises makespan for a fixed set.
const ORDER   = gs("VBC_ORDER","sjf")
#  ERROR rows record a study that did not converge. They are "done" for bookkeeping
#  but NOT for science: with VBC_RETRY_ERRORS=1 they are re-queued (e.g. after a dt
#  ladder change) instead of being skipped.
const RETRY_E = gs("VBC_RETRY_ERRORS","0") != "0"
#  Per-basis dt0. A single ladder across bases was the wrong specification: at
#  dt=0.15 the Nσ=5 basis diverges outright (MMS_VBASIS_CAMPAIGN.md §2.3b), because
#  enriching the vertical basis widens the resolved band and stiffens the nonlinear
#  solve. Rates stay comparable across bases under different ladders — a rate is
#  scale-invariant — provided each basis is measured inside its own feasible regime.
#  VBC_DT0 overrides the rule outright — used to re-run a case that FAILED at its
#  ladder's coarsest step, so a non-convergence can be converted into a measurement
#  and attributed (too-large dt) rather than left as a bare error.
dt0_for(Nsg) = (v = gs("VBC_DT0","")) != "" ? parse(Float64, v) : (Nsg <= 3 ? 0.15 : 0.05)
mkpath(joinpath(OUT,"checkpoints"))

const LX, LY, DEPTH, A_BED, KBX = 1.7, 1.1, 2.5, 0.2, 1.3
const P_U, P_ETA = 3, 2
#  VBC_BASES="M:p,M:p,..." restricts the basis set — used to re-run a single basis
#  on a finer mesh ladder without disturbing the rest of the campaign.
const BASES = let spec = gs("VBC_BASES","2:1,1:2,3:1,4:1,2:2")
    [ (M=parse(Int,split(t,":")[1]), p=parse(Int,split(t,":")[2])) for t in split(spec,",") ]
end
const MODELS = Dict(
    1=>(regime=:linear,   flat_bed=true,  nlp=:none,   rate_u=true),
    2=>(regime=:linear,   flat_bed=false, nlp=:none,   rate_u=true),
    3=>(regime=:nonlinear,flat_bed=true,  nlp=:none,   rate_u=true),
    4=>(regime=:nonlinear,flat_bed=false, nlp=:none,   rate_u=true),
    5=>(regime=:nonlinear,flat_bed=true,  nlp=:native, rate_u=true),
    6=>(regime=:nonlinear,flat_bed=false, nlp=:native, rate_u=true),
    7=>(regime=:nonlinear,flat_bed=true,  nlp=:full,   rate_u=false),
    8=>(regime=:nonlinear,flat_bed=false, nlp=:full,   rate_u=false))
const COLS = [:phase,:basis,:M,:p_vert,:Nsigma,:c_bdy,:kd_app,:kd_per_prop,:model,:regime,
              :flat_bed,:nl_pressure,:integrator,:p_u,:p_eta,:level,:h,:dt,:ndofs,:e_eta,
              :e_u,:pw_eta,:pw_u,:fit_eta,:fit_u,:opt_eta,:opt_u,:rate_gated_u,:verdict,:note]
csv_row(; kwargs...) = begin
    d = Dict(kwargs)
    for k in keys(d); k in COLS || error("csv_row: unknown column :$k"); end
    join((replace(string(get(d, c, "")), "," => " ") for c in COLS), ",")
end
cbstr(c) = "[" * join(map(x -> @sprintf("%.4f", x), c), " ") * "]"

#  Phase-1 meshes, read from the file Phase 1 wrote. If it is missing the shard
#  FAILS rather than silently falling back to resolve_cbdy: a campaign whose rows
#  claim optimised meshes but used default ones is exactly the kind of mislabelled
#  result the A2 defect produced.
isfile(MESHCSV) || error("run_vbasis_shard: $MESHCSV not found — run Phase 1 first.")
const MESH = Dict{Tuple{Int,Int},Any}()
for (i,ln) in enumerate(eachline(MESHCSV))
    i == 1 && continue
    f = split(ln, ",")
    M, p = parse(Int,f[1]), parse(Int,f[2])
    cb = parse.(Float64, split(strip(f[4], ['[',']']), " "))
    MESH[(M,p)] = (c_bdy=cb, kd_app=parse(Float64,f[5]),
                   per=@sprintf("C=%.2f Cg=%.2f g=%.2f",
                        parse(Float64,f[6]), parse(Float64,f[7]), parse(Float64,f[8])))
end

#  VBC_KINDS restricts to `space` or `time`. Used to run the `:full` tier's SPATIAL
#  half on its own: that is the half tier 3 actually asks for (the e_u FLOOR as a
#  function of Nσ), and the temporal half of `:full` would cost ~4 h per study for
#  a number that is not a rate and that no gate reads.
const KINDS = gs("VBC_KINDS","both")
cases = NamedTuple[]
for b in BASES, mno in MODELSEL
    KINDS in ("both","space") &&
        push!(cases, (kind=:space, M=b.M, p=b.p, model=mno, integrator=:sdirk))
    if KINDS in ("both","time")
        push!(cases, (kind=:time,  M=b.M, p=b.p, model=mno, integrator=:theta))
        (b.M==2 && b.p==1) && push!(cases, (kind=:time, M=b.M, p=b.p, model=mno, integrator=:sdirk))
    end
end
costrank(m) = (m.regime===:linear ? 0 : 1) + (m.nlp===:none ? 0 : m.nlp===:native ? 2 : 4)
sort!(cases, by = c -> (c.M*c.p+1, costrank(MODELS[c.model]), c.model, String(c.kind)))
out = joinpath(OUT, "checkpoints", "shard_$(SHARD).csv")

#  ---- RESUME: skip anything already recorded by ANY shard --------------------
#  Makes a relaunch cheap, which is what allows the re-balance below to be worth
#  doing at all: the 10 studies already finished are kept, not recomputed.
#  VBC_DONEDIRS: extra output dirs to scan for already-completed studies. Needed
#  when a batch is re-tasked into a NEW directory (to avoid two processes appending
#  to the same checkpoint file) but must not redo what the original batch finished.
donekeys = Set{Tuple{String,String,Int,String}}()
for d in [OUT; split(gs("VBC_DONEDIRS",""), ",", keepempty=false)]
    ckd = joinpath(String(d), "checkpoints")
    isdir(ckd) || continue
    for f in readdir(ckd)
        endswith(f, ".csv") || continue
        for ln in eachline(joinpath(ckd,f))
            fs = split(ln, ",")
            length(fs) >= 29 || continue
            RETRY_E && fs[29] == "ERROR" && continue
            push!(donekeys, (fs[1], fs[2], parse(Int,fs[9]), fs[13]))
        end
    end
end
key(c) = (String(c.kind), "P$(c.p)LFE-$(c.M)", c.model, String(c.integrator))
todo = filter(c -> !(key(c) in donekeys), cases)

#  ---- CLAIMS: make IN-FLIGHT work visible to the other slots -------------------
#  A case only becomes "done" when it writes its rows, so a slot starting while
#  another is mid-study sees that study as available and runs it too. Harmless for
#  correctness (the report dedupes) but pure waste — and acutely so at the tail of a
#  campaign, where slots outnumber remaining studies and every one is expensive.
#  A claim is a file named for the case, holding the owner's PID; it is honoured
#  only while that PID is alive, so a killed or crashed worker releases its claim
#  automatically rather than blocking the case forever.
#  ⚠ CLAIMS MUST BE GLOBAL, not per-output-dir. Batches writing to different
#  VBC_OUT directories otherwise cannot see each other's claims, and the tail of the
#  campaign duly had one study running on three slots at once. VBC_CLAIMDIR defaults
#  to a shared path so every batch arbitrates against the same set.
const CLAIMDIR = gs("VBC_CLAIMDIR", "output/local/mms/_claims"); mkpath(CLAIMDIR)
claimfile(c) = joinpath(CLAIMDIR, replace("$(c.kind)_P$(c.p)LFE-$(c.M)_M$(c.model)_$(c.integrator)", "/"=>"_") * ".claim")
pid_alive(pid) = isdir("/proc/$pid")
function claimed_by_other(c)
    f = claimfile(c)
    isfile(f) || return false
    owner = tryparse(Int, strip(read(f, String)))
    owner === nothing && return false
    owner != getpid() && pid_alive(owner)
end
function take_claim!(c)
    #  re-check at the moment of taking, not only at assignment time
    claimed_by_other(c) && return false
    write(claimfile(c), string(getpid()))
    return true
end
nclaimed = count(claimed_by_other, todo)
todo = filter(c -> !claimed_by_other(c), todo)

#  ---- COST-BALANCED ASSIGNMENT (LPT), not a contiguous slice ------------------
#  ⚠ Contiguous slicing of a COST-SORTED list is the worst possible partition for
#  makespan: it concentrates every expensive case in the last shards. Measured on
#  the first attempt — shard 0 held ~1.7 h of work and shard 7 held ~21 h, a 12x
#  imbalance, so the machine would have sat mostly idle for the last ~15 h while
#  two shards ground on. Contiguity was chosen to amortise JIT, but JIT is ~3-5
#  min against studies of 30-200 min, so it is second-order and balance is not.
#
#  Longest-Processing-Time-first: sort by estimated cost descending, hand each
#  case to the least-loaded shard. Deterministic, so every shard computes the same
#  assignment independently and they cannot collide.
estmin(c) = begin
    m = MODELS[c.model]; lin = m.regime === :linear
    if c.kind === :space; lin ? 5 : (m.nlp === :none ? 30 : m.nlp === :native ? 90 : 120)
    else                ; lin ? 20 : (m.nlp === :none ? 90 : m.nlp === :native ? 200 : 260) end
end
order = sort(todo; by = estmin, rev = (ORDER == "lpt"))
load  = zeros(Float64, NSHARD); bins = [NamedTuple[] for _ in 1:NSHARD]
for c in order
    i = argmin(load); push!(bins[i], c); load[i] += estmin(c)
end
mine = bins[SHARD+1]

@printf("[shard %d/%d] %d cases assigned (%d done, %d in flight elsewhere, of %d), est %.1f h -> %s   %s\n",
        SHARD, NSHARD, length(mine), length(donekeys ∩ Set(key.(cases))), nclaimed, length(cases),
        load[SHARD+1]/60, out, Dates.now())
for c in mine
    @printf("[shard %d]   queued: P%dLFE-%d M%d %s/%s (~%d min)\n",
            SHARD, c.p, c.M, c.model, c.kind, c.integrator, estmin(c))
end
flush(stdout)

for (n, c) in enumerate(mine)
    m   = MODELS[c.model]
    Nsg = c.M*c.p + 1
    mesh = MESH[(c.M,c.p)]
    cb  = mesh.c_bdy
    tag = "P$(c.p)LFE-$(c.M)"
    nliter = m.regime === :linear ? 50 : 400
    #  1e-12 is UNREACHABLE for the nonlinear models: quasi-Newton by design ⇒
    #  linear convergence, stalls near 1e-10. See MMS_VBASIS_CAMPAIGN.md §2.4.
    nltol  = m.regime === :linear ? 1e-12 : 1e-9
    hf = m.flat_bed ? nothing : bathymetry_field(; d0=DEPTH, a_b=A_BED, kbx=KBX, kby=0.0)
    t0 = time()
    if !take_claim!(c)
        @printf("[shard %d] (%d/%d) SKIP %s M%d %s/%s — claimed by another slot\n",
                SHARD, n, length(mine), tag, c.model, c.kind, c.integrator)
        flush(stdout); continue
    end
    @printf("[shard %d] (%d/%d) START %s M%d %s/%s  %s\n",
            SHARD, n, length(mine), tag, c.model, c.kind, c.integrator, Dates.format(now(),"HH:MM:SS"))
    flush(stdout)
    hs = Float64[]; ee = Float64[]; eu = Float64[]; nd = Int[]; dts = Float64[]
    failed = ""
    try
        for l in 0:(c.kind === :space ? LEVELS_S : LEVELS_T) - 1
            if c.kind === :space
                nx, ny = NX0*2^l, 3
                dt, Tf, ω = 1e-5, 1e-5*NSTEPS_S, 0.0
            else
                nx, ny = NXT, 3
                #  T_final = 2.4 s is HALF the manufactured period (2π/1.3 ≈ 4.83 s)
                #  and must not shrink: a short window is how G7 came to measure the
                #  spatial floor sitting still and call it a temporal rate.
                dt = dt0_for(Nsg)/2^l; Tf = 2.4; ω = 1.3
            end
            f = MMSField(Nsg; Lx=LX, Ly=LY, omega=ω, ky=0.0)
            r = run_mms_case(; nx=nx, ny=ny, dt=dt, T_final=Tf, Lx=LX, Ly=LY, d=DEPTH,
                    M=c.M, p_vert=c.p, c_bdy=cb, p_u=P_U, p_eta=P_ETA, field=f,
                    regime=m.regime, nl_pressure=m.nlp, flat_bed=m.flat_bed, hfun=hf,
                    solver_type=c.integrator, theta=0.5,
                    nl_tol=nltol, nl_iter=nliter, verbose=false)
            #  Lx/nx, NOT r.h: r.h = max(Lx/nx, Ly/ny) and ny is PINNED at 3 here,
            #  so r.h would be constant across the whole sequence.
            push!(hs, c.kind === :space ? LX/nx : dt); push!(dts, dt)
            push!(ee, r.e_eta); push!(eu, r.e_u); push!(nd, r.ndofs)
            @printf("[shard %d]   L%d nx=%-3d dt=%-9.5g e_eta=%.4e e_u=%.4e  (%.0f s)\n",
                    SHARD, l, nx, dt, r.e_eta, r.e_u, time()-t0); flush(stdout)
        end
    catch e
        failed = first(split(sprint(showerror, e), '\n'))
    end
    base = (; phase=String(c.kind), basis=tag, M=c.M, p_vert=c.p, Nsigma=Nsg,
              c_bdy=cbstr(cb), kd_app=@sprintf("%.4f", mesh.kd_app), kd_per_prop=mesh.per,
              model=c.model, regime=m.regime, flat_bed=m.flat_bed, nl_pressure=m.nlp,
              integrator=c.integrator, p_u=P_U, p_eta=P_ETA, rate_gated_u=m.rate_u)
    open(out, "a") do io
        if !isempty(failed)
            println(io, csv_row(; base..., verdict="ERROR", note=first(failed,160)))
        else
            fe, pwe = convergence_rate(hs, ee); fu, pwu = convergence_rate(hs, eu)
            opte, optu = c.kind === :space ? (Float64(P_ETA+1), Float64(P_U+1)) : (2.0, 2.0)
            okη = abs(pwe[end]-opte) < 0.3
            oku = m.rate_u ? abs(pwu[end]-optu) < 0.3 : true
            verdict = (okη && oku) ? "PASS" : "CHECK"
            for i in eachindex(hs)
                println(io, csv_row(; base..., level=i-1,
                    h = c.kind === :space ? hs[i] : LX/NXT, dt=dts[i], ndofs=nd[i],
                    e_eta=@sprintf("%.10e",ee[i]), e_u=@sprintf("%.10e",eu[i]),
                    #  level 0 has no pairwise rate; "" not 0.0 — slope 0 is the
                    #  signature of a wrong coefficient and must never be fabricated.
                    pw_eta = i==1 ? "" : @sprintf("%.4f", pwe[i-1]),
                    pw_u   = i==1 ? "" : @sprintf("%.4f", pwu[i-1]),
                    fit_eta=@sprintf("%.4f",fe), fit_u=@sprintf("%.4f",fu),
                    opt_eta=opte, opt_u=optu, verdict=verdict))
            end
            @printf("[shard %d] (%d/%d) DONE  %s M%d %s/%s fit_eta=%.3f fit_u=%.3f %s (%.0f s)\n",
                    SHARD, n, length(mine), tag, c.model, c.kind, c.integrator, fe, fu, verdict, time()-t0)
        end
    end
    isempty(failed) || @printf("[shard %d] (%d/%d) ERROR %s M%d %s/%s after %.0f s: %s\n",
            SHARD, n, length(mine), tag, c.model, c.kind, c.integrator, time()-t0, first(failed,90))
    flush(stdout)
    rssgb = self_rss_gb()
    if time() - T_START > MAX_LIFE || rssgb > MAX_RSS_GB
        @printf("[shard %d] LIFETIME %.1f h / RSS %.1f GB EXCEEDED after %d/%d cases — exiting for a fresh process\n",
                SHARD, (time()-T_START)/3600, rssgb, n, length(mine))
        flush(stdout); exit(0)
    end
end
@printf("[shard %d] ALL %d CASES FINISHED  %s\n", SHARD, length(mine), Dates.now())
