# ==============================================================================
#  run_phaseB_shard.jl — the post-campaign MMS batch, sharded over independent
#  OS processes.
#
#  WHY INDEPENDENT PROCESSES, not Distributed/pmap: measured in the vertical-basis
#  campaign (MMS_VBASIS_CAMPAIGN.md §8.1) — three pmap launches each completed ONE
#  study in >4 h while the identical calls ran at full speed in a plain process,
#  and worker prints are relayed through a master whose stdout is never flushed,
#  so slow is indistinguishable from hung. Each shard here writes its OWN csv and
#  flushes every row.
#
#  SJF ORDERING: the deliverable is coverage, so cheap studies run first
#  (rule 43). Contiguous slicing of a cost-sorted queue is the WORST partition;
#  studies are dealt round-robin instead.
#
#  USAGE:  julia --project=. examples/local_mms/run_phaseB_shard.jl <shard> <nshard>
# ==============================================================================
using GridapBALFEM, Printf, Dates

const SHARD  = parse(Int, ARGS[1])
const NSHARD = parse(Int, ARGS[2])
const OUT    = joinpath(@__DIR__, "..", "..", "output", "local", "mms_phaseB")
mkpath(OUT)
const CSV = joinpath(OUT, @sprintf("shard_%02d.csv", SHARD))

MODELS = Dict(
 1=>(regime=:linear,    flat_bed=true,  nlp=:none),
 2=>(regime=:linear,    flat_bed=false, nlp=:none),
 3=>(regime=:nonlinear, flat_bed=true,  nlp=:none),
 4=>(regime=:nonlinear, flat_bed=false, nlp=:none),
 5=>(regime=:nonlinear, flat_bed=true,  nlp=:native),
 6=>(regime=:nonlinear, flat_bed=false, nlp=:native),
 7=>(regime=:nonlinear, flat_bed=true,  nlp=:full),
 8=>(regime=:nonlinear, flat_bed=false, nlp=:full))

#  ---- the queue. `cost` is a rough ordering key only (SJF). ------------------
struct Job; task::String; basis::String; M::Int; p::Int; model::Int
       p_u::Int; domain::Symbol; levels::Int; nx0::Int; a_eta::Float64; cost::Int; end

jobs = Job[]
#  TASK 7 — P1LFE-4 extended spatial ladder to nx=128 (1-D, Q3/Q2)
for m in 1:6
    push!(jobs, Job("T7_extladder","P1LFE-4",4,1,m,3,:d1,5,8,0.8, 40))
end
#  TASK 8 — the :full pair with the FIXED driver, at a solvable amplitude
for (nm,M,p) in (("P1LFE-2",2,1),("P2LFE-1",1,2),("P1LFE-3",3,1),("P1LFE-4",4,1),("P2LFE-2",2,2))
    for m in (7,8)
        push!(jobs, Job("T8_full_projection",nm,M,p,m,3,:d1,4,8,0.4, 60))
    end
end
#  TASK 9 — tier 2: the horizontal pairings. 1-D on a LONG ladder (this is what
#  actually answers the pre-asymptotic question, and it is affordable); 2-D at the
#  scheduled 4 levels.
for pu in (2,3,4), m in (1,3)
    push!(jobs, Job("T9_tier2_1d","P1LFE-2",2,1,m,pu,:d1,5,8,0.8, 30))
    push!(jobs, Job("T9_tier2_2d","P1LFE-2",2,1,m,pu,:d2,4,8,0.8, 90))
end
#  TASK 10 — high-order vertical bases on the NEW optimised nodes
for (nm,M,p) in (("P2LFE-1",1,2),("P2LFE-2",2,2)), m in 1:6
    push!(jobs, Job("T10_highorder_newnodes",nm,M,p,m,3,:d1,4,8,0.8, 35))
end

sort!(jobs; by = j -> j.cost)                    # SJF
mine = [j for (i,j) in enumerate(jobs) if (i-1) % NSHARD == SHARD]   # round-robin

#  ---- RESUME (rule 41: bound worker lifetime; a restart must cost one JIT, not
#  the whole shard). A study is DONE if any shard csv already carries its key.
#  Without this a recycled worker re-runs finished studies, which is exactly what
#  made killing a degraded process expensive the first time round.
jobkey(j) = (j.task, j.basis, j.model, j.p_u, string(j.domain))
done = Set{Tuple{String,String,Int,Int,String}}()
for f in readdir(OUT; join=true)
    endswith(f, ".csv") || continue
    for (i,ln) in enumerate(eachline(f))
        i == 1 && continue
        c = split(ln, ","); length(c) < 20 && continue
        push!(done, (c[1], c[2], parse(Int,c[5]), parse(Int,c[9]), c[10]))
    end
end
todo = [j for j in mine if !(jobkey(j) in done)]
@printf("[shard %02d] %d assigned, %d already done, %d to run\n",
        SHARD, length(mine), length(mine)-length(todo), length(todo)); flush(stdout)

#  Exit after MAX_STUDIES so the supervisor gets a FRESH process. Measured in the
#  vertical-basis campaign: RSS grew 1.5 -> 4.6 GB and CPU fell to 7-14 % on
#  long-lived Julia+Gridap workers. Bounding lifetime is the fix; resume makes it cheap.
const MAX_STUDIES = parse(Int, get(ENV, "PHASEB_MAX_STUDIES", "2"))
mine = todo[1:min(end, MAX_STUDIES)]
isempty(mine) && (println("[shard $SHARD] nothing to do"); exit(0))

newfile = !isfile(CSV)
open(CSV, "a") do io
    newfile && println(io,"task,basis,M,p,model,regime,flat_bed,nl_pressure,p_u,domain,levels,a_eta,",
               "pw_eta,pw_u,fit_eta,fit_u,e_eta_fine,e_u_fine,nx_fine,status,seconds,note")
    flush(io)
    for (n,j) in enumerate(mine)
        mo = MODELS[j.model]
        @printf("[shard %02d] %d/%d  %s %s M%d %s Q%d/Q%d %s\n",
                SHARD,n,length(mine),j.task,j.basis,j.model,string(mo.nlp),j.p_u,j.p_u-1,
                string(j.domain)); flush(stdout)
        t0 = time(); status="OK"; note=""
        local r
        try
            r = run_conv_study(p_u=j.p_u, domain=j.domain, mode=:static,
                               levels=j.levels, nx0=j.nx0,
                               M=j.M, p_vert=j.p, a_eta=j.a_eta,
                               regime=mo.regime, flat_bed=mo.flat_bed,
                               nl_pressure=mo.nlp,
                               a_b = mo.flat_bed ? 0.0 : 0.2,
                               verbose=false)
        catch e
            status="ERROR"; note=first(split(sprint(showerror,e),"\n"))
            note=replace(note, ","=>";")[1:min(end,120)]
            r=nothing
        end
        dt = round(time()-t0; digits=1)
        if r === nothing
            @printf(io,"%s,%s,%d,%d,%d,%s,%s,%s,%d,%s,%d,%.3f,,,,,,,%s,%.1f,%s\n",
                    j.task,j.basis,j.M,j.p,j.model,mo.regime,mo.flat_bed,mo.nlp,
                    j.p_u,j.domain,j.levels,j.a_eta,status,dt,note)
        else
            pwe = length(r.pw_eta)>0 ? r.pw_eta[end] : NaN
            pwu = length(r.pw_u)>0   ? r.pw_u[end]   : NaN
            @printf(io,"%s,%s,%d,%d,%d,%s,%s,%s,%d,%s,%d,%.3f,%.4f,%.4f,%.4f,%.4f,%.6e,%.6e,%d,%s,%.1f,%s\n",
                    j.task,j.basis,j.M,j.p,j.model,mo.regime,mo.flat_bed,mo.nlp,
                    j.p_u,j.domain,j.levels,j.a_eta,pwe,pwu,r.fit_eta,r.fit_u,
                    r.e_eta[end],r.e_u[end],j.nx0*2^(j.levels-1),status,dt,note)
        end
        flush(io)
    end
end
@printf("[shard %02d] DONE %d studies\n", SHARD, length(mine)); flush(stdout)
