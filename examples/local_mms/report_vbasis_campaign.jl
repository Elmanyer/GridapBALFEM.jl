# ==============================================================
#  report_vbasis_campaign.jl — merge the campaign into ONE results file
#  and render the summary tables.
#
#  Plan: building_files/MMS_VBASIS_CAMPAIGN.md
#
#  Merges: the Phase-1 mesh rows, the models 1-6 shard checkpoints, and the
#  models 7-8 `:full` batch (a separate output dir, run concurrently) into a single
#  `campaign_results.csv` — one row per refinement level, keyed by `phase`.
#
#  ⚠ REPORTING RULES BUILT IN HERE, not left to the reader:
#    * the LAST PAIRWISE rate is the asymptotic one; the FITTED slope averages over
#      pre-asymptotic levels. Both are printed, pairwise first. A saturated study and
#      a wrong coefficient give the SAME fitted number.
#    * models 7/8 (`:full`) are NOT rate-gated on u — their `e_u` FLOOR is reported
#      instead, and compared across Nσ, which is what tier 3 actually asks.
#    * `:sdirk` temporal rows are reported SEPARATELY from `:theta` and never merged:
#      the production integrator is L-stable, i.e. dissipative by construction.
#
#  RUN:  julia --project=. examples/local_mms/report_vbasis_campaign.jl
# ==============================================================
using Printf

const MAIN = get(ENV,"VBC_OUT","output/local/mms_campaign")
const FULL = get(ENV,"VBC_OUT_FULL","output/local/mms_campaign_full")
const SPAT = get(ENV,"VBC_OUT_SPATIAL","output/local/mms_campaign_spatial")
const HDR  = "phase,basis,M,p_vert,Nsigma,c_bdy,kd_app,kd_per_prop,model,regime,flat_bed," *
             "nl_pressure,integrator,p_u,p_eta,level,h,dt,ndofs,e_eta,e_u,pw_eta,pw_u," *
             "fit_eta,fit_u,opt_eta,opt_u,rate_gated_u,verdict,note"

rows = String[]
#  Phase 1 mesh rows, rebuilt from phase1_meshes.csv so the single results file is
#  self-contained: every Phase-2 rate is traceable to the mesh it ran on.
mfile = joinpath(MAIN,"phase1_meshes.csv")
meshinfo = Tuple[]
if isfile(mfile)
    for (i,ln) in enumerate(eachline(mfile))
        i == 1 && continue
        f = split(ln,",")
        push!(meshinfo, (parse(Int,f[1]), parse(Int,f[2]), parse(Int,f[3]), f[4],
                         parse(Float64,f[5]), parse(Float64,f[6]),
                         parse(Float64,f[7]), parse(Float64,f[8])))
        push!(rows, join(["vmesh","P$(f[2])LFE-$(f[1])",f[1],f[2],f[3],f[4],
              @sprintf("%.4f",parse(Float64,f[5])),
              @sprintf("C=%.2f Cg=%.2f g=%.2f",parse(Float64,f[6]),parse(Float64,f[7]),parse(Float64,f[8])),
              #  fields 9..28 empty, then verdict(29), note(30). `fill("",21)` here
              #  produced 31-field rows — the exact positional-miscount hazard that
              #  `csv_row` exists to prevent; the assertion below now catches it.
              fill("",20)..., "OK",""], ","))
    end
end
for d in (MAIN, FULL, SPAT)
    ck = joinpath(d,"checkpoints")
    isdir(ck) || continue
    for f in sort(readdir(ck))
        endswith(f,".csv") || continue
        append!(rows, [ln for ln in eachline(joinpath(ck,f)) if !isempty(strip(ln))])
    end
end
merged = joinpath(MAIN,"campaign_results.csv")
#  ⚠ EVERY row must have exactly length(COLS) fields. A miscount produces a file
#  that parses cleanly and means something else from that column on — the failure
#  mode this whole campaign has been bitten by twice. Fail loudly instead.
const NCOL = 30
for (i,r) in enumerate(rows)
    n = count(==(','), r) + 1
    n == NCOL || error("campaign_results row $i has $n fields, expected $NCOL:\n  $r")
end

# ---- parse back for the summary ------------------------------------------------
#  DEDUPE by (phase, basis, model, integrator, level), keeping the FIRST occurrence.
#  Shards are re-partitioned during a long run to fix priority inversions, and a
#  re-partition can hand the same case to two shards. A duplicate is wasted compute,
#  never a wrong answer — but it must not appear twice in the results file, and a
#  study must not be counted twice in the summary.
recs0 = [split(r,",") for r in rows if !startswith(r,"vmesh")]
seenrow = Set(); recs = Vector{Vector{SubString{String}}}()
for r in recs0
    k = (r[1], r[2], r[9], r[13], r[16])
    k in seenrow && continue
    push!(seenrow, k); push!(recs, r)
end
ndup = length(recs0) - length(recs)
ndup > 0 && @printf("  (deduped %d duplicate row(s) from overlapping shard assignments)\n", ndup)
key(r) = (r[1], r[2], r[9], r[13])          # phase, basis, model, integrator
studies = Dict{Any,Vector{Vector{SubString{String}}}}()
for r in recs; push!(get!(studies, key(r), []), r); end

#  written AFTER dedup, so the single results file carries no duplicates
open(merged,"w") do io
    println(io,HDR)
    for r in rows; startswith(r,"vmesh") && println(io,r); end
    for r in recs; println(io, join(r, ",")); end
end

fmt(x) = isempty(strip(x)) ? "  -  " : @sprintf("%5.3f", parse(Float64,x))
println("\n" * "="^118)
println("  VERTICAL-BASIS CONVERGENCE CAMPAIGN — SUMMARY")
println("="^118)

if !isempty(meshinfo)
    println("\n### PHASE 1 — optimised σ-meshes (multi-property: C, C_g, γ; γ binds throughout)\n")
    @printf("  %-10s %3s  %-38s %9s | %8s %8s %8s\n","basis","Nσ","c_bdy","kd_app","C","Cg","gamma")
    for (M,p,Ns,cb,ka,kc,kg,kga) in meshinfo
        @printf("  P%dLFE-%-4d %3d  %-38s %9.2f | %8.2f %8.2f %8.2f\n", p,M,Ns,cb,ka,kc,kg,kga)
    end
end

for (ph, ttl, oe, ou) in (("space","SPATIAL  (fix dt, refine h)  optimal: p_eta 3, p_u 4", 3, 4),
                          ("time", "TEMPORAL (fix h, refine dt)  optimal: 2 and 2",        2, 2))
    for integ in ("sdirk","theta")
        sel = sort([k for k in keys(studies) if k[1]==ph && k[4]==integ], by=k->(k[2],k[3]))
        isempty(sel) && continue
        println("\n### $ttl   —   integrator :$integ")
        integ == "sdirk" && ph == "time" && println(
            "    ⚠ L-stable ⇒ DISSIPATIVE BY CONSTRUCTION. Reported separately from :theta,\n" *
            "      never merged with it. A depressed rate here is the scheme, not the operator.")
        @printf("    %-10s %3s %5s  %-9s %-9s %-9s %-9s %-8s %s\n",
                "basis","Nσ","model","pw_eta","pw_u","fit_eta","fit_u","e_u(fin)","verdict")
        for k in sel
            rs = sort(studies[k], by=r->parse(Int,r[16]=="" ? "0" : r[16]))
            last = rs[end]
            if last[29] == "ERROR"
                @printf("    %-10s %3s %5s  %s\n", k[2], last[5], k[3],
                        "ERROR: " * String(last[30]))
            else
                gated = last[28] == "true"
                @printf("    %-10s %3s %5s  %-9s %-9s %-9s %-9s %-8.2e %s%s\n",
                        k[2], last[5], k[3], fmt(last[22]), gated ? fmt(last[23]) : " floor",
                        fmt(last[24]), fmt(last[25]), parse(Float64,last[21]), last[29],
                        gated ? "" : "   (:full — u not rate-gated)")
            end
        end
    end
end

#  Tier 3: the :full floor as a function of Nσ. THE question, and it is a value
#  comparison across bases, never a rate.
fl = [r for r in recs if r[12]=="full" && r[29]!="ERROR" && r[1]=="space"]
if !isempty(fl)
    println("\n### TIER 3 — the `:full` e_u FLOOR vs Nσ  (a VALUE comparison, NOT a rate)")
    println("    Does the frozen-projection/omission floor depend on the vertical resolution?")
    seen = Set()
    for r in sort(fl, by=r->(parse(Int,r[5]), r[9], parse(Int,r[16])))
        k = (r[2], r[9])
        if !(k in seen) || parse(Int,r[16]) == 3
            parse(Int,r[16]) == 3 && @printf("    %-10s Nσ=%s model %s  e_u(finest)=%.4e   p_u=%s\n",
                    r[2], r[5], r[9], parse(Float64,r[21]), r[25])
            push!(seen,k)
        end
    end
end

nerr = count(r -> r[29]=="ERROR", recs)
println("\n" * "="^118)
@printf("  %d studies, %d error(s).  Single results file: %s\n",
        length(studies), nerr, merged)
println("="^118)
println("""
  Reading this table:
    * pw_* is the LAST PAIRWISE rate — the asymptotic one. fit_* averages over
      pre-asymptotic levels; where they disagree, believe the pairwise sequence.
    * check the ERROR MAGNITUDE before trusting any fine-level high-order rate.
    * :full rows carry a FLOOR, not a rate, and the floor is compared across Nσ.
    * Phase 1 (dispersion) and Phase 2 (order) answer DIFFERENT questions and
      neither is evidence for the other.""")
