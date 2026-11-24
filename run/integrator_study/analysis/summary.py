#!/usr/bin/env python3
"""summary.py — first-pass table of the Snellius integrator study (2026-11-23).

One row per case: outcome class, how far it got, energy, u_max, Newton count, cost.
Outcome classes:
  DONE       FINISHED marker, reached the target time, no divergence-guard warning
  GUARD      FINISHED but stopped by the divergence guard (|η| > 20 η_ref): a blow-up result
  NEWTON     FAILED with "Newton solver did not converge": a blow-up result (or stiff stall)
  QUOTA      killed by the Snellius home quota (LightXML write error / Disk quota / SIGTERM after it):
             NOT a result; the run was healthy up to t_end (read its diagnostics)
Usage: python3 -I summary.py <study_dir> run/integrator_study/cases.tsv > $IS_OUT/summary.tsv
(the 2026-11-23 analysis passed the then separate cases.tsv and cases_extra.tsv)
"""
import csv, os, re, sys, glob, math

study, *tables = sys.argv[1:]
meta = {}
for tab in tables:
    with open(tab) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["set"] = "extra" if "extra" in tab else "main"
            meta[r["name"]] = r

def read_diag(case):
    rows = {}
    for seg in sorted(glob.glob(os.path.join(study, case, "seg*"))):
        p = os.path.join(seg, "diagnostics.csv")
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for r in csv.DictReader(f):
                try:
                    rows[int(r["step"])] = r
                except (ValueError, KeyError):
                    pass
    return [rows[k] for k in sorted(rows)]

def fl(r, k):
    try:
        return float(r[k])
    except (ValueError, KeyError, TypeError):
        return float("nan")

def classify(case):
    d = os.path.join(study, case)
    logs = sorted(glob.glob(os.path.join(d, "seg*", "run.log")))
    txt = open(logs[-1], errors="replace").read() if logs else ""
    fail = open(os.path.join(d, "FAILED"), errors="replace").read() if os.path.exists(os.path.join(d, "FAILED")) else ""
    m = re.search(r"Diverged at t=([0-9.]+)", txt)
    if os.path.exists(os.path.join(d, "FINISHED")):
        return ("GUARD", float(m.group(1))) if m else ("DONE", None)
    both = fail + txt
    if "Newton solver did not converge" in both:
        mm = re.search(r"Final residual norm: ([0-9.eE+-]+|NaN)", both)
        return ("NEWTON", mm.group(1) if mm else "")
    if "XMLWriteError" in both or "Disk quota" in both or "signal (15)" in both:
        return ("QUOTA", None)
    return ("UNKNOWN", None)

cols = ["name", "set", "tier", "model", "T", "kA", "ncell", "dt", "integ", "gamma", "outcome", "detail",
        "t_end", "t_target", "frac", "E_end", "E_max", "t_Emax", "eta0", "eta_end", "eta_max", "u0", "u_end",
        "u_max", "t_umax", "nl_mean", "nl_last", "s_step", "cip_rate_end"]
print("\t".join(cols))
for case in sorted(meta):
    r = meta[case]
    integ = r["solver"] if r["solver"] != "sdirk" else r["tableau"]
    if r["solver"] == "gen_alpha":
        integ = f"galpha{r['rho']}"
    elif r["solver"] == "theta":
        integ = f"theta{r['theta']}"
    D = read_diag(case)
    oc, det = classify(case)
    T = float(r["T"]); tgt = float(r["periods"]) * T
    if not D:
        print("\t".join([case, r["set"], r["tier"], r["model"], r["T"], r["kA"], r["ncell"], r["dt"], integ,
                         r["gamma"], oc, str(det)] + ["nan"] * (len(cols) - 12)))
        continue
    t = [fl(x, "t") for x in D]
    E = [fl(x, "energy_ratio") for x in D]
    eta = [fl(x, "eta_max") for x in D]
    u = [fl(x, "u_max") for x in D]
    nl = [fl(x, "nl_iters") for x in D]
    ts = [fl(x, "t_solve") for x in D]
    finite = [i for i in range(len(E)) if math.isfinite(E[i])]
    iE = max(finite, key=lambda i: E[i]); iu = max(finite, key=lambda i: u[i])
    # s/step: t_solve is per diagnostics row (one row per save_every steps? check: steps column)
    steps = [int(x["step"]) for x in D]
    out = [case, r["set"], r["tier"], r["model"], r["T"], r["kA"], r["ncell"], r["dt"], integ, r["gamma"], oc,
           "" if det is None else str(det),
           f"{t[-1]:.2f}", f"{tgt:.1f}", f"{t[-1]/tgt:.3f}", f"{E[-1]:.4g}", f"{E[iE]:.4g}", f"{t[iE]:.1f}",
           f"{eta[0]:.4g}", f"{eta[-1]:.4g}", f"{max(eta[i] for i in finite):.4g}", f"{u[0]:.4g}",
           f"{u[-1]:.4g}", f"{u[iu]:.4g}", f"{t[iu]:.1f}",
           f"{sum(nl)/len(nl):.2f}", f"{nl[-1]:.0f}", f"{sorted(ts)[len(ts)//2]:.2f}",
           f"{fl(D[-1], 'cip_rate'):.3g}"]
    print("\t".join(out))
