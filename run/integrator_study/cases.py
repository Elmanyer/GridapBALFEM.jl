#!/usr/bin/env python3
"""cases.py — the case matrix of the time-integrator stability study (STABILITY.md §6b).

    python3 run/integrator_study/cases.py            writes cases.tsv and packs.tsv next to this file

The closed x-periodic box (examples/local_1d/run_periodic_1d.jl), P1LFE-2, Q3/Q2, one model
wavelength, seed 1e-8, 200 periods, 40 steps per period (dt = T/40, so ωΔt = 0.157 and every
integrator damps the carrier by the same amount per period at every T). The study asks how
much the full nonlinear solver's growth is held by integrator dissipation, so it leans on high
steepness: kA ∈ {0.158 (the campaign), 0.25, 0.32 (72 % of the Miche limit in deep water)}.

TIERS (submit separately; the pilot first):
  P  pilot: m5 n32 T=1.6 kA=0.25, every integrator, 100 steps — s/step on the node
  A  core:  m5 (𝓝, flat) n32, T ∈ {1.6, 2.2, 3.0} (kd ≈ 5.5, 3.0, 1.7), kA ∈ {0.158, 0.25, 0.32}
  B  resolution: m5 T=1.6 kA ∈ {0.25, 0.32} at n16 and n64 (n32 is in A)
  C  variable bed: m6 (cosine bed, max slope 0.1) n32 T=1.6 kA ∈ {0.25, 0.32}
  D  time step: m5 n32 T=1.6 kA=0.32, every integrator at dt = T/20 and T/80
  E  null control: m3 (𝓝 off, no CIP) n32, T ∈ {1.6, 3.0}, every kA — cheap
Every configuration runs under every integrator, Crank–Nicolson included: CN is the diagnosis
pair of rule 15 (what grows without dissipation), the others are the remedies.

PACKS: a rome 1/8 node (16 cores, 28 GiB) runs 8 cases at once (≈ 3 GiB each for 𝓝). Cases of
one tier and model class are sorted by expected cost and grouped 8 at a time, so the cases of a
pack finish at similar times (a pack holds its allocation until its slowest case ends). One SLURM
array task = one pack. The cost weights are relative guesses for the grouping only; the pilot
measures the real s/step.

Columns of cases.tsv:
  name tier model nlp flat slope T kA ncell dt periods solver tableau theta rho cost
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PACK_SIZE = 8

#  tag: (BALFEM_SOLVER, BALFEM_TABLEAU, θ, ρ∞, relative cost per step vs CN)
#  The matrix runs five (author, 2026-11-11): TR-BDF2, Crouzeix 3/4 and generalised-α were dropped,
#  explicit RK3 added next to RK4.
INTEGRATORS = {
    "cn":      ("theta",     "SDIRK_2_2",           0.5, 0.5, 1.0),   # diagnosis: no damping
    "sdirk22": ("sdirk",     "SDIRK_2_2",           0.5, 0.5, 2.2),   # 2 stages, 0.36 y⁴, R(∞)=−½
    "sdirk43": ("sdirk",     "SDIRK_4_3",           0.5, 0.5, 4.4),   # L-stable, 3rd order, 4 stages
    "rk4":     ("sdirk",     "EXRK_RungeKutta_4_4", 0.5, 0.5, 4.0),   # Yang & Liu's scheme; explicit
    "rk3":     ("sdirk",     "EXRK_SSP_3_3",        0.5, 0.5, 3.0),   # explicit, y⁴/24, limit ωΔt < √3
}
PILOT_EXTRA = {}                 # integrators run in the pilot only (none at present)
MODELS = {   # BALFEM_NL_PRESSURE, BALFEM_FLAT_BED, BALFEM_BED_SLOPE, relative cost per step
    "m3": (0, 1, 0.0, 1.0),
    "m5": (1, 1, 0.0, 10.0),
    "m6": (1, 0, 0.1, 11.0),
}
PERIODS, SPP = 200, 40          # duration in periods, steps per period


def case(tier, model, T, kA, ncell, itag, spp=SPP, periods=PERIODS, suffix=""):
    nlp, flat, slope, cm = MODELS[model]
    solver, tableau, theta, rho, ci = {**INTEGRATORS, **PILOT_EXTRA}[itag]
    name = f"{model}_T{round(T * 100):03d}_kA{round(kA * 1000):03d}_n{ncell}_{itag}{suffix}"
    cost = cm * ci * periods * spp / 8000.0
    return dict(name=name, tier=tier, model=model, nlp=nlp, flat=flat, slope=slope, T=T, kA=kA,
                ncell=ncell, dt=round(T / spp, 6), periods=periods, solver=solver, tableau=tableau,
                theta=theta, rho=rho, cost=round(cost, 3))


def matrix():
    C = []
    for itag in {**INTEGRATORS, **PILOT_EXTRA}:                         # P: pilot
        C.append(case("P", "m5", 1.6, 0.25, 32, itag, periods=2.5, suffix="_pilot"))
    for T in (1.6, 2.2, 3.0):                                           # A: core
        for kA in (0.158, 0.25, 0.32):
            for itag in INTEGRATORS:
                C.append(case("A", "m5", T, kA, 32, itag))
    for n in (16, 64):                                                  # B: resolution
        for kA in (0.25, 0.32):
            for itag in INTEGRATORS:
                C.append(case("B", "m5", 1.6, kA, n, itag))
    for kA in (0.25, 0.32):                                             # C: variable bed
        for itag in INTEGRATORS:
            C.append(case("C", "m6", 1.6, kA, 32, itag))
    for spp, sfx in ((20, "_dt20"), (80, "_dt80")):                     # D: time step
        for itag in INTEGRATORS:
            C.append(case("D", "m5", 1.6, 0.32, 32, itag, spp=spp, suffix=sfx))
    for T in (1.6, 3.0):                                                # E: null control
        for kA in (0.158, 0.25, 0.32):
            for itag in INTEGRATORS:
                C.append(case("E", "m3", T, kA, 32, itag))
    names = [c["name"] for c in C]
    assert len(names) == len(set(names)), "duplicate case names"
    return C


def packs(C):
    """Within each (tier, model): the cost-sorted cases cut into the fewest packs of at most
    PACK_SIZE, as EQUAL in size as possible (10 cases → 5 + 5, not 8 + 2: same bill, since each
    pack rents a 1/8 node, but the packs end together). (pack id, tier, case names)."""
    out, pid = [], 0
    for tier in "PABCDE":
        for model in MODELS:
            group = sorted((c for c in C if c["tier"] == tier and c["model"] == model),
                           key=lambda c: -c["cost"])
            if not group:
                continue
            npk = -(-len(group) // PACK_SIZE)
            base, extra = divmod(len(group), npk)
            i = 0
            for k in range(npk):
                size = base + (1 if k < extra else 0)
                pid += 1
                out.append((pid, tier, [c["name"] for c in group[i:i + size]]))
                i += size
    return out


if __name__ == "__main__":
    C = matrix()
    cols = ["name", "tier", "model", "nlp", "flat", "slope", "T", "kA", "ncell", "dt", "periods",
            "solver", "tableau", "theta", "rho", "cost"]
    with open(os.path.join(HERE, "cases.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for c in C:
            f.write("\t".join(str(c[k]) for k in cols) + "\n")
    P = packs(C)
    with open(os.path.join(HERE, "packs.tsv"), "w") as f:
        f.write("pack\ttier\tcases\n")
        for pid, tier, names in P:
            f.write(f"{pid}\t{tier}\t{','.join(names)}\n")
    for tier in "PABCDE":
        nc = sum(c["tier"] == tier for c in C)
        npk = sum(p[1] == tier for p in P)
        cost = sum(c["cost"] for c in C if c["tier"] == tier)
        print(f"tier {tier}: {nc:3d} cases, {npk:2d} packs, relative cost {cost:7.1f}")
    print(f"total: {len(C)} cases, {len(P)} packs")
