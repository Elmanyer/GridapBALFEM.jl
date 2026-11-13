#!/usr/bin/env python3
"""cases.py — the case matrix of the time-integrator stability study (STABILITY.md §6b).

    python3 run/integrator_study/cases.py            writes cases.tsv and packs.tsv next to this file

The closed x-periodic box (examples/local_1d/run_periodic_1d.jl), P1LFE-2, Q3/Q2, one model
wavelength, seed 1e-8, 200 periods, 40 steps per period (dt = T/40, so ωΔt = 0.157 and every
integrator damps the carrier by the same amount per period at every T). The study asks how
much the full nonlinear solver's growth is held by integrator dissipation, so it leans on high
steepness: kA ∈ {0.158 (the campaign), 0.25, 0.32 (72 % of the Miche limit in deep water)}.

TIERS (submit separately; the pilot first). γ is the C⁰-CIP penalty (γ_u = γ_η); its default is
0.3 for m5/m6 (mandatory with 𝓝) and 0 for m3/m4 (no CIP: the weak form of models 3–4 has no
non-conforming term to stabilise; author, 2026-11-12).
  P  pilot: m5 n32 T=1.6 kA=0.25, the implicit integrators, 100 steps — s/step on the node
     (+ CN at γ = 1, the stiffest penalty in the matrix: does it slow Newton?)
  A  core:  m5 (𝓝, flat) n32, T ∈ {1.6, 2.2, 3.0} (kd ≈ 5.5, 3.0, 1.7), kA ∈ {0.158, 0.25, 0.32}
  B  resolution: m5 T=1.6 kA ∈ {0.25, 0.32} at n16 and n64 (n32 is in A)
  C  variable bed: m6 (cosine bed, max slope 0.1) n32 T=1.6 kA ∈ {0.25, 0.32}; kA = 0.25 also at
     γ ∈ {0.5, 1}
  D  time step: m5 n32 T=1.6 kA=0.32, every implicit integrator at dt = T/20 and T/80
  E  null control: m3 (𝓝 off, no CIP) n32, T ∈ {1.6, 3.0}, every kA — cheap; the only tier
     with the explicit pair rk3/rk4 at T/40
  F  m3 selective-damping ladder (hypothesis H1): m3 n32 T=1.6, every kA, the 4th-order schemes at
     a larger dt to aim their damping at the 2ω₀ member of the decay triad: Crouzeix 3/4 at T/40,
     T/20, T/10 and RK4 at T/20, T/10 (RK4 is affordable without CIP)
  G  γ map (hypothesis H2): m5 n32 T=1.6, every kA × γ ∈ GAMMAS_G × the four implicit integrators
     (γ = 0.3 is tier A)
  H  γ × dt (hypothesis H3): m5 n32 T=1.6 kA=0.32 × γ ∈ {0.5, 1} × {SDIRK_2_2, Crouzeix 3/4} at T/20
     and T/80 (γ = 0.3 is tier D)
The hypotheses, the reasons and the expected results: markdown_files/CH8_CAMPAIGN_TODO.md
("Integrator study: hypotheses"). Submission order: P → E F → A → G → C → B → D H.
Every configuration runs under every integrator it can afford, Crank–Nicolson included: CN is the
diagnosis pair of rule 15 (what grows without dissipation), the others are the remedies.
EXPLICIT SCHEMES ONLY WITHOUT CIP (author, 2026-11-12). In the smoke test every explicit tableau
failed on m5 at Δt = T/40 (energy ×10⁶ in one RK4 step), and RK4 failed on m6. (The Δt limits once
quoted here, from CIP damping modes down to −9.5e3 s⁻¹, are void: stability_eig.jl's FD Jacobian
missed 17 % of the CIP entries until 2026-11-13; to be re-measured.) So rk3/rk4 run in tier E
only (m3, no CIP, fastest discrete mode ≈ 9 s⁻¹, imaginary: T/40 is well inside both limits). In
their place m5/m6 run SDIRK_Crouzeix_3_4
(4th order like RK4, damping rising ≈ y⁶; author, 2026-11-12).

PACKS: a rome 1/8 node (16 cores, 28 GiB) runs 8 cases at once (≈ 3 GiB each for 𝓝). Cases of
one tier and model class are sorted by expected cost and grouped 8 at a time, so the cases of a
pack finish at similar times (a pack holds its allocation until its slowest case ends). One SLURM
array task = one pack. The cost weights are relative guesses for the grouping only; the pilot
measures the real s/step.

OUTPUT CADENCE follows the period, not the step: VTK and diagnostics every max(1, spp ÷ 8) steps
(8–10 samples per period, as the campaign's temporal-spectrum analysis needs: at 2 samples per
period, 2ω₀ would alias to 0) and a checkpoint every period (spp steps): a segment cut by the
120 h wall loses at most one period of work (≈ 1.5–3 h for the costliest cases).

Columns of cases.tsv:
  name tier model nlp flat slope T kA ncell dt periods solver tableau theta rho gamma save_every
  ckpt_every cost
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
    "crz34":   ("sdirk",     "SDIRK_Crouzeix_3_4",  0.5, 0.5, 3.3),   # 4th order, 3 stages, R(∞)=0.63:
}                                                                     #   the implicit stand-in for RK4
PILOT_EXTRA = {}                 # integrators run in the pilot only (none at present)
EXPLICIT = {"rk3", "rk4"}        # unaffordable with the C⁰-CIP (m5/m6): see the docstring
FULL_ONLY = {"crz34"}            # m5/m6 only (author, 2026-11-12): RK4's implicit stand-in where RK4
                                 #   cannot run; m3 has RK4 itself
MODELS = {   # BALFEM_NL_PRESSURE, BALFEM_FLAT_BED, BALFEM_BED_SLOPE, relative cost per step
    "m3": (0, 1, 0.0, 1.0),
    "m5": (1, 1, 0.0, 10.0),
    "m6": (1, 0, 0.1, 11.0),
}
TIERS = "PABCDEFGH"
PERIODS, SPP = 200, 40          # duration in periods, steps per period
GAMMA_DEFAULT = {0: 0.0, 1: 0.3}  # by nlp: no CIP without 𝓝 (m3/m4), DEFAULT_CIP_GAMMA with it
GAMMAS_G = (0.1, 0.5, 1.0)      # tier G (+ 0.3 in tier A). Set from pre-cluster check 2 (2026-11-13):
                                #   the steady grid-scale mode leads below γ ≈ 0.2 (kA 0.25) / 0.5 (kA 0.32)
                                #   and nothing changes above γ ≈ 1, so γ = 3 was dropped (author)


def integrators(model, extra=False):
    """The integrator tags a model runs: with 𝓝 (hence the CIP) the implicit ones, Crouzeix 3/4
    included; without it everything but the m5/m6-only Crouzeix 3/4."""
    tags = {**INTEGRATORS, **PILOT_EXTRA} if extra else INTEGRATORS
    full = MODELS[model][0] == 1
    return [t for t in tags if not ((full and t in EXPLICIT) or (not full and t in FULL_ONLY))]


def case(tier, model, T, kA, ncell, itag, spp=SPP, periods=PERIODS, suffix="", gamma=None):
    """One case. `gamma` None = the model's default (GAMMA_DEFAULT); a non-default γ adds `_g<100γ>`
    to the name, a non-default spp adds `_dt<spp>`."""
    nlp, flat, slope, cm = MODELS[model]
    solver, tableau, theta, rho, ci = {**INTEGRATORS, **PILOT_EXTRA}[itag]
    g0 = GAMMA_DEFAULT[nlp]
    gamma = g0 if gamma is None else gamma
    assert nlp == 0 or gamma > 0, "with 𝓝 the C⁰-CIP is mandatory (γ > 0)"
    gsfx = "" if gamma == g0 else f"_g{round(gamma * 100):03d}"
    dsfx = "" if spp == SPP else f"_dt{spp}"
    name = f"{model}_T{round(T * 100):03d}_kA{round(kA * 1000):03d}_n{ncell}_{itag}{dsfx}{gsfx}{suffix}"
    cost = cm * ci * periods * spp / 8000.0
    return dict(name=name, tier=tier, model=model, nlp=nlp, flat=flat, slope=slope, T=T, kA=kA,
                ncell=ncell, dt=round(T / spp, 6), periods=periods, solver=solver, tableau=tableau,
                theta=theta, rho=rho, gamma=gamma, save_every=max(1, spp // 8),
                ckpt_every=spp, cost=round(cost, 3))


def matrix():
    C = []
    for itag in integrators("m5", extra=True):                          # P: pilot
        C.append(case("P", "m5", 1.6, 0.25, 32, itag, periods=2.5, suffix="_pilot"))
    C.append(case("P", "m5", 1.6, 0.25, 32, "cn", periods=2.5, suffix="_pilot", gamma=1.0))
    for T in (1.6, 2.2, 3.0):                                           # A: core
        for kA in (0.158, 0.25, 0.32):
            for itag in integrators("m5"):
                C.append(case("A", "m5", T, kA, 32, itag))
    for n in (16, 64):                                                  # B: resolution
        for kA in (0.25, 0.32):
            for itag in integrators("m5"):
                C.append(case("B", "m5", 1.6, kA, n, itag))
    for kA in (0.25, 0.32):                                             # C: variable bed
        for itag in integrators("m6"):
            C.append(case("C", "m6", 1.6, kA, 32, itag))
    for g in (0.5, 1.0):
        for itag in integrators("m6"):
            C.append(case("C", "m6", 1.6, 0.25, 32, itag, gamma=g))
    for spp in (20, 80):                                                # D: time step
        for itag in integrators("m5"):
            C.append(case("D", "m5", 1.6, 0.32, 32, itag, spp=spp))
    for T in (1.6, 3.0):                                                # E: null control
        for kA in (0.158, 0.25, 0.32):
            for itag in integrators("m3"):
                C.append(case("E", "m3", T, kA, 32, itag))
    for kA in (0.158, 0.25, 0.32):                                      # F: m3 selective damping
        for itag, spp in (("crz34", 40), ("crz34", 20), ("crz34", 10), ("rk4", 20), ("rk4", 10)):
            C.append(case("F", "m3", 1.6, kA, 32, itag, spp=spp))
    for kA in (0.158, 0.25, 0.32):                                      # G: γ map
        for g in GAMMAS_G:
            for itag in integrators("m5"):
                C.append(case("G", "m5", 1.6, kA, 32, itag, gamma=g))
    for g in (0.5, 1.0):                                                # H: γ × dt
        for spp in (20, 80):
            for itag in ("sdirk22", "crz34"):
                C.append(case("H", "m5", 1.6, 0.32, 32, itag, spp=spp, gamma=g))
    names = [c["name"] for c in C]
    assert len(names) == len(set(names)), "duplicate case names"
    return C


def packs(C):
    """Within each (tier, model): the cost-sorted cases cut into the fewest packs of at most
    PACK_SIZE, as EQUAL in size as possible (10 cases → 5 + 5, not 8 + 2: same bill, since each
    pack rents a 1/8 node, but the packs end together). (pack id, tier, case names)."""
    out, pid = [], 0
    for tier in TIERS:
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
            "solver", "tableau", "theta", "rho", "gamma", "save_every", "ckpt_every", "cost"]
    with open(os.path.join(HERE, "cases.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for c in C:
            f.write("\t".join(str(c[k]) for k in cols) + "\n")
    P = packs(C)
    with open(os.path.join(HERE, "packs.tsv"), "w") as f:
        f.write("pack\ttier\tcases\n")
        for pid, tier, names in P:
            f.write(f"{pid}\t{tier}\t{','.join(names)}\n")
    for tier in TIERS:
        nc = sum(c["tier"] == tier for c in C)
        npk = sum(p[1] == tier for p in P)
        cost = sum(c["cost"] for c in C if c["tier"] == tier)
        print(f"tier {tier}: {nc:3d} cases, {npk:2d} packs, relative cost {cost:7.1f}")
    print(f"total: {len(C)} cases, {len(P)} packs")
