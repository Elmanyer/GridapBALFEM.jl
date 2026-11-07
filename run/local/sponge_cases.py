#!/usr/bin/env python3
"""sponge_cases.py — case table of the sponge-reflection sweep (CH8_CAMPAIGN_TODO.md §7).

    sponge_cases.py A                     stage A: 80 matrix cases + 4 floor controls
    sponge_cases.py B <W> <mu>            stage B: relaxation width 12 m at the chosen (W, mu_max)

Writes TSV to stdout, one case per line, longest first (LPT, rule 43):
    name  T  W  mu  relax_w  Lx  nx  t_final  x_fit0  x_fit1  role

Geometry: relaxation zone relax_w, then a FIXED 32 m interior, then the outflow sponge W, so
Lx = relax_w + 32 + W and the fit stretch [relax_w + 1, relax_w + 31] is the same 30 m in every
case. dx = 0.25 m (the Snellius suite).
Duration (rule 14): T_ramp (2T) + the round trip relax edge → sponge FAR end → back,
2 (32 + W) / c_g, with c_g from linear (Airy) theory at d = 3.5 m, + 10 analysis periods + 5 T
margin, rounded up to a whole number of save intervals.
"""
import math, sys

G, D, DX, DT, SAVE_EVERY = 9.81, 3.5, 0.25, 0.02, 5
INTERIOR = 32.0
PERIODS = [1.5, 2.0, 2.5, 2.9]
MUS = [4, 8, 16, 25, 40]
WIDTHS = [6, 12, 18, 24]

def airy_k(T):
    w2 = (2 * math.pi / T) ** 2 / G
    k = max(w2, math.sqrt(w2 / D))
    for _ in range(100):
        f = k * math.tanh(k * D) - w2
        df = math.tanh(k * D) + k * D / math.cosh(k * D) ** 2
        k -= f / df
    return k

def c_group(T):
    k = airy_k(T)
    c = (2 * math.pi / T) / k
    return 0.5 * c * (1 + 2 * k * D / math.sinh(2 * k * D))

def case(T, W, mu, relax_w, role):
    Lx = relax_w + INTERIOR + W
    nx = int(round(Lx / DX))
    t = 2 * T + 2 * (INTERIOR + W) / c_group(T) + 15 * T
    q = DT * SAVE_EVERY
    t = math.ceil(t / q) * q
    tag = "spA" if role in ("matrix", "floor", "replica") else "spB"
    name = f"{tag}_T{T:g}_W{W:g}_mu{mu:g}" + (f"_rw{relax_w:g}" if relax_w != 6 else "")
    return dict(name=name, T=T, W=W, mu=mu, relax_w=relax_w, Lx=Lx, nx=nx, t_final=round(t, 4),
                x_fit0=relax_w + 1.0, x_fit1=relax_w + 31.0, role=role)

def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "A"
    cases = []
    if stage == "A":
        for T in PERIODS:
            for W in WIDTHS:
                for mu in MUS:
                    role = "replica" if (W == 12 and mu == 40) else "matrix"
                    cases.append(case(T, W, mu, 6.0, role))
            cases.append(case(T, 48, 4, 6.0, "floor"))
    elif stage == "B":
        W, mu = float(sys.argv[2]), float(sys.argv[3])
        for T in PERIODS:
            cases.append(case(T, W, mu, 12.0, "stageB"))
    else:
        sys.exit("usage: sponge_cases.py A | B <W> <mu>")
    cases.sort(key=lambda c: -(c["t_final"] / DT) * c["nx"])
    cols = ["name", "T", "W", "mu", "relax_w", "Lx", "nx", "t_final", "x_fit0", "x_fit1", "role"]
    print("\t".join(cols))
    for c in cases:
        print("\t".join(str(c[k]) for k in cols))

main()
