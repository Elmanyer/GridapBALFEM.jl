#!/usr/bin/env python3
"""make_pvd.py — rebuild a ParaView .pvd index from the snapshots on disk.

WHY THIS EXISTS. The solver writes its collection through WriteVTK's
`createpvd(...) do pvd ... end`, and that serialises the .pvd only when the block
EXITS. A run that is still going -- or one that was killed, or diverged, or whose
job hit the wall clock -- therefore leaves every snapshot on disk with no index,
and ParaView cannot open the series. This script reconstructs the index from the
filenames, so a running simulation can be inspected without waiting for it, and a
crashed one can be inspected at all.

The timestep is recovered from the filename: the solver formats it as
`@sprintf("%.4f", t)` with the decimal point replaced by an underscore to avoid a
WriteVTK extension warning, so `sol_t_22_4000.pvtu` is t = 22.4000.

USAGE
    python3 postprocessing/make_pvd.py <run-dir> [-o solution.pvd] [--watch SEC]

`--watch` regenerates every SEC seconds, which is what you want against a live
run: ParaView re-reads the .pvd on "Reload Files".
"""
import argparse, glob, os, re, sys, time

def timestep_of(path):
    stem = os.path.basename(path).rsplit(".", 1)[0]      # sol_t_22_4000
    m = re.match(r"^sol_t_(.+)$", stem)
    if not m:
        return None
    tok = m.group(1)                                      # 22_4000
    if "_" not in tok:
        return None
    whole, frac = tok.rsplit("_", 1)
    try:
        return float(f"{whole}.{frac}")
    except ValueError:
        return None

def build(run_dir, out_name):
    # .pvtu is the distributed index; .vtu the sequential single-piece file.
    # Prefer .pvtu when both exist -- the .vtu files then live inside the piece dirs.
    files = sorted(glob.glob(os.path.join(run_dir, "sol_t_*.pvtu")))
    if not files:
        files = [f for f in sorted(glob.glob(os.path.join(run_dir, "sol_t_*.vtu")))]
    entries = []
    for f in files:
        t = timestep_of(f)
        if t is not None:
            entries.append((t, os.path.basename(f)))
    entries.sort(key=lambda e: e[0])
    if not entries:
        return 0, None, None
    out = os.path.join(run_dir, out_name)
    with open(out, "w") as fh:
        #  header matched to what WriteVTK emits for a solver-written collection,
        #  so ParaView cannot treat a reconstructed index differently from a native
        #  one. The compressor attribute is inert on a Collection (a .pvd holds no
        #  data) but is kept for byte-level similarity.
        fh.write('<?xml version="1.0" encoding="utf-8"?>\n')
        fh.write('<VTKFile type="Collection" version="1.0" byte_order="LittleEndian"'
                 ' compressor="vtkZLibDataCompressor">\n')
        fh.write('  <Collection>\n')
        for t, name in entries:
            fh.write(f'    <DataSet timestep="{t:.4f}" part="0" file="{name}"/>\n')
        fh.write('  </Collection>\n</VTKFile>\n')
    return len(entries), entries[0][0], entries[-1][0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--out", default="solution.pvd")
    ap.add_argument("--watch", type=float, default=0.0,
                    help="regenerate every N seconds (for a live run)")
    a = ap.parse_args()
    if not os.path.isdir(a.run_dir):
        sys.exit(f"not a directory: {a.run_dir}")
    while True:
        n, t0, t1 = build(a.run_dir, a.out)
        if n == 0:
            sys.exit(f"no sol_t_*.pvtu or sol_t_*.vtu found in {a.run_dir}")
        print(f"{a.out}: {n} snapshots, t = {t0:.4f} .. {t1:.4f}", flush=True)
        if a.watch <= 0:
            break
        time.sleep(a.watch)

if __name__ == "__main__":
    main()
