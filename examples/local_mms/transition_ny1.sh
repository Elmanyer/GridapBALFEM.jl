#!/bin/bash
#  Waits for the ny_1d=3 Phase-B batch to drain, then switches every :d1 study to
#  ny_1d=1 and re-runs them. Measured equivalence beforehand (static, Q3/Q2, nx=16):
#      ny=3  ndofs=2817  e_eta=3.365575179792e-05  e_u=1.302487479913e-06
#      ny=1  ndofs= 957  e_eta=3.365575179792e-05  e_u=1.302487479915e-06
#  i.e. agreement to 14 and 12 significant figures -- LU round-off, not a difference.
#  mms.jl says the same analytically: the :d1 field has ky=0, so u*_y == 0 and the
#  solution is y-invariant -- "refine nx only; refining ny changes nothing".
#
#  The 6 :d2 studies are NOT affected by ny_1d and are the expensive ones, so their
#  rows are PRESERVED and they are not recomputed.
cd "$(dirname "$0")/../.." || exit 1
OUT=output/local/mms_phaseB

# 1. wait for the supervisor to exit (it breaks on COMPLETE or on a 3-round stall)
while pgrep -f 'supervise_phaseB2[.]sh' >/dev/null; do sleep 60; done
while pgrep -f 'run_phaseB_shard[.]jl [0-9]' >/dev/null; do sleep 60; done
echo "[trans] batch drained at $(date +%H:%M:%S)"

# 2. archive everything as it stands
ARCH=$OUT/ny3_archive; mkdir -p $ARCH
cp $OUT/shard_*.csv $ARCH/ 2>/dev/null
echo "[trans] archived $(ls $ARCH/*.csv 2>/dev/null | wc -l) csv files to $ARCH"

# 3. keep ONLY the :d2 rows, so the :d1 studies are seen as not-done and re-run
for f in $OUT/shard_*.csv; do
  [ -f "$f" ] || continue
  head -1 "$f" > "$f.tmp"
  #  FNR, not NR: with one file per shard, NR>1 would skip only the FIRST header
  #  and leak the other ten into the filtered output.
  awk -F, 'FNR>1 && $10=="d2"' "$f" >> "$f.tmp"
  mv "$f.tmp" "$f"
done
echo "[trans] kept $(cat $OUT/shard_*.csv | grep -v '^task,' | grep -c . ) :d2 rows; :d1 rows cleared for re-run"

# 4. flip the default
python3 - <<'PY'
import pathlib
p=pathlib.Path("src/mms_driver.jl"); s=p.read_text()
a="levels::Int = 4, nx0::Int = 8, ny0::Int = 8, ny_1d::Int = 3,"
b="""levels::Int = 4, nx0::Int = 8, ny0::Int = 8,
                          #  ny_1d = 1, NOT 3. The :d1 manufactured field sets ky=0, so
                          #  u*_y == 0 identically and the solution is y-invariant --
                          #  mms.jl: "refine nx only; refining ny changes nothing".
                          #  Measured: ny=1 and ny=3 agree to 14 (e_eta) and 12 (e_u)
                          #  significant figures at 2.94x fewer DOFs, and far more than
                          #  that in wall time, because a 3-cell-wide Q3 mesh widens the
                          #  LU front. ny>=3 is a PERIODIC-direction requirement only;
                          #  run_mms_case builds y_wall_bc=:wall, so ny=1 is legal here.
                          ny_1d::Int = 1,"""
assert s.count(a)==1, "ny_1d anchor"
p.write_text(s.replace(a,b))
print("[trans] src/mms_driver.jl: ny_1d 3 -> 1")
PY

# 5. re-arm and restart
rm -f $OUT/exhausted_*
setsid nohup ./examples/local_mms/supervise_phaseB2.sh >> $OUT/supervisor2.log 2>&1 </dev/null &
echo "[trans] supervisor restarted for the ny=1 re-run at $(date +%H:%M:%S)"
