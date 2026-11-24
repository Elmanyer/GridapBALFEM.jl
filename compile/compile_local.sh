#!/usr/bin/env bash
# ==============================================================
#  compile_local.sh — build GridapBALFEM_sysimage_local.so on this workstation (see compile_local.jl)
#
#  1. dry run of warmup.jl (seconds: every case's keywords against its driver);
#  2. the build, one MPI rank via mpiexecjl (hours: the cold compile of the complete solver);
#  3. the src/ stamp <image>.src.sha256, same hash as run/balfem_env.sh balfem_src_hash.
#  Run it detached, in its own scope, so an OOM kill elsewhere cannot take it down:
#      systemd-run --user --scope --unit=balfem-sysimage -p MemoryMax=16G \
#          nohup compile/compile_local.sh > output/sysimage_build/build.log 2>&1 &
# ==============================================================
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJ"
IMAGE="$PROJ/GridapBALFEM_sysimage_local.so"

echo "=== $(date '+%F %T')  warmup dry run"
BALFEM_WARMUP_DRYRUN=1 julia --project="$PROJ" "$PROJ/compile/warmup.jl"

echo "=== $(date '+%F %T')  build (src/ at $(git -C "$PROJ" rev-parse --short HEAD)$(git -C "$PROJ" diff --quiet -- src || echo ' + uncommitted src edits'))"
"$HOME/.julia/bin/mpiexecjl" --project="$PROJ" -n 1 julia --project="$PROJ" "$PROJ/compile/compile_local.jl"

( cd "$PROJ/src" && find . -name '*.jl' -type f -print0 \
    | LC_ALL=C sort -z | xargs -0 -r sha256sum ) | sha256sum | awk '{print $1}' > "$IMAGE.src.sha256"
echo "=== $(date '+%F %T')  built: $IMAGE  ($(du -h "$IMAGE" | cut -f1)); src hash $(cat "$IMAGE.src.sha256")"
