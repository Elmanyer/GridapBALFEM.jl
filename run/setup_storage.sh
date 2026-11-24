#!/bin/bash
# ============================================================================
#  setup_storage.sh — ONE-TIME move of the cluster output/ to the project storage (run/storage.sh)
#
#      export STORAGE=/gpfs/admin/projectlinks/prjs1960      # in ~/.bashrc (Snellius, prjs1960)
#      bash $HOME/GridapBALFEM.jl/run/setup_storage.sh
#
#  1. copies $BALFEM_PROJ/output/ into $BALFEM_OUTPUT_ROOT (= $STORAGE/GridapBALFEM_output) with
#     rsync -a (resumable: run it again if it is interrupted);
#  2. renames the home copy to output.home_backup_<date> (nothing is deleted);
#  3. creates the link $BALFEM_PROJ/output -> $BALFEM_OUTPUT_ROOT.
#  Afterwards check the copy (du -sh on both, a few files), then delete the backup YOURSELF to
#  free the home quota:  rm -rf $BALFEM_PROJ/output.home_backup_<date>
#  Do not run it while a job is writing to output/ (squeue -u $USER must be empty).
# ============================================================================
set -euo pipefail
: "${BALFEM_PROJ:=$HOME/GridapBALFEM.jl}"
source "$BALFEM_PROJ/run/storage.sh"
balfem_storage_root

src="$BALFEM_PROJ/output"
if [ -L "$src" ]; then
    echo "already set up: $src -> $(readlink "$src")"; exit 0
fi
if squeue -h -u "$USER" 2>/dev/null | grep -q .; then
    echo "jobs of $USER are queued or running (squeue): wait until none writes to output/" >&2
    exit 1
fi
mkdir -p "$BALFEM_OUTPUT_ROOT"
if [ -d "$src" ]; then
    echo "copying $src/ → $BALFEM_OUTPUT_ROOT/ ($(du -sh "$src" | cut -f1)) …"
    rsync -a "$src/" "$BALFEM_OUTPUT_ROOT/"
    bak="$src.home_backup_$(date +%F)"
    mv "$src" "$bak"
    echo "home copy kept as $bak (delete it yourself once the copy is checked)"
fi
ln -s "$BALFEM_OUTPUT_ROOT" "$src"
echo "linked $src -> $BALFEM_OUTPUT_ROOT"
