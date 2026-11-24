# ============================================================================
#  storage.sh — where cluster runs write their output (sourced; never executed)
#
#  WHY. On 2026-11-21 the integrator study filled the Snellius HOME quota with VTK output and
#  161 of its 275 cases died writing a file (markdown_files/INTEGRATOR_STUDY_RESULTS.md §0). Home
#  is for code and the sysimage only. Every cluster output goes to the project storage.
#
#  THE LAYOUT. One output root on the storage space, with the same layout as a local output/:
#      $STORAGE/GridapBALFEM_output/            (= BALFEM_OUTPUT_ROOT)
#          integrator_study/<case>/seg<NN>/     slurm_logs/<campaign>/     <run dirs of the drivers>
#  and the repository's output/ on the cluster is a SYMLINK to it:
#      $BALFEM_PROJ/output -> $BALFEM_OUTPUT_ROOT
#  so every driver that writes to <repo>/output/… (the default of every example and test, and
#  BALFEM_OUTDIR values built from $BALFEM_PROJ/output) lands on the storage space without
#  touching a single Julia file.
#
#  ONE-TIME SET-UP per cluster account (Snellius project space of the author: prjs1960):
#      echo 'export STORAGE=/gpfs/admin/projectlinks/prjs1960' >> ~/.bashrc && source ~/.bashrc
#      bash $HOME/GridapBALFEM.jl/run/setup_storage.sh      # moves an existing output/ there
#
#  KNOBS
#      STORAGE              the storage root. Required on Snellius (no default: a silent fall-back
#                           to $HOME is exactly the failure this file exists to prevent).
#                           DelftBlue default: /scratch/$USER.
#      BALFEM_OUTPUT_ROOT   default $STORAGE/GridapBALFEM_output
#
#  `balfem_storage_check` (called by balfem_env.sh, submit.sh and the other submit scripts, so it
#  fails at submission, before the queue): resolves the root, creates it, creates the symlink when
#  output/ does not exist, and REFUSES when output/ is a real directory in home (run
#  setup_storage.sh) or a link pointing elsewhere.
# ============================================================================

balfem_storage_root() {
    case "${BALFEM_CLUSTER:-snellius}" in
        blue) : "${STORAGE:=/scratch/$USER}" ;;
    esac
    if [ -z "${STORAGE:-}" ]; then
        echo "[storage] STORAGE is not set. Cluster output must not go to \$HOME (quota)." >&2
        echo "[storage]   Snellius: echo 'export STORAGE=/gpfs/admin/projectlinks/prjs1960' >> ~/.bashrc && source ~/.bashrc" >&2
        return 1
    fi
    [ -d "$STORAGE" ] || { echo "[storage] STORAGE=$STORAGE does not exist" >&2; return 1; }
    : "${BALFEM_OUTPUT_ROOT:=$STORAGE/GridapBALFEM_output}"
    export STORAGE BALFEM_OUTPUT_ROOT
}

balfem_storage_check() {
    balfem_storage_root || return 1
    mkdir -p "$BALFEM_OUTPUT_ROOT" || return 1
    local link="${BALFEM_PROJ:?BALFEM_PROJ not set}/output"
    if [ -L "$link" ]; then
        if [ "$(readlink -f "$link")" != "$(readlink -f "$BALFEM_OUTPUT_ROOT")" ]; then
            echo "[storage] $link -> $(readlink "$link"), expected $BALFEM_OUTPUT_ROOT" >&2
            echo "[storage]   fix the link (or set BALFEM_OUTPUT_ROOT to its target)" >&2
            return 1
        fi
    elif [ -e "$link" ]; then
        echo "[storage] $link is a real directory in home; cluster output must go to the storage." >&2
        echo "[storage]   run once: bash $BALFEM_PROJ/run/setup_storage.sh" >&2
        return 1
    else
        ln -s "$BALFEM_OUTPUT_ROOT" "$link" || return 1
        echo "[storage] linked $link -> $BALFEM_OUTPUT_ROOT"
    fi
    return 0
}
