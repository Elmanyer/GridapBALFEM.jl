# ============================================================================
#  cluster.sh — the per-cluster part of the integrator study (sourced by submit.sh and by the
#  self-resubmission in run_pack.sh, so a pack's later segments land where the first one did)
#
#  run_pack.sh's #SBATCH header is the Snellius set-up (rome 1/8 node: 16 cores, 28 GiB, 119:59).
#  Flags on the sbatch command line override the header, so another cluster only adds them here.
#  BALFEM_CLUSTER (snellius | blue, default snellius) also selects the modules in run/balfem_env.sh.
#    blue: DelftBlue compute-p2 (64 cores, 4 GB per core, 5-day wall like Snellius). A pack's
#          8 single-threaded cases need 8 cores and 28 GiB. Account: BALFEM_ACCOUNT (default
#          innovation). Output on /scratch: the home directory is too small for ≈ 30 GB.
# ============================================================================
case "${BALFEM_CLUSTER:-snellius}" in
    snellius)
        STUDY_SBATCH_FLAGS=() ;;
    blue)
        STUDY_SBATCH_FLAGS=(--partition=compute-p2 --account="${BALFEM_ACCOUNT:-innovation}" --cpus-per-task=8)
        : "${BALFEM_STUDY_OUT:=/scratch/$USER/integrator_study}" ;;
    *)
        echo "unknown BALFEM_CLUSTER=${BALFEM_CLUSTER} (snellius | blue)" >&2
        return 1 ;;
esac
