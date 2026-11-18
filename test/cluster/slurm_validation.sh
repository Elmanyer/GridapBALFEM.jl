#!/bin/bash
# ==============================================================
#  slurm_validation.sh — SLURM template for the cluster VALIDATION runs
#  (mirror of ../../run/run_snellius.sh; edit for your cluster).
#
#  Submit one validation at a time by setting VALIDATION=, or duplicate the
#  final block. px*py MUST equal the total task count (-n).
# ==============================================================
#SBATCH --job-name="BALFEM-validate"
#SBATCH --partition=rome
#SBATCH --time=24:00:00
#SBATCH --nodes=1                     # increase for multi-node; keep ntasks = nodes*ntasks-per-node
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --output=BALFEM-validate.%j.out
#SBATCH --error=BALFEM-validate.%j.err

# --- environment (edit paths/modules for your machine) --------------------------
ROOT=$HOME/GridapBABALFEM.jl
source $ROOT/compile/load_modules_snellius.sh    # or your cluster's module setup

NP=$SLURM_NTASKS                                  # total ranks (= -n)
export BABALFEM_PX=${BABALFEM_PX:-$NP}                    # process grid; PX*PY must equal NP
export BABALFEM_PY=${BABALFEM_PY:-1}

# --- which validation to run ----------------------------------------------------
VALIDATION=${VALIDATION:-conservation}            # conservation | mms

case "$VALIDATION" in
  conservation)
    # long-run mass/energy conservation on a big closed basin
    export BABALFEM_NX=${BABALFEM_NX:-2000}  BABALFEM_NY=${BABALFEM_NY:-200}
    export BABALFEM_LX=${BABALFEM_LX:-200.0} BABALFEM_LY=${BABALFEM_LY:-20.0}
    export BABALFEM_DT=${BABALFEM_DT:-0.02}  BABALFEM_TFINAL=${BABALFEM_TFINAL:-40.0}
    export BABALFEM_A=${BABALFEM_A:-0.01}    BABALFEM_ADVECTION=${BABALFEM_ADVECTION:-1}
    export BABALFEM_OUTDIR=${BABALFEM_OUTDIR:-$ROOT/output/cluster_conservation}
    SCRIPT=$ROOT/test/cluster/cluster_conservation.jl ;;
  mms)
    # unsteady nonlinear manufactured-solution recovery at scale (all 𝓝 comps)
    export BABALFEM_NX=${BABALFEM_NX:-400}   BABALFEM_NY=${BABALFEM_NY:-200}
    export BABALFEM_NSTEPS=${BABALFEM_NSTEPS:-200} BABALFEM_DT=${BABALFEM_DT:-0.05}
    export BABALFEM_NLPFULL=${BABALFEM_NLPFULL:-1}  BABALFEM_AMP=${BABALFEM_AMP:-1.0}
    export BABALFEM_OUTDIR=${BABALFEM_OUTDIR:-$ROOT/output/cluster_selfconsistency}
    SCRIPT=$ROOT/test/cluster/cluster_selfconsistency.jl ;;
  *) echo "unknown VALIDATION=$VALIDATION (use conservation|mms)"; exit 1 ;;
esac

echo "### BALFEM validation: $VALIDATION | -n $NP | grid ${BABALFEM_PX}x${BABALFEM_PY}"
mpiexecjl -n $NP julia --project=$ROOT $SCRIPT
