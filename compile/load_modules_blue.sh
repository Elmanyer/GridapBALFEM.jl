
# DelftBlue (2026 software stack). OpenMPI from the modules; Julia from juliaup, NOT the
# julia/1.11.6 module: Manifest.toml is resolved for Julia 1.10.10, and 1.11 would re-resolve it
# (other package versions than Snellius and the development machine run). One-time install:
#     curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel 1.10.10
module unload julia 2>/dev/null || true
module load 2026 cpu openmpi/5.0.10 slurm/current

export PATH="$HOME/.juliaup/bin:$HOME/.julia/bin:$PATH"      # juliaup's julia, then mpiexecjl
julia --version 2>/dev/null | grep -q '^julia version 1\.10\.' || {
    echo "[load_modules_blue] need Julia 1.10 from juliaup (found: $(julia --version 2>&1 | head -1))." >&2
    echo "[load_modules_blue] install: curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel 1.10.10" >&2
    return 1
}

export MPI_LIB_DIR=$(mpicc --showme:libdirs | awk '{print $1}')   # read by compile/set_preferences.jl
export LD_LIBRARY_PATH=$MPI_LIB_DIR:${LD_LIBRARY_PATH:-}
