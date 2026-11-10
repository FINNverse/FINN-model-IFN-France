#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | run ONE v0.2 pilot script on a SLURM CPU core
#   sbatch hpc/run_pilot.sbatch calibration/v02/pilot_B_thetafix.R
# Resources only (one CPU, one torch thread, 8 GB, 2 h). Every scientific
# choice is inside the pilot script. Expects the repository copy and the data
# as listed in calibration/v02/CLUSTER.md, FINN 0.2.0 (0848b60) in
# $HOME/R/lib-finn-0848b60 and the other packages in $R_LIBS_USER.
# ===========================================================================
#SBATCH --job-name=ifnfr-pilot
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/pilot_%j.out
#SBATCH --error=logs/pilot_%j.err
set -euo pipefail
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/site-library}"
# libtorch of the torch that the pilots load (lib-finn-0848b60 is first on .libPaths());
# the libtorch of an older torch in $R_LIBS_USER fails with "Lantern is not loaded".
export LD_LIBRARY_PATH="$HOME/R/lib-finn-0848b60/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$HOME/FINN-zoo/frnfi}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch hpc/run_pilot.sbatch calibration/v02/pilot_<name>.R}"
echo "[pilot] $SCRIPT $(date)"
Rscript "$SCRIPT"
echo "[pilot] done $(date)"
