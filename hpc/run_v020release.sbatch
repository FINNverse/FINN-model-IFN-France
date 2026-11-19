#!/bin/bash
# ===========================================================================
# FINN-model-IFN-France | run ONE v020-release script on a SLURM CPU core
#   sbatch --time=24:00:00 --mem=16G hpc/run_v020release.sbatch calibration/v020-release/fit.R
# A copy of hpc/run_v020soil.sbatch for calibration/v020-release/ (the v020-soil
# fit refitted unchanged on the released FINN 0.2.0, commit 6c9156f; data
# fiveyear-v3-soil). Resources only (one CPU, one torch thread; time and memory
# on the sbatch line). Every scientific choice is inside the R script. Expects
# the repository at $HOME/finn-zoo/FINN-zoo-models/FINN-model-IFN-France and the
# data at $HOME/finn-zoo/FINN-zoo-data/frnfi (the ../../FINN-zoo-data/frnfi
# layout the scripts use), the released FINN 0.2.0 (6c9156f) in
# $HOME/R/lib-finn-release-6c9156f (FINN only; first on the library path here
# and in the scripts' .libPaths()), torch 0.17 in $HOME/R/site-library-torch017
# and the remaining packages in $HOME/R/site-library.
# ===========================================================================
#SBATCH --job-name=zoo-fr-v020rel
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/v020rel_%j.out
#SBATCH --error=logs/v020rel_%j.err
set -euo pipefail
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
# the release library (FINN only), then torch 0.17, then the general site library;
# torch 0.17's libtorch on LD_LIBRARY_PATH: the libtorch of the older torch in the
# site library fails with "Lantern is not loaded" (the France finding of 2026-11-10).
export R_LIBS_USER="$HOME/R/lib-finn-release-6c9156f:$HOME/R/site-library-torch017:$HOME/R/site-library"
export LD_LIBRARY_PATH="$HOME/R/site-library-torch017/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd "${SLURM_SUBMIT_DIR:-$HOME/finn-zoo/FINN-zoo-models/FINN-model-IFN-France}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch hpc/run_v020release.sbatch calibration/v020-release/<script>.R}"
echo "[v020rel] $SCRIPT on $(hostname), $OMP_NUM_THREADS thread(s), $(date)"
Rscript "$SCRIPT"
echo "[v020rel] done $(date)"
